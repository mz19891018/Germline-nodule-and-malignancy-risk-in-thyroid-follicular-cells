# -*- coding: utf-8 -*-
"""
Step 01  GWAS 汇总统计整理为 MAGMA 输入格式（SNP  P  N）
- 若文件含 rsID：直接使用 rsID，对应原版 1000G 欧洲参考
- 若无 rsID：一次性把 1000G 参考面板坐标 hg19 -> hg38 转换，用“染色体:位置:等位基因”匹配
"""
import os, shutil
from pathlib import Path
import numpy as np
import pandas as pd
import config as C
from utils import ensure, get_logger, find_col, trait_file_map, all_traits, load_json, save_json

log = get_logger("step01_gwas_prep")
OUT = ensure(C.WORK / "magma_input")
REF = ensure(C.WORK / "magma_ref")

COLS = dict(
    rsid=["rsid", "rs_id", "hm_rsid", "variant_id", "snp", "markername", "snpid", "id"],
    chr=["chromosome", "hm_chrom", "chr", "chrom", "#chrom"],
    pos=["base_pair_location", "hm_pos", "pos", "bp", "position", "genpos"],
    ea=["effect_allele", "hm_effect_allele", "ea", "alt", "a1", "allele1"],
    oa=["other_allele", "hm_other_allele", "nea", "ref", "a2", "allele0", "allele2"],
    p=["p_value", "pvalue", "p", "pval", "p-value", "p_val"],
    n=["n", "n_total", "sample_size", "neff", "n_eff"],
    eaf=["effect_allele_frequency", "hm_effect_allele_frequency", "eaf", "af", "freq", "a1freq"],
)


def detect_columns(path):
    hdr = pd.read_csv(path, sep="\t", nrows=2000)
    cols = list(hdr.columns)
    m = {k: find_col(cols, v, contains=False) for k, v in COLS.items()}
    if m["p"] is None:
        m["p"] = find_col(cols, ["p_value", "pval", "p"], contains=True)
    use_rsid = False
    if m["rsid"] is not None:
        frac = hdr[m["rsid"]].astype(str).str.startswith("rs").mean()
        use_rsid = frac > 0.5
    return m, use_rsid, cols


def build_hg38_reference(g1000_prefix):
    """1000G 欧洲参考面板：生成以 chr:pos(hg38):A:B 为ID的bim（坐标仍为hg19用于注释）。"""
    new_prefix = REF / "g1000_eur_hg38id"
    if Path(str(new_prefix) + ".bim").exists():
        return str(new_prefix)
    from pyliftover import LiftOver
    log.info("首次运行：把1000G参考面板坐标 hg19 -> hg38（约5-15分钟，仅需一次；需联网下载chain文件）")
    chain = list(C.MAGMA_DIR.glob("hg19ToHg38.over.chain*"))
    lo = LiftOver(str(chain[0])) if chain else LiftOver("hg19", "hg38")   # 无本地chain文件则自动联网下载
    bim = pd.read_csv(g1000_prefix + ".bim", sep=r"\s+", header=None,
                      names=["chr", "snp", "cm", "pos", "a1", "a2"], dtype={"chr": str, "snp": str})
    ids = np.empty(len(bim), dtype=object)
    for i, (c, p, a1, a2) in enumerate(zip(bim.chr.values, bim.pos.values, bim.a1.values, bim.a2.values)):
        r = lo.convert_coordinate("chr" + c, int(p) - 1)
        if r and r[0][0] == "chr" + c:
            x, y = sorted([str(a1).upper(), str(a2).upper()])
            ids[i] = f"{c}:{r[0][1] + 1}:{x}:{y}"
        else:
            ids[i] = f"unmapped_{i}"
        if i % 2_000_000 == 0:
            log.info(f"  已转换 {i:,} / {len(bim):,}")
    s = pd.Series(ids)
    dup = s.duplicated(keep="first")      # 保留第一个，其余重复ID加后缀
    s[dup] = s[dup] + "_dup" + s[dup].groupby(s[dup]).cumcount().astype(str)
    bim["snp"] = s.values
    bim.to_csv(str(new_prefix) + ".bim", sep="\t", header=False, index=False)
    for ext in [".bed", ".fam"]:
        dst = Path(str(new_prefix) + ext)
        if not dst.exists():
            try:
                os.link(g1000_prefix + ext, dst)
            except Exception:
                shutil.copy(g1000_prefix + ext, dst)
    # 注释用的 SNP 位置文件（ID + hg19 坐标）
    bim[["snp", "chr", "pos"]].to_csv(REF / "g1000_eur_hg38id.snploc", sep="\t", header=False, index=False)
    log.info(f"参考面板转换完成，成功映射 {(~s.str.startswith('unmapped')).mean():.1%}")
    return str(new_prefix)


def norm_chr(x):
    return x.astype(str).str.replace("chr", "", case=False).str.replace("^23$", "X", regex=True)


def prepare_trait(key, rec, mode_info):
    out = OUT / f"{key}.txt"
    if out.exists() and key in mode_info:
        log.info(f"{key}: 已存在，跳过")
        return
    final_out = out
    out = OUT / f"{key}.tmp"
    path = rec.file
    m, use_rsid, cols = detect_columns(path)
    log.info(f"{key}: 文件 {Path(path).name}；列映射 {m}；使用rsID={use_rsid}")
    need = [m["p"]] + ([m["rsid"]] if use_rsid else [m["chr"], m["pos"], m["ea"], m["oa"]])
    if None in need:
        raise SystemExit(f"{key}: 无法识别必要的列，已有列为 {cols}")
    usecols = [c for c in set(need + [m["n"], m["eaf"]]) if c]
    n_const = None
    if m["n"] is None:
        n_const = rec.N
        if pd.isna(n_const):
            raise SystemExit(f"{key}: 文件无样本量列，且GWAS Catalog未返回样本量，请在 gwas_map_manual.csv 的 N 列填写")
    total, kept = 0, 0
    first = True
    for ch in pd.read_csv(path, sep="\t", usecols=usecols, chunksize=2_000_000, dtype=str):
        total += len(ch)
        p = pd.to_numeric(ch[m["p"]], errors="coerce")
        keep = p.notna() & (p > 0) & (p <= 1)
        if m["eaf"]:
            f = pd.to_numeric(ch[m["eaf"]], errors="coerce")
            maf = np.minimum(f, 1 - f)
            keep &= (maf >= C.MIN_MAF) | f.isna()
        ch = ch[keep]
        p = p[keep].clip(lower=1e-300)
        if use_rsid:
            snp = ch[m["rsid"]].astype(str)
            ok = snp.str.startswith("rs")
        else:
            a = ch[m["ea"]].astype(str).str.upper()
            b = ch[m["oa"]].astype(str).str.upper()
            ok = a.str.fullmatch("[ACGT]") & b.str.fullmatch("[ACGT]") & (a != b)
            lo = np.where(a < b, a, b)
            hi = np.where(a < b, b, a)
            snp = norm_chr(ch[m["chr"]]) + ":" + ch[m["pos"]].astype(str) + ":" + lo + ":" + hi
        n = pd.to_numeric(ch[m["n"]], errors="coerce") if m["n"] else pd.Series(n_const, index=ch.index)
        df = pd.DataFrame({"SNP": snp, "P": p, "N": n.round().astype("Int64")})[ok.values]
        df = df.dropna()
        kept += len(df)
        df.to_csv(out, sep="\t", index=False, header=first, mode="w" if first else "a")
        first = False
        log.info(f"  {key}: 已处理 {total:,} 行，保留 {kept:,}")
    os.replace(out, final_out)
    mode_info[key] = {"mode": "rsid" if use_rsid else "chrpos", "n_total": total, "n_kept": kept,
                      "file": Path(path).name, "accession": rec.accession}


def main():
    paths = load_json(C.RES / "00_check" / "magma_paths.json")
    fmap = trait_file_map()
    info_path = C.RES / "00_check" / "gwas_prep_info.json"
    mode_info = load_json(info_path) if info_path.exists() else {}
    for key in all_traits():
        if key not in fmap:
            log.warning(f"{key}: 未找到对应GWAS文件，跳过（如为主分析性状请补充 gwas_map_manual.csv）")
            continue
        prepare_trait(key, fmap[key], mode_info)
        save_json(mode_info, info_path)
    if any(v["mode"] == "chrpos" for v in mode_info.values()):
        mode_info["_ref_hg38id"] = build_hg38_reference(paths["g1000_prefix"])
    save_json(mode_info, info_path)
    log.info("Step01 完成")


if __name__ == "__main__":
    main()
