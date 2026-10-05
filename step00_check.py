# -*- coding: utf-8 -*-
"""
Step 00  数据清点与准备
1. 检查 Python 包
2. 解压 MAGMA、单细胞、cBioPortal 压缩包
3. 读取40个GWAS文件的表头；联网查询 GWAS Catalog，自动识别每个编号对应的疾病和人群
4. 输出 results/00_check/ 下的清点报告
"""
import re, json, importlib
from pathlib import Path
import pandas as pd
import config as C
from utils import ensure, get_logger, head_lines, extract_archive, save_json

log = get_logger("step00_check")
OUT = ensure(C.RES / "00_check")

PKGS = ["numpy", "pandas", "scipy", "scanpy", "anndata", "harmonypy", "leidenalg",
        "scdrs", "pyliftover", "statsmodels", "lifelines", "requests", "matplotlib",
        "seaborn", "gseapy", "mygene"]


def check_packages():
    missing = []
    for p in PKGS:
        try:
            m = importlib.import_module(p)
            log.info(f"包 {p:12s} OK  版本 {getattr(m, '__version__', '?')}")
        except Exception:
            missing.append(p)
            log.warning(f"包 {p} 未安装")
    if missing:
        raise SystemExit("缺少Python包: " + ", ".join(missing) + "。请按 README 安装环境后重试。")


def prepare_archives():
    w = ensure(C.WORK)
    for z in C.MAGMA_DIR.glob("*.zip"):
        extract_archive(z, w / "magma_files" / z.stem, log)
    for t in C.SC_DIR.glob("*_RAW.tar"):
        extract_archive(t, w / t.name.replace(".tar", ""), log)
    for t in C.CLIN_DIR.glob("*.tar.gz"):
        extract_archive(t, w / "cbioportal", log)

    exe = list((w / "magma_files").rglob("magma.exe")) + list((w / "magma_files").rglob("magma"))
    exe = [p for p in exe if p.is_file()]
    bed = list((w / "magma_files").rglob("g1000_eur.bed"))
    gl = list((w / "magma_files").rglob("NCBI37.3.gene.loc"))
    paths = {"magma_exe": str(exe[0]) if exe else None,
             "g1000_prefix": str(bed[0].with_suffix("")) if bed else None,
             "gene_loc": str(gl[0]) if gl else None}
    for k, v in paths.items():
        log.info(f"{k}: {v}")
        if v is None:
            raise SystemExit(f"未找到 {k}，请检查 02_magma 中的压缩包是否完整")
    save_json(paths, OUT / "magma_paths.json")


def _walk_strings(obj, keys_like, found):
    if isinstance(obj, dict):
        for k, v in obj.items():
            if any(s in k.lower() for s in keys_like):
                found.append(json.dumps(v, ensure_ascii=False) if not isinstance(v, str) else v)
            _walk_strings(v, keys_like, found)
    elif isinstance(obj, list):
        for v in obj:
            _walk_strings(v, keys_like, found)


def query_catalog(acc):
    import requests
    for url in [f"https://www.ebi.ac.uk/gwas/rest/api/v2/studies/{acc}",
                f"https://www.ebi.ac.uk/gwas/rest/api/studies/{acc}"]:
        try:
            r = requests.get(url, timeout=40)
            if r.ok:
                return r.json()
        except Exception as e:
            log.warning(f"{acc} 查询失败 {url}: {e}")
    return None


def classify_trait(t):
    t = t.lower()
    is_cancer = ("cancer" in t) or ("carcinoma" in t) or ("malignan" in t)
    is_goiter = ("goit" in t) or ("benign" in t) or ("nodul" in t)
    if is_cancer and is_goiter:
        return "ThCvsBNG"
    if "graves" in t:
        return "GD"
    if "thyroiditis" in t or "hashimoto" in t:
        return "LT"
    if "hypothyroid" in t:
        return "HYPO"
    if is_goiter:
        return "BNG"
    if is_cancer:
        return "ThC"
    return None


def classify_ancestry(t):
    t = t.lower()
    groups = {"EUR": ["european"], "AFR": ["african"], "EAS": ["east asian"],
              "AMR": ["hispanic", "latin", "admixed american"], "SAS": ["south asian"]}
    hit = [g for g, kws in groups.items() if any(k in t for k in kws)]
    if len(hit) == 1:
        return hit[0]
    if len(hit) > 1 or "multi" in t or "mixed" in t:
        return "MIX"
    return None


def parse_n(texts):
    best = None
    for t in texts:
        nums = [int(x.replace(",", "")) for x in re.findall(r"\d[\d,]{2,}", t)]
        if nums:
            s = sum(nums)
            best = s if best is None else max(best, s)
    return best


def inventory_gwas():
    inv = OUT / "gwas_inventory.csv"
    need = C.MAIN_TRAITS + C.SENS_TRAITS + C.CTRL_TRAITS
    if (OUT / "gwas_map_manual.csv").exists():
        log.info("检测到手动对应表 gwas_map_manual.csv，跳过自动识别")
        return
    if inv.exists() and set(need) <= set(pd.read_csv(inv)["key"].dropna()):
        log.info("GWAS 识别结果已存在，跳过")
        return
    rows = []
    files = sorted(C.GWAS_DIR.glob("GCST*.tsv*"))
    log.info(f"发现 {len(files)} 个 GWAS 文件")
    for f in files:
        acc = re.match(r"(GCST\d+)", f.name).group(1)
        hdr = head_lines(f, 3)
        cols = hdr[0].split("\t") if hdr else []
        js = query_catalog(acc)
        trait_txt, anc_txt, n_txt = [], [], []
        if js is not None:
            _walk_strings(js, ["trait"], trait_txt)
            _walk_strings(js, ["ancestr", "sample"], anc_txt)
            _walk_strings(js, ["sample_size", "samplesize"], n_txt)
            (OUT / "catalog_json").mkdir(exist_ok=True)
            with open(OUT / "catalog_json" / f"{acc}.json", "w", encoding="utf-8") as fh:
                json.dump(js, fh, ensure_ascii=False, indent=1)
        trait = " | ".join(dict.fromkeys(trait_txt))[:300]
        anc = " | ".join(dict.fromkeys(anc_txt))[:500]
        dis = classify_trait(trait) if trait else None
        an = classify_ancestry(anc) if anc else None
        key = f"{dis}_{an}" if dis and an else ""
        N = parse_n(n_txt or anc_txt)
        rows.append(dict(accession=acc, file=str(f), size_MB=round(f.stat().st_size / 1e6),
                         trait_catalog=trait, ancestry_catalog=anc, disease=dis, ancestry=an,
                         key=key, N=N, columns=",".join(cols),
                         example_row=hdr[1] if len(hdr) > 1 else ""))
        log.info(f"{acc}: {key or '未识别'} | N={N} | {trait[:80]}")
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "gwas_inventory.csv", index=False, encoding="utf-8-sig")

    need = C.MAIN_TRAITS + C.SENS_TRAITS + C.CTRL_TRAITS
    have = set(df["key"])
    miss = [k for k in need if k not in have]
    if miss:
        tmpl = df[["accession", "file", "trait_catalog", "ancestry_catalog", "N", "key"]]
        tmpl.to_csv(OUT / "gwas_map_manual_TEMPLATE.csv", index=False, encoding="utf-8-sig")
        log.warning("以下性状未自动识别: " + ", ".join(miss))
        log.warning("请打开 results/00_check/gwas_map_manual_TEMPLATE.csv，在 key 列手动填写"
                    "（如 ThCvsBNG_EUR），另存为 gwas_map_manual.csv 后重新运行。")
        if not (OUT / "gwas_map_manual.csv").exists():
            raise SystemExit("GWAS 性状识别不完整，需要手动补充（见上方提示）。")


def inventory_files():
    lines = []
    for d in [C.WORK / "GSE184362_RAW", C.WORK / "GSE193581_RAW", C.WORK / "cbioportal",
              C.TCGA_DIR, C.CLIN_DIR, C.SC_DIR]:
        lines.append(f"\n===== {d} =====")
        if d.exists():
            fs = sorted([p for p in d.rglob("*") if p.is_file()])
            lines.append(f"共 {len(fs)} 个文件")
            for p in fs[:80]:
                lines.append(f"{p.relative_to(d)}\t{p.stat().st_size // 1024} KB")
    # 关键表格文件的表头
    lines.append("\n===== 表头预览 =====")
    cands = list((C.WORK / "cbioportal").rglob("data_clinical_*.txt")) + \
        list(C.CLIN_DIR.glob("Survival_*")) + list(C.SC_DIR.glob("*annotation*"))
    for p in cands:
        lines.append(f"\n--- {p.name} ---")
        try:
            lines += [l[:400] for l in head_lines(p, 8)]
        except Exception as e:
            lines.append(f"读取失败: {e}")
    (OUT / "file_inventory.txt").write_text("\n".join(lines), encoding="utf-8")
    log.info("文件清单已写入 results/00_check/file_inventory.txt")


def main():
    check_packages()
    prepare_archives()
    inventory_files()
    inventory_gwas()
    log.info("Step00 完成")


if __name__ == "__main__":
    main()
