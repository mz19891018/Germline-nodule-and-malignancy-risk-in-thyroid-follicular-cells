# -*- coding: utf-8 -*-
"""
Step 06  独立单细胞数据集 GSE193581 验证
- 用 GSE184362 得到的 GIMP / BNP 基因程序，在新数据的滤泡细胞中打分
- 在新数据中重新计算 scDRS 恶性特异残差
- 检验：GIMP 分数是否与恶性特异残差正相关（细胞水平 + 样本水平）
"""
import re
import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
from scipy import stats
import config as C
from utils import ensure, get_logger, save_json
import sc_io
from step04_scdrs import run_scdrs

log = get_logger("step06_replication")
OUT = ensure(C.RES / "06_replication")
RAW = C.WORK / "GSE193581_RAW"
H5 = C.WORK / "thy_GSE193581.h5ad"


def load_annotation(a):
    f = list(C.SC_DIR.glob("GSE193581_celltype_annotation*"))
    if not f:
        return None
    ann = pd.read_csv(f[0], sep=None, engine="python")
    ann.to_csv(OUT / "annotation_preview.csv", index=False) if len(ann) < 50 else ann.head(30).to_csv(OUT / "annotation_preview.csv", index=False)
    raw_bc = pd.Series(a.obs_names.str.split("_", n=1).str[1].str.replace(r"-\d+$", "", regex=True), index=a.obs_names)
    best, best_col, best_key = 0, None, None
    cand_cols = [None] + list(ann.columns)
    for col in cand_cols:
        vals = ann.index.astype(str) if col is None else ann[col].astype(str)
        v2 = vals.str.replace(r"-\d+$", "", regex=True)
        for keyname, ref in [("full", pd.Series(a.obs_names, index=a.obs_names)), ("barcode", raw_bc)]:
            ov = ref.isin(set(v2) | set(vals)).mean()
            if ov > best:
                best, best_col, best_key = ov, col, keyname
    log.info(f"注释文件匹配率 {best:.1%}（列={best_col}，方式={best_key}）")
    if best < 0.2:
        return None
    vals = (ann.index.astype(str) if best_col is None else ann[best_col].astype(str)).str.replace(r"-\d+$", "", regex=True)
    ann = ann.set_index(vals)
    ann = ann[~ann.index.duplicated()]
    ref = pd.Series(a.obs_names, index=a.obs_names) if best_key == "full" else raw_bc
    label_col = [c for c in ann.columns if re.search("type|cluster|annot|label", c, re.I) and c != best_col]
    if not label_col:
        return None
    return ref.map(ann[label_col[0]]).astype(str)


def build():
    if H5.exists():
        return sc.read_h5ad(H5)
    groups = sc_io.discover_samples(RAW)
    rows, adatas = [], []
    for s, fs in groups.items():
        a = sc_io.read_sample(fs)
        if a is None:
            continue
        name = " ".join(f.name for f in fs)
        tis = sc_io.guess_tissue(name)
        a.obs_names = [f"{s}_{x}" for x in a.obs_names]
        a.obs["sample"], a.obs["tissue"] = s, tis
        a.obs["patient"] = sc_io.guess_patient(name.replace(s, ""))
        rows.append(dict(sample=s, files=name[:300], tissue_guess=tis, n_cells=a.n_obs))
        adatas.append(a)
    pd.DataFrame(rows).to_csv(OUT / "sample_metadata_auto.csv", index=False, encoding="utf-8-sig")
    a = ad.concat(adatas, join="outer", fill_value=0)
    a.obs_names_make_unique()
    a = sc_io.qc_filter(a, log)
    a = sc_io.normalize(a)
    lab = load_annotation(a)
    a = sc_io.embed_cluster(a, C.LEIDEN_RES_ALL)
    sc_io.annotate_by_markers(a)
    if lab is not None:
        a.obs["author_label"] = lab.reindex(a.obs_names).values
        is_thy = a.obs.author_label.str.contains("thyro|epith|tumou?r|malig|follic|cancer", case=False, regex=True)
        log.info("使用作者注释识别滤泡细胞")
    else:
        is_thy = a.obs.cell_type == "Thyrocyte"
        log.info("使用标志基因注释识别滤泡细胞")
    pd.crosstab(a.obs.cell_type, a.obs.tissue).to_csv(OUT / "celltype_by_tissue.csv")
    t = a[is_thy.values & (a.obs.tissue != "ATC").values].copy()
    log.info(f"验证集滤泡细胞 {t.n_obs} 个（已排除ATC样本）")
    t.write_h5ad(H5, compression="gzip")
    return t


def main():
    log.info("Step06 代码版本 v3")
    t = build()
    # 按样本名重新判定组织来源（修正旧版把 NORM 样本误判为 Tumor 的问题）
    mf = OUT / "sample_metadata_auto.csv"
    if not mf.exists():
        groups = sc_io.discover_samples(RAW)
        pd.DataFrame([dict(sample=s, files=" ".join(f.name for f in fs)) for s, fs in groups.items()]) \
            .to_csv(mf, index=False, encoding="utf-8-sig")
    meta = pd.read_csv(mf, encoding="utf-8-sig")
    tis = meta.set_index("sample")["files"].map(sc_io.guess_tissue)
    t.obs["tissue"] = t.obs["sample"].map(tis).astype(str)
    t.obs["patient"] = t.obs["sample"].map(meta.set_index("sample")["files"].str.extract(r"(?:PTC|NORM|ATC)(\d+)", expand=False)).astype(str)
    log.info("验证集组织来源：" + t.obs.groupby("tissue")["sample"].nunique().to_string())
    t = t[t.obs.tissue != "ATC"].copy()
    gimp = pd.read_csv(C.RES / "05_program" / "GIMP_genes.csv").gene
    bnp = pd.read_csv(C.RES / "05_program" / "BNP_genes.csv").gene
    t = sc_io.stratified_subsample(t, "sample", C.MAX_CELLS_THY, C.SEED)
    grp, scores = run_scdrs(t, [C.MAL, C.NOD], ["tissue"], "rep")
    df = t.obs[["sample", "tissue", "patient"]].copy()
    df["S_mal"] = scores[C.MAL]["norm_score"].reindex(df.index)
    df["S_nod"] = scores[C.NOD]["norm_score"].reindex(df.index)
    X = np.c_[np.ones(len(df)), df.S_nod]
    df["resid_mal"] = df.S_mal - X @ np.linalg.lstsq(X, df.S_mal, rcond=None)[0]
    X = np.c_[np.ones(len(df)), df.S_mal]
    df["resid_nod"] = df.S_nod - X @ np.linalg.lstsq(X, df.S_nod, rcond=None)[0]
    tt = t.copy()
    tt.X = tt.layers["counts"].copy(); sc.pp.normalize_total(tt, target_sum=1e4); sc.pp.log1p(tt)
    for nm, gs in [("GIMP", gimp), ("BNP", bnp)]:
        g = [x for x in gs if x in tt.var_names]
        sc.tl.score_genes(tt, g, score_name=nm, random_state=C.SEED)
        df[nm] = tt.obs[nm].reindex(df.index)
        log.info(f"{nm}: 验证集中可用 {len(g)}/{len(gs)} 个基因")
    df.to_csv(OUT / "replication_cell_scores.csv.gz")
    res = {}
    for prog, target in [("GIMP", "resid_mal"), ("GIMP", "resid_nod"), ("BNP", "resid_nod"), ("BNP", "resid_mal")]:
        r, p = stats.spearmanr(df[prog], df[target])
        per = []
        for s, d in df.groupby("sample"):
            if len(d) >= 50:
                per.append(stats.spearmanr(d[prog], d[target])[0])
        per = np.array(per)
        p_sign = stats.wilcoxon(per).pvalue if len(per) >= 3 else np.nan
        res[f"{prog}~{target}"] = dict(cell_rho=r, cell_p=p, n_samples=len(per),
                                       median_sample_rho=float(np.median(per)) if len(per) else None,
                                       frac_samples_positive=float((per > 0).mean()) if len(per) else None,
                                       p_across_samples=p_sign)
        log.info(f"{prog} vs {target}: 细胞 rho={r:.3f}；样本中位 rho={res[f'{prog}~{target}']['median_sample_rho']}")
    for c in ["sample", "patient", "tissue"]:
        df[c] = df[c].astype(str)
    pb = df.groupby(["sample", "patient", "tissue"], observed=True)[["GIMP", "BNP", "resid_mal", "resid_nod"]] \
        .mean().dropna(how="all").reset_index()
    pb.to_csv(OUT / "replication_pseudobulk.csv", index=False)
    rows = []
    for v in ["GIMP", "BNP", "resid_mal", "resid_nod"]:
        a, b = pb[pb.tissue == "Tumor"][v].dropna(), pb[pb.tissue == "Paratumor"][v].dropna()
        if len(a) >= 3 and len(b) >= 3:
            rows.append(dict(variable=v, n_tumor_samples=len(a), n_normal_samples=len(b),
                             mean_tumor=a.mean(), mean_normal=b.mean(),
                             p_mannwhitney=stats.mannwhitneyu(a.dropna(), b.dropna()).pvalue))
        w = pb.pivot_table(index="patient", columns="tissue", values=v).dropna()
        if {"Tumor", "Paratumor"} <= set(w.columns) and len(w) >= 1:
            rows.append(dict(variable=v, paired_patients=",".join(w.index.astype(str)),
                             paired_diffs=",".join(f"{x:.3f}" for x in (w["Tumor"] - w["Paratumor"]))))
    pd.DataFrame(rows).to_csv(OUT / "replication_tissue_tests.csv", index=False)
    save_json(res, OUT / "replication_summary.json")
    log.info("Step06 完成")


if __name__ == "__main__":
    main()
