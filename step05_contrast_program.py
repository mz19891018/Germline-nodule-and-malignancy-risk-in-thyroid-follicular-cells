# -*- coding: utf-8 -*-
"""
Step 05  核心分析（v2：全部采用“样本内/患者内”设计，排除患者间差异的混杂）
1. 恶性特异风险残差：ThCvsBNG 的 scDRS 分数对 BNG 分数回归取残差（resid_mal）
2. 组织来源：同一患者 Tumor / LN / DistantMet vs Paratumor 配对比较 + 细胞水平混合模型
3. 分化程度：resid_mal ~ TDS，样本固定效应 + 按患者聚类稳健标准误
4. 细胞状态：样本内比较（描述性；肿瘤细胞状态多为患者特异）
5. GIMP 基因程序：样本内 Spearman（在每个样本内部排秩后合并），并要求跨样本方向一致；
   排除 GWAS 输入基因以避免循环论证
6. 敏感性分析：用混合人群 GWAS（MIX）重复 2、3
"""
import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
from scipy import stats
import statsmodels.formula.api as smf
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
import config as CFG   # 不能命名为C：会与公式里的 C() 冲突
from utils import ensure, get_logger, bh_fdr, save_json

log = get_logger("step05_contrast_program")
IN = CFG.RES / "04_scdrs"
OUT = ensure(CFG.RES / "05_program")
MIN_CELLS = 20
MIN_CELLS_GENE = 100
ORDER = ["Paratumor", "Tumor", "LN", "DistantMet"]


def residual(y, x):
    X = np.c_[np.ones(len(x)), x]
    b = np.linalg.lstsq(X, y, rcond=None)[0]
    return y - X @ b, b[1]


def load_scores():
    meta = pd.read_csv(IN / "thy_scored_cells_meta.csv.gz", index_col=0)
    obs = ad.read_h5ad(CFG.WORK / "thy_GSE184362.h5ad", backed="r").obs
    meta["TDS"] = obs["TDS"].reindex(meta.index).values
    for tr, nm in [(CFG.MAL, "S_mal"), (CFG.NOD, "S_nod"), ("ThC_EUR", "S_thc"),
                   ("ThCvsBNG_MIX", "S_mal_mix"), ("BNG_MIX", "S_nod_mix")]:
        f = IN / f"thy_cellscore_{tr}.csv.gz"
        if f.exists():
            meta[nm] = pd.read_csv(f, index_col=0)["norm_score"].reindex(meta.index)
    meta = meta.dropna(subset=["S_mal", "S_nod"])
    meta["resid_mal"], b = residual(meta.S_mal.values, meta.S_nod.values)
    meta["resid_nod"], _ = residual(meta.S_nod.values, meta.S_mal.values)
    if "S_mal_mix" in meta and "S_nod_mix" in meta:
        ok = meta[["S_mal_mix", "S_nod_mix"]].notna().all(1)
        meta.loc[ok, "resid_mal_mix"] = residual(meta.loc[ok, "S_mal_mix"].values, meta.loc[ok, "S_nod_mix"].values)[0]
    r, p = stats.spearmanr(meta.S_mal, meta.S_nod)
    save_json({"cell_rho_mal_nod": r, "cell_p": p, "beta_mal_on_nod": b, "n_cells": len(meta),
               "n_samples": meta["sample"].nunique(), "n_patients": meta["patient"].nunique(),
               "cells_by_tissue": meta.tissue.value_counts().to_dict()},
              OUT / "residual_model.json")
    log.info(f"细胞水平 S_mal 与 S_nod：rho={r:.3f}；回归系数={b:.3f}；细胞数={len(meta)}")
    return meta


# ---------------- 2. 组织来源 ----------------
def tissue_tests(meta, var):
    pb = meta.groupby(["patient", "tissue"])[var].agg(["mean", "size"]).reset_index()
    pb = pb[pb["size"] >= MIN_CELLS]
    pb.to_csv(OUT / f"pseudobulk_{var}_by_patient_tissue.csv", index=False)
    wide = pb.pivot(index="patient", columns="tissue", values="mean")
    rows = []
    if "Paratumor" in wide.columns:
        for t in [x for x in wide.columns if x != "Paratumor"]:
            d = wide[[t, "Paratumor"]].dropna()
            if len(d) >= 3:
                diff = d[t] - d["Paratumor"]
                rows.append(dict(variable=var, comparison=f"{t} vs Paratumor (paired patients)",
                                 n_patients=len(d), mean_diff=diff.mean(),
                                 frac_patients_positive=(diff > 0).mean(),
                                 p_wilcoxon=stats.wilcoxon(diff).pvalue,
                                 p_paired_t=stats.ttest_rel(d[t], d["Paratumor"]).pvalue))
    d = meta[["patient", "tissue", var]].dropna()
    if "Paratumor" in set(d.tissue):
        try:
            m = smf.mixedlm(f"{var} ~ C(tissue, Treatment('Paratumor'))", d, groups=d["patient"]).fit(reml=True)
            for k in m.params.index:
                if "tissue" in k:
                    rows.append(dict(variable=var, comparison=k.split("T.")[-1].rstrip("]") + " vs Paratumor (mixed model)",
                                     n_patients=d.patient.nunique(), mean_diff=m.params[k], p_mixed=m.pvalues[k]))
        except Exception as e:
            log.warning(f"混合模型失败：{e}")
    return pd.DataFrame(rows)


# ---------------- 3. 分化程度 ----------------
def differentiation_test(meta, var):
    d = meta[["sample", "patient", "TDS", var]].dropna()
    d = d.assign(TDS_z=(d.TDS - d.TDS.mean()) / d.TDS.std())
    m = smf.ols(f"{var} ~ TDS_z + C(sample)", d).fit(
        cov_type="cluster", cov_kwds={"groups": d["patient"].astype("category").cat.codes})
    per = [stats.spearmanr(g.TDS, g[var])[0] for s, g in d.groupby("sample") if len(g) >= MIN_CELLS_GENE]
    per = np.array(per)
    return dict(variable=var, n_cells=len(d), n_samples=d["sample"].nunique(),
                beta_per_SD_TDS=m.params["TDS_z"], se_cluster=m.bse["TDS_z"], p_cluster=m.pvalues["TDS_z"],
                n_samples_rho=len(per), median_within_sample_rho=float(np.median(per)) if len(per) else np.nan,
                frac_samples_negative=float((per < 0).mean()) if len(per) else np.nan,
                p_sign_across_samples=stats.wilcoxon(per).pvalue if len(per) >= 5 else np.nan)


# ---------------- 4. 细胞状态（描述性） ----------------
def state_tests(meta, var):
    rows = []
    stt = meta.thy_state.astype(str)
    for st in sorted(stt.unique()):
        diffs = []
        for s, d in meta.groupby("sample"):
            ds = d.thy_state.astype(str)
            a, b = d[ds == st][var], d[ds != st][var]
            if len(a) >= MIN_CELLS and len(b) >= MIN_CELLS:
                diffs.append(a.mean() - b.mean())
        diffs = np.array(diffs)
        sub = meta[stt == st]
        dom = sub["patient"].value_counts(normalize=True)
        rows.append(dict(thy_state=st, variable=var, n_cells=len(sub),
                         n_patients_present=int((sub["patient"].value_counts() >= MIN_CELLS).sum()),
                         top_patient=dom.index[0] if len(dom) else None,
                         top_patient_frac=dom.iloc[0] if len(dom) else np.nan,
                         main_tissue=sub.tissue.value_counts().index[0] if len(sub) else None,
                         mean_TDS=sub.TDS.mean(),
                         n_samples_tested=len(diffs), mean_within_sample_diff=diffs.mean() if len(diffs) else np.nan,
                         p_wilcoxon=stats.wilcoxon(diffs).pvalue if len(diffs) >= 5 else np.nan))
    df = pd.DataFrame(rows)
    df["fdr"] = bh_fdr(df.p_wilcoxon.values)
    return df


# ---------------- 5. 基因程序 ----------------
def gwas_input_genes():
    gs = pd.read_csv(CFG.RES / "02_magma" / "scdrs_genesets.gs", sep="\t")
    genes = set()
    for tr in (CFG.MAL, CFG.NOD, "ThC_EUR", "ThCvsBNG_MIX", "BNG_MIX"):
        row = gs[gs.TRAIT == tr]
        if len(row):
            genes |= {x.split(":")[0] for x in row.GENESET.iloc[0].split(",")}
    return genes


# 组织解离诱导的即刻早期/热休克应激基因（参考 van den Brink et al., Nat Methods 2017 核心基因），
# 从基因程序候选中剔除，避免技术性应激信号混入
DISSOCIATION_GENES = {
    "FOS", "FOSB", "JUN", "JUNB", "JUND", "ATF3", "EGR1", "EGR2", "EGR3", "IER2", "IER3", "IER5",
    "DUSP1", "NR4A1", "NR4A2", "NR4A3", "BTG2", "CYR61", "CTGF", "KLF2", "KLF4", "KLF6",
    "ZFP36", "ZFP36L1", "ZFP36L2", "PPP1R15A", "GADD45B", "GADD45G", "SOCS3", "MYC", "RHOB",
    "NFKBIA", "TNFAIP3", "IRF1", "HSPA1A", "HSPA1B", "HSPA8", "HSPB1", "HSP90AA1", "HSPH1",
    "DNAJB1", "HSPE1", "HSPD1", "CEBPB", "CEBPD", "SRF", "MAFF", "SAT1", "TPPP3", "UBC"}


def _zrank(M):
    R = np.apply_along_axis(stats.rankdata, 0, M)
    return (R - R.mean(0)) / (R.std(0) + 1e-12)


def load_thy_expr(cells):
    t = sc.read_h5ad(CFG.WORK / "thy_GSE184362.h5ad")
    t = t[pd.Index(cells).intersection(t.obs_names)].copy()
    t.X = t.layers["counts"].copy()
    sc.pp.normalize_total(t, target_sum=1e4); sc.pp.log1p(t)
    return t


def within_sample_gene_corr(meta, score_col, max_cells=30000):
    cells = meta.index
    if len(cells) > max_cells:
        cells = pd.Index(np.random.default_rng(CFG.SEED).choice(cells, max_cells, replace=False))
    t = load_thy_expr(cells)
    frac = np.asarray((t.X > 0).mean(0)).ravel()
    excl = gwas_input_genes()
    keep = (frac >= CFG.MIN_EXPR_FRAC) & ~t.var_names.str.match(r"^(RP[SL]\d|MT-|HB[AB])") \
        & ~t.var_names.isin(list(excl)) & ~t.var_names.isin(list(DISSOCIATION_GENES))
    t = t[:, keep].copy()
    samp = meta.loc[t.obs_names, "sample"].values
    y_all = meta.loc[t.obs_names, score_col].values
    X = t.X.tocsc()
    G = t.n_vars
    num = np.zeros(G); ntot = 0; per_rho = []
    for s in pd.unique(samp):
        idx = np.where(samp == s)[0]
        if len(idx) < MIN_CELLS_GENE:
            continue
        y = _zrank(y_all[idx][:, None])[:, 0]
        rho_s = np.zeros(G)
        for i in range(0, G, 500):
            blk = X[idx, i:i + 500].toarray()
            rho_s[i:i + 500] = _zrank(blk).T @ y / len(idx)
        rho_s = np.nan_to_num(rho_s)
        num += rho_s * len(idx); ntot += len(idx); per_rho.append(rho_s)
    P = np.vstack(per_rho)
    rho = num / ntot
    n_eff = ntot - len(per_rho)
    tval = rho * np.sqrt((n_eff - 2) / (1 - rho ** 2 + 1e-12))
    df = pd.DataFrame({"gene": t.var_names, "within_sample_rho": rho, "expr_frac": frac[keep],
                       "n_samples": P.shape[0], "frac_samples_positive": (P > 0).mean(0),
                       "median_sample_rho": np.median(P, 0)})
    df["p"] = 2 * stats.t.sf(np.abs(tval), n_eff - 2)
    df["fdr"] = bh_fdr(df.p.values)
    log.info(f"{score_col}: {len(per_rho)} 个样本、{ntot} 个细胞、{G} 个基因（已排除 {len(excl)} 个GWAS输入基因）")
    return df.sort_values("within_sample_rho", ascending=False)


def select_program(df):
    sel = df[(df.within_sample_rho > 0) & (df.fdr < 0.05) & (df.frac_samples_positive >= 0.6)]
    return sel.head(CFG.PROGRAM_SIZE)


def enrich(genes, name):
    try:
        import gseapy as gp
        r = gp.enrichr(gene_list=list(genes), gene_sets=["MSigDB_Hallmark_2020", "GO_Biological_Process_2023",
                                                          "KEGG_2021_Human"], organism="human", outdir=None)
        r.results.to_csv(OUT / f"enrichment_{name}.csv", index=False)
    except Exception as e:
        log.warning(f"富集分析失败（需要联网）：{e}")


def program_scores(meta, gimp, bnp):
    t = load_thy_expr(meta.index)
    for nm, gs in [("GIMP", gimp), ("BNP", bnp)]:
        g = [x for x in gs if x in t.var_names]
        if len(g) >= 5:
            sc.tl.score_genes(t, g, score_name=nm, random_state=CFG.SEED)
            meta[nm] = t.obs[nm].reindex(meta.index)
    return meta


def main():
    log.info(f"Step05 代码版本 v3（剔除 {len(DISSOCIATION_GENES)} 个组织解离应激基因）")
    meta = load_scores()

    tv = [v for v in ["resid_mal", "S_mal", "S_nod", "resid_mal_mix"] if v in meta]
    pd.concat([tissue_tests(meta, v) for v in tv]).to_csv(OUT / "tissue_tests.csv", index=False)

    dv = [v for v in ["resid_mal", "S_mal", "S_nod", "resid_nod", "resid_mal_mix"] if v in meta]
    diff = pd.DataFrame([differentiation_test(meta, v) for v in dv])
    diff.to_csv(OUT / "differentiation_tests.csv", index=False)
    log.info("\n" + diff.to_string())

    st = state_tests(meta, "resid_mal")
    grp = pd.read_csv(IN / "thy_group_results.csv")
    g = grp[grp.group_col == "thy_state"].pivot(index="group", columns="trait", values="assoc_mcp")
    g.columns = [f"scDRS_mcp_{c}" for c in g.columns]
    st = st.merge(g, left_on="thy_state", right_index=True, how="left")
    st["malignancy_specific_state"] = (st.mean_within_sample_diff > 0) & (st.fdr < 0.05)
    st.to_csv(OUT / "malignancy_specific_states.csv", index=False)

    gm = within_sample_gene_corr(meta, "resid_mal"); gm.to_csv(OUT / "gene_corr_resid_mal.csv", index=False)
    gimp = select_program(gm); gimp.to_csv(OUT / "GIMP_genes.csv", index=False)
    gn = within_sample_gene_corr(meta, "resid_nod"); gn.to_csv(OUT / "gene_corr_resid_nod.csv", index=False)
    bnp = select_program(gn); bnp.to_csv(OUT / "BNP_genes.csv", index=False)
    log.info(f"GIMP {len(gimp)} 个基因；BNP {len(bnp)} 个基因；重叠 {len(set(gimp.gene) & set(bnp.gene))}")
    enrich(gimp.gene, "GIMP"); enrich(bnp.gene, "BNP")

    meta = program_scores(meta, gimp.gene.tolist(), bnp.gene.tolist())
    meta.to_csv(OUT / "thy_cell_scores_with_residual.csv.gz")
    pv = [v for v in ["GIMP", "BNP"] if v in meta]
    pd.concat([tissue_tests(meta, v) for v in pv]).to_csv(OUT / "program_tissue_tests.csv", index=False)
    pd.DataFrame([differentiation_test(meta, v) for v in pv]).to_csv(OUT / "program_differentiation_tests.csv", index=False)

    fig, ax = plt.subplots(1, 3, figsize=(14, 4))
    pb = meta.groupby(["patient", "tissue"])["resid_mal"].agg(["mean", "size"]).reset_index()
    pb = pb[pb["size"] >= MIN_CELLS]
    order = [x for x in ORDER if x in set(pb.tissue)]
    sns.boxplot(data=pb, x="tissue", y="mean", ax=ax[0], color="lightsteelblue", showfliers=False, order=order)
    sns.stripplot(data=pb, x="tissue", y="mean", ax=ax[0], color="k", size=4, order=order)
    ax[0].set_ylabel("Malignancy-specific risk (patient mean)")
    d = meta.sample(min(20000, len(meta)), random_state=0)
    d = d.assign(resid_c=d.resid_mal - d.groupby("sample").resid_mal.transform("mean"),
                 TDS_c=d.TDS - d.groupby("sample").TDS.transform("mean"))
    ax[1].hexbin(d.TDS_c, d.resid_c, gridsize=40, cmap="Blues", mincnt=1)
    ax[1].set_xlabel("TDS (within-sample centred)"); ax[1].set_ylabel("Malignancy-specific risk (centred)")
    top = gimp.head(20)
    ax[2].barh(top.gene[::-1], top.within_sample_rho[::-1], color="crimson")
    ax[2].set_xlabel("Within-sample Spearman rho (GIMP top 20)")
    fig.tight_layout(); fig.savefig(OUT / "fig_core_results.png", dpi=300); plt.close(fig)
    log.info("Step05 完成")


if __name__ == "__main__":
    main()
