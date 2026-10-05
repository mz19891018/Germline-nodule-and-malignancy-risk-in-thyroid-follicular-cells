# -*- coding: utf-8 -*-
"""
Step 04  scDRS：把各 GWAS 的基因水平风险投射到单个细胞
A. 全图谱（各细胞类型）：全部性状，含自身免疫对照（预期落在免疫细胞）
B. 甲状腺滤泡细胞（各细胞状态）：主分析 + 敏感性分析性状
"""
import numpy as np
import pandas as pd
import scanpy as sc
import scdrs
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
import config as C
from utils import ensure, get_logger, all_traits, bh_fdr
import sc_io

log = get_logger("step04_scdrs")
OUT = ensure(C.RES / "04_scdrs")
GS = C.RES / "02_magma" / "scdrs_genesets.gs"


def run_scdrs(a, traits, group_cols, tag, gene_traits=()):
    a = a.copy()
    if "counts" in a.layers:          # 确保 X 为 log 标准化表达
        a.X = a.layers["counts"].copy()
        sc.pp.normalize_total(a, target_sum=1e4); sc.pp.log1p(a)
    sc.pp.filter_genes(a, min_cells=max(10, int(0.001 * a.n_obs)))
    n_genes = np.asarray((a.X > 0).sum(1)).ravel()
    cov = pd.DataFrame({"const": 1.0, "n_genes": n_genes}, index=a.obs_names)
    scdrs.preprocess(a, cov=cov, n_mean_bin=20, n_var_bin=20, copy=False)
    gs = scdrs.util.load_gs(str(GS), to_intersect=a.var_names)
    group_rows, cell_scores = [], {}
    for tr in traits:
        if tr not in gs:
            log.warning(f"{tr} 不在基因集中，跳过")
            continue
        genes, w = gs[tr]
        log.info(f"[{tag}] scDRS 打分 {tr}：{len(genes)} 个基因，{a.n_obs} 个细胞")
        df = scdrs.score_cell(a, genes, gene_weight=w, ctrl_match_key="mean_var",
                              n_ctrl=C.SCDRS_N_CTRL, weight_opt="vs", return_ctrl_raw_score=False,
                              return_ctrl_norm_score=True, random_seed=C.SEED, verbose=False)
        keep = ["raw_score", "norm_score", "mc_pval", "pval", "nlog10_pval", "zscore"]
        cell_scores[tr] = df[[c for c in keep if c in df.columns]]
        cell_scores[tr].to_csv(OUT / f"{tag}_cellscore_{tr}.csv.gz")
        res = scdrs.method.downstream_group_analysis(adata=a, df_full_score=df, group_cols=group_cols)
        for gcol, g in res.items():
            g = g.copy(); g["trait"] = tr; g["group_col"] = gcol; g.index.name = "group"
            group_rows.append(g.reset_index())
        if tr in gene_traits:
            ga = scdrs.method.downstream_gene_analysis(a, df)
            ga.to_csv(OUT / f"{tag}_gene_corr_{tr}.csv")
        del df
    grp = pd.concat(group_rows, ignore_index=True)
    grp["assoc_fdr"] = np.nan
    for (tr, gc), idx in grp.groupby(["trait", "group_col"]).groups.items():
        grp.loc[idx, "assoc_fdr"] = bh_fdr(grp.loc[idx, "assoc_mcp"].values)
    grp.to_csv(OUT / f"{tag}_group_results.csv", index=False)
    return grp, cell_scores


def heatmap(grp, group_col, tag):
    g = grp[grp.group_col == group_col]
    mat = g.pivot(index="group", columns="trait", values="assoc_mcp")
    mat = -np.log10(mat.clip(lower=1e-3))
    sig = g.pivot(index="group", columns="trait", values="assoc_fdr") < 0.05
    fig, ax = plt.subplots(figsize=(1 + 0.7 * mat.shape[1], 1 + 0.35 * mat.shape[0]))
    sns.heatmap(mat, cmap="Reds", annot=sig.replace({True: "*", False: ""}), fmt="",
                cbar_kws={"label": "-log10 MC P"}, ax=ax)
    ax.set_title(f"scDRS association ({tag})")
    fig.tight_layout(); fig.savefig(OUT / f"fig_{tag}_{group_col}_heatmap.png", dpi=300); plt.close(fig)


def main():
    atlas = sc.read_h5ad(C.WORK / "atlas_GSE184362.h5ad")
    thy = sc.read_h5ad(C.WORK / "thy_GSE184362.h5ad")

    # A. 全图谱
    if not (OUT / "atlas_group_results.csv").exists():
        a = sc_io.stratified_subsample(atlas, "cell_type", C.MAX_CELLS_ATLAS, C.SEED)
        grp, _ = run_scdrs(a, all_traits(), ["cell_type", "tissue"], "atlas")
        heatmap(grp, "cell_type", "atlas")
    # B. 甲状腺滤泡细胞
    if not (OUT / "thy_group_results.csv").exists():
        t = sc_io.stratified_subsample(thy, "thy_state", C.MAX_CELLS_THY, C.SEED)
        t.obs[["thy_state", "tissue", "sample", "patient"]].to_csv(OUT / "thy_scored_cells_meta.csv.gz")
        grp, scores = run_scdrs(t, C.MAIN_TRAITS + C.SENS_TRAITS, ["thy_state", "tissue"], "thy",
                                gene_traits=(C.MAL, C.NOD))
        heatmap(grp, "thy_state", "thy")
        heatmap(grp, "tissue", "thy")
        for tr in (C.MAL, C.NOD):
            if tr in scores:
                t.obs[f"scDRS_{tr}"] = scores[tr]["norm_score"].reindex(t.obs_names).values
        cols = [c for c in t.obs.columns if c.startswith("scDRS_")]
        if cols:
            sc.pl.umap(t, color=cols, cmap="RdBu_r", show=False, frameon=False)
            plt.savefig(OUT / "fig_thy_umap_scDRS.png", dpi=300, bbox_inches="tight"); plt.close()
    log.info("Step04 完成")


if __name__ == "__main__":
    main()
