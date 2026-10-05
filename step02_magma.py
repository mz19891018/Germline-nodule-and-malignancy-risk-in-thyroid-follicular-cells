# -*- coding: utf-8 -*-
"""
Step 02  MAGMA 基因水平分析 + 生成 scDRS 基因集
输出：
  results/02_magma/<trait>_genes.csv        每个性状的基因水平 Z 值
  results/02_magma/scdrs_genesets.gs        scDRS 输入（每性状 top1000 基因及权重）
  results/02_magma/gene_overlap_summary.csv 恶性特异 vs 长结节 风险基因的重叠
  results/02_magma/fig_Z_scatter.png
"""
import subprocess
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats
import config as C
from utils import ensure, get_logger, load_json, all_traits, bh_fdr

log = get_logger("step02_magma")
OUT = ensure(C.RES / "02_magma")
MOUT = ensure(C.WORK / "magma_out")
REF = C.WORK / "magma_ref"


def run(cmd):
    log.info("运行: " + " ".join(cmd))
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        log.error(r.stdout[-3000:] + "\n" + r.stderr[-3000:])
        raise SystemExit("MAGMA 运行失败，请把 results/logs/step02_magma.log 发给我")
    return r


def annotate(mode, paths, info):
    prefix = REF / f"annot_{mode}"
    if Path(str(prefix) + ".genes.annot").exists():
        return str(prefix) + ".genes.annot"
    snploc = paths["g1000_prefix"] + ".bim" if mode == "rsid" else str(REF / "g1000_eur_hg38id.snploc")
    run([paths["magma_exe"], "--annotate", f"window={C.MAGMA_WINDOW}", "--snp-loc", snploc,
         "--gene-loc", paths["gene_loc"], "--out", str(prefix)])
    return str(prefix) + ".genes.annot"


def gene_analysis(key, mode, paths, info, annot):
    prefix = MOUT / key
    if Path(str(prefix) + ".genes.out").exists():
        return str(prefix) + ".genes.out"
    bfile = paths["g1000_prefix"] if mode == "rsid" else info["_ref_hg38id"]
    run([paths["magma_exe"], "--bfile", bfile, "--pval", str(C.WORK / "magma_input" / f"{key}.txt"),
         "use=SNP,P", "ncol=N", "--gene-annot", annot, "--out", str(prefix)])
    return str(prefix) + ".genes.out"


def read_genes(path, gl):
    df = pd.read_csv(path, sep=r"\s+")
    df = df.merge(gl, left_on="GENE", right_on="ENTREZ", how="left")
    df = df.dropna(subset=["SYMBOL"]).drop_duplicates("SYMBOL")
    df["FDR"] = bh_fdr(df["P"].values)
    df["BONF_SIG"] = df["P"] < 0.05 / len(df)
    return df[["SYMBOL", "GENE", "CHR", "START", "STOP", "NSNPS", "N", "ZSTAT", "P", "FDR", "BONF_SIG"]]


def main():
    paths = load_json(C.RES / "00_check" / "magma_paths.json")
    info = load_json(C.RES / "00_check" / "gwas_prep_info.json")
    gl = pd.read_csv(paths["gene_loc"], sep=r"\s+", header=None,
                     names=["ENTREZ", "CHR_", "S", "E", "STRAND", "SYMBOL"])[["ENTREZ", "SYMBOL"]]
    tables = {}
    for key in all_traits():
        if key not in info:
            continue
        mode = info[key]["mode"]
        annot = annotate(mode, paths, info)
        g = read_genes(gene_analysis(key, mode, paths, info, annot), gl)
        g.to_csv(OUT / f"{key}_genes.csv", index=False)
        tables[key] = g
        log.info(f"{key}: {len(g)} 个基因；Bonferroni显著 {g.BONF_SIG.sum()} 个")

    # scDRS 基因集：每个性状取 Z 最大的 top N 基因，权重为 Z（截断在 0.01~10）
    with open(OUT / "scdrs_genesets.gs", "w") as f:
        f.write("TRAIT\tGENESET\n")
        for key, g in tables.items():
            top = g.sort_values("ZSTAT", ascending=False).head(C.SCDRS_TOP_GENES)
            w = top["ZSTAT"].clip(0.01, 10)
            f.write(key + "\t" + ",".join(f"{s}:{x:.4f}" for s, x in zip(top.SYMBOL, w)) + "\n")
    log.info("已生成 scDRS 基因集 scdrs_genesets.gs")

    # 恶性特异 vs 长结节 风险基因比较
    rows = []
    pairs = [(C.MAL, C.NOD), (C.MAL, "ThC_EUR"), ("ThC_EUR", C.NOD), (C.MAL, "GD_EUR"), (C.MAL, "LT_EUR")]
    for a, b in pairs:
        if a not in tables or b not in tables:
            continue
        m = tables[a].merge(tables[b], on="SYMBOL", suffixes=("_a", "_b"))
        rho, p = stats.spearmanr(m.ZSTAT_a, m.ZSTAT_b)
        ta = set(tables[a].nlargest(C.SCDRS_TOP_GENES, "ZSTAT").SYMBOL)
        tb = set(tables[b].nlargest(C.SCDRS_TOP_GENES, "ZSTAT").SYMBOL)
        k, M = len(ta & tb), len(m)
        p_hyper = stats.hypergeom.sf(k - 1, M, len(ta), len(tb))
        rows.append(dict(trait_a=a, trait_b=b, n_genes=M, spearman_rho=rho, spearman_p=p,
                         top_overlap=k, jaccard=k / len(ta | tb), hypergeom_p=p_hyper))
    pd.DataFrame(rows).to_csv(OUT / "gene_overlap_summary.csv", index=False)

    if C.MAL in tables and C.NOD in tables:
        m = tables[C.MAL].merge(tables[C.NOD], on="SYMBOL", suffixes=("_MAL", "_NOD"))
        spec = m[(m.FDR_MAL < 0.05) & (m.P_NOD > 0.05)].sort_values("P_MAL")
        spec.to_csv(OUT / "malignancy_specific_genes.csv", index=False)
        log.info(f"恶性特异基因（MAL FDR<0.05 且 NOD P>0.05）: {len(spec)} 个")
        import matplotlib; matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(5, 5))
        ax.scatter(m.ZSTAT_NOD, m.ZSTAT_MAL, s=3, alpha=.3, c="grey")
        ax.scatter(spec.ZSTAT_NOD, spec.ZSTAT_MAL, s=6, c="crimson", label="malignancy-specific")
        for r in spec.head(15).itertuples():
            ax.text(r.ZSTAT_NOD, r.ZSTAT_MAL, r.SYMBOL, fontsize=6)
        ax.set_xlabel("MAGMA Z: benign nodular goiter"); ax.set_ylabel("MAGMA Z: thyroid cancer vs BNG")
        ax.legend(fontsize=7); fig.tight_layout(); fig.savefig(OUT / "fig_Z_scatter.png", dpi=300)
    log.info("Step02 完成")


if __name__ == "__main__":
    main()
