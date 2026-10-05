# -*- coding: utf-8 -*-
"""
Step 08  汇总全部关键结果 -> results/08_report/summary_report.md
写文章时只需把这个文件（以及 results 文件夹打包）发给我。
"""
import json, platform, importlib
from pathlib import Path
import pandas as pd
import config as C
from utils import ensure, get_logger, save_json

log = get_logger("step08_report")
OUT = ensure(C.RES / "08_report")
R = C.RES


def tbl(path, n=30, cols=None):
    p = Path(path)
    if not p.exists():
        return f"_（缺失：{p.relative_to(C.BASE)}）_\n"
    try:
        df = pd.read_csv(p)
    except pd.errors.EmptyDataError:
        return "_（无结果）_\n"
    if cols:
        df = df[[c for c in cols if c in df.columns]]
    return df.head(n).to_markdown(index=False, floatfmt=".4g") + "\n"


def js(path):
    p = Path(path)
    return "```\n" + p.read_text(encoding="utf-8") + "\n```\n" if p.exists() else f"_（缺失：{p.name}）_\n"


def main():
    md = ["# 结果汇总报告（自动生成）\n"]
    md += ["## 0. GWAS 数据\n", tbl(R / "00_check/gwas_inventory.csv", 50,
                                    ["accession", "key", "N", "trait_catalog"]),
           js(R / "00_check/gwas_prep_info.json")]
    md += ["## 1. MAGMA 基因水平结果\n", tbl(R / "02_magma/gene_overlap_summary.csv"),
           "### 恶性特异风险基因（前30）\n",
           tbl(R / "02_magma/malignancy_specific_genes.csv", 30, ["SYMBOL", "ZSTAT_MAL", "P_MAL", "FDR_MAL", "ZSTAT_NOD", "P_NOD"])]
    md += ["## 2. 单细胞图谱\n", tbl(R / "03_scrna/atlas_summary.csv"), tbl(R / "03_scrna/celltype_by_tissue.csv"),
           tbl(R / "03_scrna/thy_state_by_tissue_fraction.csv"), tbl(R / "03_scrna/thy_state_TDS.csv")]
    md += ["## 3. scDRS：全图谱细胞类型\n",
           tbl(R / "04_scdrs/atlas_group_results.csv", 200, ["trait", "group_col", "group", "n_cell", "assoc_mcp", "assoc_fdr", "hetero_mcp"]),
           "## 4. scDRS：滤泡细胞状态\n",
           tbl(R / "04_scdrs/thy_group_results.csv", 200, ["trait", "group_col", "group", "n_cell", "assoc_mcp", "assoc_fdr", "hetero_mcp"])]
    md += ["## 5. 恶性特异状态与组织来源（核心）\n", js(R / "05_program/residual_model.json"),
           "### 组织来源（患者内配对 + 混合模型）\n", tbl(R / "05_program/tissue_tests.csv"),
           "### 分化程度（样本固定效应）\n", tbl(R / "05_program/differentiation_tests.csv"),
           "### 细胞状态（描述性）\n", tbl(R / "05_program/malignancy_specific_states.csv", 50),
           "### GIMP/BNP 组织来源与分化\n", tbl(R / "05_program/program_tissue_tests.csv"),
           tbl(R / "05_program/program_differentiation_tests.csv"),
           "### GIMP 基因（前40）\n", tbl(R / "05_program/GIMP_genes.csv", 40),
           "### BNP 基因（前40）\n", tbl(R / "05_program/BNP_genes.csv", 40)]
    for nm in ["GIMP", "BNP"]:
        md += [f"### {nm} 富集（前15）\n", tbl(R / f"05_program/enrichment_{nm}.csv", 15,
                                             ["Gene_set", "Term", "Overlap", "Adjusted P-value", "Genes"])]
    md += ["## 6. GSE193581 独立验证\n", js(R / "06_replication/replication_summary.json"),
           tbl(R / "06_replication/replication_tissue_tests.csv"), tbl(R / "06_replication/sample_metadata_auto.csv", 30),
           tbl(R / "06_replication/celltype_by_tissue.csv")]
    md += ["## 7. TCGA 临床验证\n", tbl(R / "07_tcga/clinical_column_mapping.csv"),
           tbl(R / "07_tcga/tcga_associations.csv", 50), js(R / "07_tcga/random_geneset_null_N1.json"),
           js(R / "07_tcga/tumor_vs_normal.json")]
    (OUT / "summary_report.md").write_text("\n".join(md), encoding="utf-8")

    vers = {"python": platform.python_version()}
    for p in ["scanpy", "anndata", "scdrs", "harmonypy", "statsmodels", "lifelines", "gseapy", "numpy", "pandas"]:
        try:
            vers[p] = importlib.import_module(p).__version__
        except Exception:
            pass
    params = {k: getattr(C, k) for k in dir(C) if k.isupper() and not isinstance(getattr(C, k), Path)}
    save_json({"versions": vers, "parameters": params}, OUT / "methods_versions_parameters.json")
    log.info("报告已生成：results/08_report/summary_report.md")


if __name__ == "__main__":
    main()
