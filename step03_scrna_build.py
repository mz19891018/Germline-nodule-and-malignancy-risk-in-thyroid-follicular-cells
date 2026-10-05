# -*- coding: utf-8 -*-
"""
Step 03  构建 GSE184362 单细胞图谱
- 读取全部样本 -> 质控 -> 去双细胞 -> 标准化 -> Harmony整合 -> 聚类 -> 细胞类型注释
- 提取甲状腺滤泡细胞（Thyrocyte），再细分为细胞状态（Thy_0, Thy_1 ...）
首次运行会生成 results/03_scrna/sample_metadata.csv 并暂停，
请核对每个样本的组织来源（Paratumor/Tumor/LN/DistantMet）和患者编号后，再重新运行。
"""
import shutil
import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import config as C
from utils import ensure, get_logger
import sc_io

log = get_logger("step03_scrna_build")
OUT = ensure(C.RES / "03_scrna")
RAW = C.WORK / "GSE184362_RAW"
ATLAS = C.WORK / "atlas_GSE184362.h5ad"
THY = C.WORK / "thy_GSE184362.h5ad"
sc.settings.figdir = str(OUT)


def sample_table(groups):
    meta_path = OUT / "sample_metadata.csv"
    if meta_path.exists():
        try:
            m = pd.read_csv(meta_path, dtype=str, encoding="utf-8-sig")
        except UnicodeDecodeError:          # Excel 另存后可能是 GBK 编码
            m = pd.read_csv(meta_path, dtype=str, encoding="gbk")
        m.columns = m.columns.str.strip()
        if (m.get("confirmed", pd.Series(["no"])).str.lower() == "yes").all():
            return m
        raise SystemExit("请打开 results/03_scrna/sample_metadata.csv，核对 tissue 和 patient 两列，"
                         "确认无误后把 confirmed 列全部改为 yes，保存后重新运行。")
    rows = []
    for s, fs in groups.items():
        name = " ".join(sorted(f.name for f in fs))
        rows.append(dict(sample=s, files=name[:300], tissue=sc_io.guess_tissue(name),
                         patient=sc_io.guess_patient(name.replace(s, "")), confirmed="no"))
    m = pd.DataFrame(rows)
    m.to_csv(meta_path, index=False, encoding="utf-8-sig")
    log.info(m.to_string())
    raise SystemExit("已自动生成 results/03_scrna/sample_metadata.csv（猜测的组织来源和患者编号）。"
                     "请对照 GEO 页面 GSE184362 的样本说明核对 tissue/patient 两列，"
                     "把 confirmed 改为 yes 后重新运行。tissue 只能填：Paratumor / Tumor / LN / DistantMet")


def load_all():
    groups = sc_io.discover_samples(RAW)
    log.info(f"发现 {len(groups)} 个样本组：{list(groups)[:20]}")
    meta = sample_table(groups).set_index("sample")
    adatas = []
    for s, fs in groups.items():
        a = sc_io.read_sample(fs)
        if a is None:
            log.warning(f"{s}: 未识别出表达矩阵，跳过")
            continue
        a.obs_names = [f"{s}_{x}" for x in a.obs_names]
        a.obs["sample"] = s
        a.obs["tissue"] = meta.loc[s, "tissue"]
        a.obs["patient"] = meta.loc[s, "patient"]
        log.info(f"{s}: {a.n_obs} 细胞 x {a.n_vars} 基因 | {meta.loc[s, 'tissue']}")
        adatas.append(a)
    a = ad.concat(adatas, join="outer", fill_value=0)
    a.obs_names_make_unique()
    return a


def thyrocyte_states(a):
    t = a[a.obs.cell_type == "Thyrocyte"].copy()
    log.info(f"甲状腺滤泡细胞：{t.n_obs} 个")
    t.X = t.layers["counts"].copy()
    sc.pp.normalize_total(t, target_sum=1e4); sc.pp.log1p(t)
    t = sc_io.embed_cluster(t, C.LEIDEN_RES_THY, key="thy_cluster")
    t.obs["thy_state"] = ("Thy_" + t.obs["thy_cluster"].astype(str)).astype("category")
    g = [x for x in sc_io.TDS_GENES if x in t.var_names]
    sc.tl.score_genes(t, g, score_name="TDS", random_state=C.SEED)
    sc.tl.rank_genes_groups(t, "thy_state", method="wilcoxon")
    res = sc.get.rank_genes_groups_df(t, None)
    res.groupby("group").head(50).to_csv(OUT / "thy_state_markers_top50.csv", index=False)
    comp = pd.crosstab(t.obs.thy_state, t.obs.tissue)
    comp.to_csv(OUT / "thy_state_by_tissue_counts.csv")
    (comp.div(comp.sum(1), axis=0)).to_csv(OUT / "thy_state_by_tissue_fraction.csv")
    t.obs.groupby("thy_state", observed=True)["TDS"].describe().to_csv(OUT / "thy_state_TDS.csv")
    return t


def main():
    if ATLAS.exists() and THY.exists():
        log.info("图谱已存在，跳过（如需重建请删除 work/atlas_GSE184362.h5ad）")
        return
    a = load_all()
    a = sc_io.qc_filter(a, log)
    a = sc_io.remove_doublets(a, log)
    a = sc_io.normalize(a)
    a = sc_io.embed_cluster(a, C.LEIDEN_RES_ALL)
    z = sc_io.annotate_by_markers(a)
    z.to_csv(OUT / "cluster_annotation_marker_z.csv")
    pd.crosstab(a.obs.cell_type, a.obs.tissue).to_csv(OUT / "celltype_by_tissue.csv")
    for col in ["cell_type", "tissue", "sample"]:
        sc.pl.umap(a, color=col, show=False, frameon=False, save=f"_atlas_{col}.png")
    sc.pl.dotplot(a, {k: [g for g in v if g in a.var_names] for k, v in sc_io.MARKERS.items()},
                  groupby="cell_type", show=False, save="_atlas_markers.png")
    t = thyrocyte_states(a)
    for col in ["thy_state", "tissue", "TDS"]:
        sc.pl.umap(t, color=col, show=False, frameon=False, save=f"_thy_{col}.png")
    a.obs[["sample", "tissue", "patient", "cell_type"]].join(
        t.obs[["thy_state", "TDS"]]).to_csv(OUT / "cell_metadata.csv.gz")
    a.write_h5ad(ATLAS, compression="gzip")
    t.write_h5ad(THY, compression="gzip")
    pd.Series({"n_cells_atlas": a.n_obs, "n_genes": a.n_vars, "n_thyrocytes": t.n_obs,
               "n_samples": a.obs["sample"].nunique(), "n_patients": a.obs["patient"].nunique(),
               "n_thy_states": t.obs.thy_state.nunique()}).to_csv(OUT / "atlas_summary.csv")
    log.info("Step03 完成")


if __name__ == "__main__":
    main()
