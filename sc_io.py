# -*- coding: utf-8 -*-
"""单细胞数据通用读取与处理函数（自动识别 10x mtx / h5 / 文本矩阵）"""
import re, tarfile, gzip
from pathlib import Path
import numpy as np
import pandas as pd
import scipy.io, scipy.sparse as sp
import anndata as ad
import scanpy as sc
import config as C

MARKERS = {
    "Thyrocyte": ["TG", "TPO", "TSHR", "PAX8", "NKX2-1", "FOXE1", "KRT18", "EPCAM", "KRT19"],
    "T_NK": ["CD3D", "CD3E", "CD2", "NKG7", "GNLY", "IL7R"],
    "B": ["MS4A1", "CD79A", "CD19", "BANK1"],
    "Plasma": ["JCHAIN", "MZB1", "IGHG1", "IGKC"],
    "Myeloid": ["CD68", "LYZ", "CD14", "C1QA", "AIF1", "FCGR3A"],
    "Mast": ["TPSAB1", "TPSB2", "CPA3", "MS4A2"],   # 不用KIT：正常甲状腺滤泡细胞也高表达KIT
    "Fibroblast": ["DCN", "COL1A1", "COL1A2", "LUM", "PDGFRA"],
    "Endothelial": ["PECAM1", "VWF", "PLVAP", "CDH5"],
    "SMC_Pericyte": ["ACTA2", "RGS5", "MYH11", "PDGFRB"],
}
# 甲状腺分化评分（TDS，TCGA 2014）基因
TDS_GENES = ["DIO1", "DIO2", "DUOX1", "DUOX2", "FOXE1", "GLIS3", "NKX2-1", "PAX8", "SLC26A4",
             "SLC26A7", "SLC5A5", "SLC5A8", "TG", "THRA", "THRB", "TPO", "TSHR"]


def _expand_nested(raw_dir):
    for p in list(raw_dir.rglob("*")):
        if p.is_file() and p.name.endswith((".tar", ".tar.gz", ".tgz")):
            dst = p.parent / (p.name.split(".tar")[0] + "_x")
            if not dst.exists():
                with tarfile.open(p) as t:
                    t.extractall(dst)


def discover_samples(raw_dir):
    raw_dir = Path(raw_dir)
    _expand_nested(raw_dir)
    files = [p for p in raw_dir.rglob("*") if p.is_file() and not p.name.startswith(".")
             and not p.name.endswith((".tar", ".tar.gz", ".tgz", ".extracted"))]
    groups = {}
    for p in files:
        m = re.match(r"(GSM\d+)", p.name)
        key = m.group(1) if m else p.parent.name
        groups.setdefault(key, []).append(p)
    return groups


def _first(files, pats):
    for f in files:
        n = f.name.lower()
        if any(x in n for x in pats):
            return f
    return None


def _read_lines(p):
    op = gzip.open if str(p).endswith(".gz") else open
    with op(p, "rt") as f:
        return [l.rstrip("\n") for l in f]


def read_sample(files):
    """返回 AnnData（细胞 x 基因，原始计数）"""
    mtx = _first(files, ["matrix.mtx", ".mtx"])
    h5 = _first(files, [".h5"])
    if mtx is not None:
        bc = _first(files, ["barcode"])
        ft = _first(files, ["features", "genes"])
        X = scipy.io.mmread(str(mtx)).T.tocsr()
        barcodes = [l.split("\t")[0] for l in _read_lines(bc)]
        feats = pd.read_csv(ft, sep="\t", header=None)
        genes = feats[1] if feats.shape[1] > 1 else feats[0]
        if X.shape[0] != len(barcodes):          # 极少数矩阵已是 细胞x基因
            X = X.T.tocsr()
        a = ad.AnnData(X=X.astype(np.float32), obs=pd.DataFrame(index=barcodes),
                       var=pd.DataFrame(index=genes.astype(str).values))
    elif h5 is not None:
        a = sc.read_10x_h5(str(h5))
    else:
        cands = [f for f in files if re.search(r"\.(txt|csv|tsv)(\.gz)?$", f.name.lower())]
        if not cands:
            return None
        f = max(cands, key=lambda x: x.stat().st_size)
        sep = "," if ".csv" in f.name.lower() else "\t"
        df = pd.read_csv(f, sep=sep, index_col=0)
        probe = {"GAPDH", "ACTB", "MALAT1", "B2M", "TG"}
        genes_in_rows = len(probe & set(map(str, df.index))) >= len(probe & set(map(str, df.columns)))
        if genes_in_rows:
            df = df.T
        df = df.apply(pd.to_numeric, errors="coerce").fillna(0)
        a = ad.AnnData(X=sp.csr_matrix(df.values.astype(np.float32)),
                       obs=pd.DataFrame(index=df.index.astype(str)),
                       var=pd.DataFrame(index=df.columns.astype(str)))
    a.var_names_make_unique()
    return a


def guess_tissue(name):
    n = name.lower()
    tok = r"(^|[_\-\.\s])"
    end = r"\d*($|[_\-\.\s])"
    if re.search(tok + r"(ln|lymph|lnm|node|leftln|rightln)" + end, n) or "lymph" in n:
        return "LN"
    if re.search(tok + r"(m|met|meta|metastasis|distant|rair|lung|bone)" + end, n):
        return "DistantMet"
    if re.search(tok + r"(p|para|paratumor|n|nt|normal|norm|adj|adjacent)" + end, n) or "normal" in n \
            or "para" in n or re.search(r"norm\d", n):
        return "Paratumor"
    if "atc" in n or "anaplastic" in n:
        return "ATC"
    return "Tumor"


def guess_patient(name):
    m = re.search(r"((?:ptc|patient|pt|p|case)[_\-]?\d+)", name.lower())
    return m.group(1).upper() if m else name.split("_")[0]


def qc_filter(a, logger=None):
    a.var["mt"] = a.var_names.str.upper().str.startswith("MT-")
    sc.pp.calculate_qc_metrics(a, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)
    n0 = a.n_obs
    a = a[(a.obs.n_genes_by_counts >= C.MIN_GENES) & (a.obs.n_genes_by_counts <= C.MAX_GENES)
          & (a.obs.pct_counts_mt <= C.MAX_PCT_MT)].copy()
    sc.pp.filter_genes(a, min_cells=10)
    if logger:
        logger.info(f"质控：{n0} -> {a.n_obs} 个细胞")
    return a


def remove_doublets(a, logger=None):
    try:
        sc.pp.scrublet(a, batch_key="sample")
        n0 = a.n_obs
        a = a[~a.obs["predicted_doublet"].astype(bool)].copy()
        if logger:
            logger.info(f"Scrublet 去除双细胞：{n0} -> {a.n_obs}")
    except Exception as e:
        if logger:
            logger.warning(f"Scrublet 未运行（{e}），跳过去双细胞")
    return a


def normalize(a):
    a.layers["counts"] = a.X.copy()
    sc.pp.normalize_total(a, target_sum=1e4)
    sc.pp.log1p(a)
    return a


def embed_cluster(a, resolution, batch="sample", n_hvg=2000, key="leiden"):
    sc.pp.highly_variable_genes(a, n_top_genes=n_hvg, batch_key=batch, flavor="seurat")
    h = a[:, a.var.highly_variable].copy()
    sc.pp.scale(h, max_value=10)
    sc.tl.pca(h, n_comps=50, random_state=C.SEED)
    a.obsm["X_pca"] = h.obsm["X_pca"]
    rep = "X_pca"
    if a.obs[batch].nunique() > 1:
        try:
            sc.external.pp.harmony_integrate(a, batch, basis="X_pca", adjusted_basis="X_pca_harmony")
            rep = "X_pca_harmony"
        except Exception:
            pass
    sc.pp.neighbors(a, use_rep=rep, n_neighbors=20, random_state=C.SEED)
    sc.tl.umap(a, random_state=C.SEED)
    try:
        sc.tl.leiden(a, resolution=resolution, key_added=key, flavor="igraph", n_iterations=2,
                     directed=False, random_state=C.SEED)
    except TypeError:
        sc.tl.leiden(a, resolution=resolution, key_added=key, random_state=C.SEED)
    return a


def annotate_by_markers(a, cluster_key="leiden"):
    scores = {}
    for ct, genes in MARKERS.items():
        g = [x for x in genes if x in a.var_names]
        if len(g) >= 2:
            sc.tl.score_genes(a, g, score_name=f"mk_{ct}", random_state=C.SEED)
            scores[ct] = f"mk_{ct}"
    df = a.obs.groupby(cluster_key)[list(scores.values())].mean()
    z = (df - df.mean()) / (df.std() + 1e-9)
    lab = z.idxmax(axis=1).str.replace("mk_", "", regex=False)
    # 安全规则：甲状腺标志明显阳性（z>1）的簇一律判为滤泡细胞，防止被其他类型抢走
    if "mk_Thyrocyte" in z.columns:
        lab[(z["mk_Thyrocyte"] > 1) & (lab != "Thyrocyte")] = "Thyrocyte"
    a.obs["cell_type"] = a.obs[cluster_key].map(lab).astype("category")
    return z.assign(assigned=lab)


def stratified_subsample(a, key, max_cells, seed=0):
    if a.n_obs <= max_cells:
        return a
    rng = np.random.default_rng(seed)
    frac = max_cells / a.n_obs
    idx = []
    for g, sub in a.obs.groupby(key, observed=True).groups.items():
        sub = np.array(list(sub))
        k = max(min(len(sub), 50), int(round(len(sub) * frac)))
        idx += list(rng.choice(sub, size=min(k, len(sub)), replace=False))
    return a[idx].copy()
