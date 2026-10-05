# -*- coding: utf-8 -*-
import logging, sys, re, json, gzip, tarfile, zipfile, shutil
from pathlib import Path
import numpy as np
import pandas as pd
import config as C

# 兼容：scDRS 1.0.x 使用了 NumPy 2.0 已移除的 np.float_
if not hasattr(np, "float_"):
    np.float_ = np.float64


def ensure(p):
    p = Path(p)
    p.mkdir(parents=True, exist_ok=True)
    return p


def get_logger(name):
    ensure(C.LOG)
    logger = logging.getLogger(name)
    if logger.handlers:
        return logger
    logger.setLevel(logging.INFO)
    fmt = logging.Formatter("%(asctime)s | %(levelname)s | %(message)s", "%H:%M:%S")
    fh = logging.FileHandler(C.LOG / f"{name}.log", encoding="utf-8")
    fh.setFormatter(fmt)
    sh = logging.StreamHandler(sys.stdout)
    sh.setFormatter(fmt)
    logger.addHandler(fh)
    logger.addHandler(sh)
    return logger


def find_col(columns, candidates, contains=True):
    """按候选名（不区分大小写）在列中查找；先精确，后包含。"""
    low = {c.lower(): c for c in columns}
    for cand in candidates:
        if cand.lower() in low:
            return low[cand.lower()]
    if contains:
        for cand in candidates:
            for c in columns:
                if cand.lower() in c.lower():
                    return c
    return None


def open_text(path):
    path = str(path)
    if path.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, "r", encoding="utf-8", errors="replace")


def head_lines(path, n=5):
    out = []
    with open_text(path) as f:
        for i, line in enumerate(f):
            if i >= n:
                break
            out.append(line.rstrip("\n"))
    return out


def extract_archive(src, dst, logger=None):
    """解压 zip / tar / tar.gz；已解压则跳过。"""
    src, dst = Path(src), Path(dst)
    flag = dst / ".extracted"
    if flag.exists():
        return dst
    ensure(dst)
    if logger:
        logger.info(f"解压 {src.name} -> {dst}")
    name = src.name.lower()
    if name.endswith(".zip"):
        with zipfile.ZipFile(src) as z:
            z.extractall(dst)
    elif name.endswith((".tar", ".tar.gz", ".tgz")):
        with tarfile.open(src) as t:
            t.extractall(dst)
    else:
        raise ValueError(f"无法识别的压缩格式: {src}")
    flag.write_text("ok")
    return dst


def bh_fdr(p):
    p = np.asarray(p, dtype=float)
    out = np.full_like(p, np.nan)
    ok = ~np.isnan(p)
    if ok.sum() == 0:
        return out
    from statsmodels.stats.multitest import multipletests
    out[ok] = multipletests(p[ok], method="fdr_bh")[1]
    return out


def save_json(obj, path):
    ensure(Path(path).parent)
    with open(path, "w", encoding="utf-8") as f:
        json.dump(obj, f, ensure_ascii=False, indent=2, default=str)


def load_json(path):
    with open(path, "r", encoding="utf-8") as f:
        return json.load(f)


def trait_file_map():
    """读取 step00 生成的性状-文件对应表。"""
    manual = C.RES / "00_check" / "gwas_map_manual.csv"
    auto = C.RES / "00_check" / "gwas_inventory.csv"
    path = manual if manual.exists() else auto
    df = pd.read_csv(path)
    df = df[df["key"].notna() & (df["key"].astype(str) != "")]
    # 同一 key 多个文件时取样本量最大的
    df["N"] = pd.to_numeric(df["N"], errors="coerce")
    df = df.sort_values("N", ascending=False).drop_duplicates("key")
    return {r.key: r for r in df.itertuples()}


def all_traits():
    return C.MAIN_TRAITS + C.SENS_TRAITS + C.CTRL_TRAITS
