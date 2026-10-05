# -*- coding: utf-8 -*-
"""
全局配置。一般只需要确认 BASE 路径；其余参数保持默认即可。
"""
from pathlib import Path

# ====== 1. 路径（与你电脑上的文件夹一致，注意文件夹名是 tyhroid） ======
BASE = Path(r"D:\tyhroid")
GWAS_DIR = BASE / "01_gwas"
MAGMA_DIR = BASE / "02_magma"
SC_DIR = BASE / "03_scrna"
TCGA_DIR = BASE / "04_tcga"
CLIN_DIR = BASE / "05_Molecular classification and progression-free survival"

WORK = BASE / "work"        # 中间文件（解压、缓存），可随时删除重算
RES = BASE / "results"      # 全部分析结果
LOG = RES / "logs"

# ====== 2. GWAS ======
GWAS_BUILD = "hg38"          # 该 GWAS 论文统一使用 hg38
MIN_MAF = 0.005              # 有频率列时过滤极低频变异
# 本研究使用的性状（键名由 step00 自动识别生成）
MAIN_TRAITS = ["ThCvsBNG_EUR", "BNG_EUR", "ThC_EUR"]     # 主分析
SENS_TRAITS = ["ThCvsBNG_MIX", "BNG_MIX", "ThC_MIX"]     # 敏感性分析（混合人群）
CTRL_TRAITS = ["GD_EUR", "LT_EUR"]                       # 内置对照（自身免疫）
MAL = "ThCvsBNG_EUR"     # 恶性特异风险
NOD = "BNG_EUR"          # 长结节风险

# ====== 3. MAGMA ======
MAGMA_WINDOW = "35,10"       # 基因上游35kb、下游10kb

# ====== 4. scDRS ======
SCDRS_TOP_GENES = 1000
SCDRS_N_CTRL = 1000
MAX_CELLS_ATLAS = 60000      # 全图谱打分时的最大细胞数（分层抽样，节省内存）
MAX_CELLS_THY = 60000        # 甲状腺滤泡细胞打分的最大细胞数
SEED = 0

# ====== 5. 单细胞质控 ======
MIN_GENES = 200
MAX_GENES = 7000
MAX_PCT_MT = 25
LEIDEN_RES_ALL = 0.8
LEIDEN_RES_THY = 0.6

# ====== 6. 基因程序 ======
PROGRAM_SIZE = 100           # 遗传驱动恶性程序（GIMP）基因数
MIN_EXPR_FRAC = 0.10         # 基因至少在10%的滤泡细胞中表达
N_RANDOM_SETS = 1000         # TCGA 随机基因集零分布次数

# ====== 7. TCGA 临床列名（自动识别失败时，可在此手动填写列名） ======
TCGA_COL_OVERRIDE = {
    # "N_STAGE": "PATH_N_STAGE",
    # "ETE": "EXTRATHYROIDAL_EXTENSION",
}

THREADS = 4
