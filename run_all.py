# -*- coding: utf-8 -*-
"""
一键运行全部分析：  python run_all.py
从某一步开始：      python run_all.py 4      （从 Step04 开始）
每一步都会自动跳过已完成的部分，中断后重新运行即可续跑。
"""
import sys, time, traceback, importlib
import config as C
from utils import ensure

STEPS = [
    ("step00_check", "数据清点、解压、识别GWAS"),
    ("step01_gwas_prep", "GWAS 整理为 MAGMA 输入"),
    ("step02_magma", "MAGMA 基因水平分析"),
    ("step03_scrna_build", "构建单细胞图谱 GSE184362"),
    ("step04_scdrs", "scDRS 单细胞遗传风险打分"),
    ("step05_contrast_program", "恶性特异状态 + 基因程序 GIMP"),
    ("step06_replication", "GSE193581 独立验证"),
    ("step07_tcga", "TCGA 临床验证"),
    ("step08_report", "汇总报告"),
]


def main():
    start = int(sys.argv[1]) if len(sys.argv) > 1 else 0
    ensure(C.RES); ensure(C.WORK)
    for i, (mod, desc) in enumerate(STEPS):
        if i < start:
            continue
        print(f"\n{'=' * 70}\n▶ Step{i:02d}  {desc}\n{'=' * 70}")
        t0 = time.time()
        try:
            importlib.import_module(mod).main()
        except SystemExit as e:
            print(f"\n⏸ Step{i:02d} 暂停/中止：{e}")
            print(f"处理后运行： python run_all.py {i}")
            sys.exit(1)
        except Exception:
            traceback.print_exc()
            print(f"\n✖ Step{i:02d} 出错。请把 results/logs/{mod}.log 的最后50行发给我。")
            print(f"修复后从这一步续跑： python run_all.py {i}")
            sys.exit(1)
        print(f"✔ Step{i:02d} 完成，用时 {(time.time() - t0) / 60:.1f} 分钟")
    print("\n全部完成！请把 results/08_report/summary_report.md 发给我。")


if __name__ == "__main__":
    main()
