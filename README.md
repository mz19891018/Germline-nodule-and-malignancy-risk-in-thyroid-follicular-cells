# Germline-nodule-and-malignancy-risk-in-thyroid-follicular-cells
Malignancy-specific germline liability is weaker and tracks the opposite, dedifferentiated BRAF-like axis. These results place the germline architecture of nodule formation and malignant transformation on a single follicular differentiation axis.
甲状腺结节恶性特异遗传风险・单细胞定位分析流程
一、放置代码
把整个 thyroid_pipeline 文件夹放到 D:\tyhroid\ 下，最终路径是 D:\tyhroid\thyroid_pipeline\。
数据文件夹保持你现在的样子，不要改名、不要手动解压（代码会自动解压到 D:\tyhroid\work\）。
二、安装运行环境（只需一次）
安装 Miniconda（https://docs.conda.io/en/latest/miniconda.html ，选 Windows 64-bit）。
打开开始菜单里的 Anaconda Prompt，依次输入：
conda create -n thyroid python=3.11 -y
conda activate thyroid
cd /d D:\tyhroid\thyroid_pipeline
pip install -r requirements.txt

三、一键运行
conda activate thyroid
cd /d D:\tyhroid\thyroid_pipeline
python run_all.py

全程需要联网（识别 GWAS 编号、下载坐标转换文件、富集分析）。建议电脑不休眠。
中途会有一次有意的暂停（Step03）：
代码会生成 D:\tyhroid\results\03_scrna\sample_metadata.csv，里面是代码猜测的每个单细胞样本的
组织来源（Paratumor 癌旁 / Tumor 原发灶 / LN 淋巴结 / DistantMet 远处转移）和患者编号。
请对照 GEO 页面 GSE184362 的样本列表核对 tissue 和 patient 两列，改正后把 confirmed 列全部填 yes，保存，然后运行：
python run_all.py 3

（数字表示从第几步续跑。任何一步出错，修好后用同样方式从该步续跑，已完成的部分会自动跳过。）
可能出现的另一次暂停（Step00）： 如果 40 个 GWAS 编号没能全部自动识别，按屏幕提示在
results\00_check\gwas_map_manual_TEMPLATE.csv 的 key 列手动填写后另存为 gwas_map_manual.csv。
key 的写法：ThCvsBNG_EUR、BNG_EUR、ThC_EUR、ThCvsBNG_MIX、BNG_MIX、ThC_MIX、GD_EUR、LT_EUR。
四、各步骤与预计用时（32GB 内存电脑）
步骤
内容
预计用时
Step00
清点数据、解压、联网识别 40 个 GWAS
5–10 分钟
Step01
8 套 GWAS 整理为 MAGMA 格式；参考面板坐标转换（仅一次）
1–2 小时
Step02
MAGMA 基因水平分析
2–4 小时
Step03
构建单细胞图谱、注释细胞类型、滤泡细胞分状态
0.5–1 小时
Step04
scDRS 单细胞遗传风险打分
1–3 小时
Step05
恶性特异状态检验 + 遗传驱动恶性基因程序（GIMP）
10–30 分钟
Step06
GSE193581 独立验证
0.5–1 小时
Step07
TCGA 临床验证（淋巴结转移、腺外侵犯、分子分型、PFI）
10–30 分钟
Step08
生成汇总报告
1 分钟

五、结果在哪里
全部在 D:\tyhroid\results\，按步骤分文件夹，图片为 300 dpi PNG。
运行完成后，请把以下内容发给我用于撰写文章：
results\08_report\summary_report.md（最重要，汇总了全部关键数字）
results\08_report\methods_versions_parameters.json
整个 results 文件夹打包（不含 work 文件夹，work 只是中间文件）
六、出错怎么办
屏幕会提示出错的步骤。把 results\logs\ 里对应步骤 .log 文件的最后 50 行发给我即可。
