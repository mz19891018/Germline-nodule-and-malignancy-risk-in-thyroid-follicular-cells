# -*- coding: utf-8 -*-
"""
Step 07  TCGA 临床验证
- cBioPortal (TCGA Cell 2014) 表达 + 临床 + 分子评分（BRAF-RAS、TDS、ERK）
- TCGA-CDR 无进展间隔（PFI）
- Xena STAR counts：肿瘤 vs 配对正常组织（可选）
检验 GIMP（及对照 BNP）与淋巴结转移、腺外侵犯、T分期、分子分型、PFI 的关系；
并以随机基因集零分布检验特异性。
"""
import re
import numpy as np
import pandas as pd
from scipy import stats
import statsmodels.formula.api as smf
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
import config as C
from utils import ensure, get_logger, find_col, save_json

log = get_logger("step07_tcga")
OUT = ensure(C.RES / "07_tcga")
CB = C.WORK / "cbioportal"

COLMAP = {
    "AGE": ["AGE", "AGE_AT_DIAGNOSIS"],
    "SEX": ["SEX", "GENDER"],
    "N_STAGE": ["PATH_N_STAGE", "AJCC_NODES_PATHOLOGIC_PN", "PATHOLOGIC_N", "N_STAGE"],
    "T_STAGE": ["PATH_T_STAGE", "AJCC_TUMOR_PATHOLOGIC_PT", "PATHOLOGIC_T", "T_STAGE"],
    "STAGE": ["DISEASE_STAGE", "AJCC_PATHOLOGIC_TUMOR_STAGE", "PATH_STAGE", "PATHOLOGIC_STAGE", "TUMOR_STAGE"],
    "M_STAGE": ["PATH_M_STAGE"],
    "TERT": ["TERT_PROMOTER_MUTATION"],
    "RISK": ["RISK_GROUP"],
    "MACIS": ["MACIS_SCORE"],
    "BRS_CLASS": ["BRAFV600E_RAS"],
    "ETE": ["EXTRATHYROIDAL_EXTENSION", "EXTRATHYROIDAL"],
    "FOCALITY": ["TUMOR_FOCALITY", "FOCALITY", "MULTIFOCALITY"],
    "BRS": ["BRAFV600E_RAF_SCORE", "BRAF_RAS_SCORE", "BRS"],   # 数值型 BRAF-RAS 评分
    "TDS": ["DIFFERENTIATION_SCORE", "THYROID_DIFFERENTIATION_SCORE", "TDS"],
    "ERK": ["ERK_SCORE", "ERK"],
    "HIST": ["HISTOLOGICAL_SUBTYPE", "HISTOLOGICAL_TYPE", "CANCER_TYPE_DETAILED"],
}


def read_cbio(name_pat):
    f = [p for p in CB.rglob("*") if p.is_file() and re.search(name_pat, p.name, re.I)]
    if not f:
        return None, None
    f = sorted(f, key=lambda p: len(p.name))[0]
    na = ["", "NA", "NaN", "nan", "[Not Available]", "[Unknown]", "[Not Applicable]",
          "[Not Evaluated]", "[Discrepancy]", "[Completed]"]
    # 注意："None" 在腺外侵犯中是有效取值（无侵犯），不能被当作缺失
    return pd.read_csv(f, sep="\t", comment="#", low_memory=False, keep_default_na=False,
                       na_values=na), f.name


def load_clinical():
    pat, _ = read_cbio(r"^data_clinical_patient\.txt$")
    smp, _ = read_cbio(r"^data_clinical_sample\.txt$")
    clin = smp.merge(pat, on="PATIENT_ID", how="left") if pat is not None else smp
    mapping = {}
    for k, cands in COLMAP.items():
        col = C.TCGA_COL_OVERRIDE.get(k) or find_col(clin.columns, cands, contains=False) \
            or find_col(clin.columns, cands, contains=True)
        mapping[k] = col
    pd.Series(mapping).to_csv(OUT / "clinical_column_mapping.csv", header=["column_used"])
    pd.Series(list(clin.columns)).to_csv(OUT / "clinical_all_columns.csv", index=False)
    log.info(f"临床列映射：{mapping}")
    d = pd.DataFrame({"SAMPLE_ID": clin.SAMPLE_ID, "PATIENT_ID": clin.PATIENT_ID})
    g = lambda k: clin[mapping[k]] if mapping[k] else pd.Series(np.nan, index=clin.index)
    d["age"] = pd.to_numeric(g("AGE"), errors="coerce")
    d["male"] = g("SEX").astype(str).str.lower().str.startswith("m").astype(float)
    n = g("N_STAGE").astype(str).str.upper()
    d["N1"] = np.where(n.str.contains("N1"), 1.0, np.where(n.str.contains("N0"), 0.0, np.nan))
    t = g("T_STAGE").astype(str).str.upper()
    d["T34"] = np.where(t.str.contains("T3|T4"), 1.0, np.where(t.str.contains("T1|T2"), 0.0, np.nan))
    s = g("STAGE").astype(str).str.upper().str.replace("STAGE", "").str.strip()
    d["stage34"] = np.where(s.str.match(r"^(III|IV)"), 1.0, np.where(s.str.match(r"^(I|II)[ABC]?$"), 0.0, np.nan))
    e = g("ETE").astype(str).str.lower()
    d["ETE"] = np.where(e.str.contains("none|^no"), 0.0,
                        np.where(e.str.contains("minimal|moderate|advanced|yes|t3|t4"), 1.0, np.nan))
    f = g("FOCALITY").astype(str).str.lower()
    d["multifocal"] = np.where(f.str.contains("multi"), 1.0, np.where(f.str.contains("uni"), 0.0, np.nan))
    for k in ["BRS", "TDS", "ERK", "MACIS"]:
        d[k] = pd.to_numeric(g(k), errors="coerce")
    bc = g("BRS_CLASS").astype(str).str.lower()
    d["braf_like"] = np.where(bc.str.contains("braf"), 1.0, np.where(bc.str.contains("ras"), 0.0, np.nan))
    m_ = g("M_STAGE").astype(str).str.upper()
    d["M1"] = np.where(m_.str.contains("M1"), 1.0, np.where(m_.str.contains("M0"), 0.0, np.nan))
    tt = g("TERT").astype(str).str.upper()
    mut = tt.str.contains("C228T|C250T|MUT|YES").values
    cl = [p for p in CB.rglob("cases_TERT_mutation.txt")]
    if cl:                      # 检测过TERT的样本：有突变=1，其余=0；未检测=缺失
        import re as _re
        assessed = set(_re.findall(r"TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}-\d{2}", cl[0].read_text()))
        d["TERT_mut"] = np.where(mut, 1.0, np.where(clin.SAMPLE_ID.isin(assessed), 0.0, np.nan))
        log.info(f"TERT：检测样本 {len(assessed)}，突变 {int(mut.sum())}")
    else:
        d["TERT_mut"] = np.where(mut, 1.0, np.where(tt.str.contains("WT|WILD|NONE|^NO"), 0.0, np.nan))
    rk = g("RISK").astype(str).str.lower()
    d["risk_high"] = np.where(rk.str.contains("high"), 1.0, np.where(rk.str.contains("low|interm"), 0.0, np.nan))
    for k in ["BRS_CLASS", "TERT", "RISK", "M_STAGE", "STAGE"]:
        log.info(f"{k} 取值：{g(k).value_counts(dropna=False).head(8).to_dict()}")
    d["hist"] = g("HIST").astype(str)
    return d


def load_expr():
    e, name = read_cbio(r"(rsem|RNA_Seq_v2_expression|mrna_seq).*\.txt$")
    if e is None or "zscore" in name.lower():
        cands = [p for p in CB.rglob("*.txt") if re.search("rsem|rna_seq", p.name, re.I) and "zscore" not in p.name.lower()]
        e = pd.read_csv(cands[0], sep="\t", comment="#", low_memory=False); name = cands[0].name
    log.info(f"表达文件：{name}")
    gcol = "Hugo_Symbol" if "Hugo_Symbol" in e.columns else e.columns[0]
    e = e.drop(columns=[c for c in e.columns if c.lower().startswith("entrez")], errors="ignore")
    e = e.dropna(subset=[gcol]).groupby(gcol).mean(numeric_only=True)
    e = np.log2(e.clip(lower=0) + 1)
    e = e[(e > 1).mean(1) >= 0.2]
    return e


def score(z, genes):
    g = [x for x in genes if x in z.index]
    return z.loc[g].mean(0), len(g)


def logit_or(df, y, x, cov=("age", "male")):
    d = df[[y, x] + list(cov)].dropna()
    if d[y].nunique() < 2 or len(d) < 30:
        return dict(n=len(d), OR=np.nan, lo=np.nan, hi=np.nan, p=np.nan)
    m = smf.logit(f"{y} ~ {x} + " + " + ".join(cov), d).fit(disp=0)
    ci = m.conf_int().loc[x]
    return dict(n=len(d), n_events=int(d[y].sum()), OR=np.exp(m.params[x]), lo=np.exp(ci[0]),
                hi=np.exp(ci[1]), p=m.pvalues[x])


def survival(df, x, outdir_tag):
    from lifelines import CoxPHFitter, KaplanMeierFitter
    from lifelines.statistics import logrank_test
    cov = ["age", "male"] + (["stage34"] if df["stage34"].notna().sum() > 100 else [])
    d = df[["PFI", "PFI.time", x] + cov].dropna()
    d = d[d["PFI.time"] > 0]
    cph = CoxPHFitter().fit(d, "PFI.time", "PFI")
    s = cph.summary.loc[x]
    cu = CoxPHFitter().fit(d[["PFI", "PFI.time", x]], "PFI.time", "PFI").summary.loc[x]
    hi = d[x] > d[x].median()
    lr = logrank_test(d["PFI.time"][hi], d["PFI.time"][~hi], d["PFI"][hi], d["PFI"][~hi])
    fig, ax = plt.subplots(figsize=(4.5, 4))
    for lab, msk in [("High", hi), ("Low", ~hi)]:
        KaplanMeierFitter().fit(d["PFI.time"][msk] / 365.25, d["PFI"][msk], label=f"{x} {lab}").plot_survival_function(ax=ax, ci_show=False)
    ax.set_xlabel("Years"); ax.set_ylabel("Progression-free interval")
    ax.text(0.05, 0.08, f"log-rank P = {lr.p_value:.3g}", transform=ax.transAxes)
    fig.tight_layout(); fig.savefig(OUT / f"fig_KM_{outdir_tag}.png", dpi=300); plt.close(fig)
    return dict(n=len(d), events=int(d.PFI.sum()), covariates="+".join(cov), HR_adj=s["exp(coef)"], lo=s["exp(coef) lower 95%"],
                hi=s["exp(coef) upper 95%"], p_adj=s["p"], HR_uni=cu["exp(coef)"], p_uni=cu["p"],
                logrank_p=lr.p_value)


def tumor_vs_normal(gimp, bnp):
    f = list(C.TCGA_DIR.glob("*star_counts*"))
    if not f:
        return None
    try:
        import mygene
        x = pd.read_csv(f[0], sep="\t", index_col=0)
        ens = x.index.str.replace(r"\.\d+$", "", regex=True)
        q = mygene.MyGeneInfo().querymany(list(set(ens)), scopes="ensembl.gene", fields="symbol",
                                          species="human", as_dataframe=True, verbose=False)
        sym = q["symbol"].dropna()
        sym = sym[~sym.index.duplicated()]
        x.index = ens.map(sym)
        x = x[x.index.notna()].groupby(level=0).mean()
        x = x[(x > 1).mean(1) >= 0.2]
        z = x.sub(x.mean(1), axis=0).div(x.std(1) + 1e-9, axis=0)
        typ = x.columns.str[13:15]
        pid = x.columns.str[:12]
        res = {}
        for nm, gs in [("GIMP", gimp), ("BNP", bnp)]:
            sc_, n = score(z, gs)
            tumor = sc_[typ == "01"].groupby(pid[typ == "01"]).mean()
            normal = sc_[typ == "11"].groupby(pid[typ == "11"]).mean()
            pair = pd.concat([tumor, normal], axis=1, keys=["T", "N"]).dropna()
            res[nm] = dict(n_genes=n, n_pairs=len(pair), mean_T=pair["T"].mean(), mean_N=pair["N"].mean(),
                           p_paired_wilcoxon=stats.wilcoxon(pair["T"], pair["N"]).pvalue if len(pair) > 5 else None,
                           p_unpaired=stats.mannwhitneyu(tumor, normal).pvalue)
        return res
    except Exception as e:
        log.warning(f"肿瘤-正常比较跳过：{e}")
        return None


def main():
    log.info("Step07 代码版本 v3（TERT 按检测样本清单编码）")
    gimp = pd.read_csv(C.RES / "05_program" / "GIMP_genes.csv").gene.tolist()
    bnp = pd.read_csv(C.RES / "05_program" / "BNP_genes.csv").gene.tolist()
    clin = load_clinical()
    e = load_expr()
    e = e[[c for c in e.columns if str(c)[13:15] in ("01", "")] or e.columns]
    z = e.sub(e.mean(1), axis=0).div(e.std(1) + 1e-9, axis=0)
    df = pd.DataFrame(index=z.columns)
    for nm, gs in [("GIMP", gimp), ("BNP", bnp)]:
        s, n = score(z, gs)
        df[nm] = (s - s.mean()) / s.std()
        log.info(f"{nm}: TCGA 中可用 {n}/{len(gs)} 个基因")
    df = df.join(clin.set_index("SAMPLE_ID"), how="left")
    cdr = pd.read_csv(next(C.CLIN_DIR.glob("Survival_SupplementalTable*")), sep="\t")
    pcol = find_col(cdr.columns, ["_PATIENT", "bcr_patient_barcode"])
    cdr = cdr.drop_duplicates(pcol).set_index(pcol)[["PFI", "PFI.time", "OS", "OS.time"]]
    df["PATIENT_ID"] = df["PATIENT_ID"].fillna(pd.Series(df.index.str[:12], index=df.index))
    df = df.join(cdr, on="PATIENT_ID")
    df.to_csv(OUT / "tcga_scores_clinical.csv")

    rows = []
    for prog in ["GIMP", "BNP"]:
        for y in ["N1", "ETE", "T34", "stage34", "M1", "risk_high", "TERT_mut", "braf_like"]:
            r = logit_or(df, y, prog); r.update(program=prog, outcome=y, test="logistic(OR per SD, adj age+sex)")
            rows.append(r)
            if y not in ("braf_like",):
                r = logit_or(df, y, prog, cov=("age", "male", "BRS"))
                r.update(program=prog, outcome=y, test="logistic(OR per SD, adj age+sex+BRAF-RAS score)")
                rows.append(r)
                r = logit_or(df, y, prog, cov=("age", "male", "BRS", "TDS"))
                r.update(program=prog, outcome=y, test="logistic(OR per SD, adj age+sex+BRS+TDS)")
                rows.append(r)
        for y in ["BRS", "TDS", "ERK", "MACIS"]:
            d = df[[prog, y]].dropna()
            if len(d) > 30:
                rho, p = stats.spearmanr(d[prog], d[y])
                rows.append(dict(program=prog, outcome=y, test="spearman", n=len(d), rho=rho, p=p))
        try:
            r = survival(df, prog, prog); r.update(program=prog, outcome="PFI", test="Cox (HR per SD)")
            rows.append(r)
        except Exception as ex:
            log.warning(f"{prog} 生存分析失败：{ex}")
        hd = df[[prog, "hist"]].dropna()
        hd = hd[hd["hist"].isin(hd["hist"].value_counts()[lambda v: v >= 10].index) & (hd["hist"] != "nan")]
        if hd["hist"].nunique() >= 2:
            p = stats.kruskal(*[g[prog] for _, g in hd.groupby("hist")]).pvalue
            means = hd.groupby("hist")[prog].mean().round(3).to_dict()
            rows.append(dict(program=prog, outcome="histological type", test="Kruskal-Wallis", n=len(hd), p=p, note=str(means)))
    assoc = pd.DataFrame(rows)
    assoc.to_csv(OUT / "tcga_associations.csv", index=False)
    log.info("\n" + assoc.to_string())

    # 随机基因集零分布（淋巴结转移 & PFI）
    rng = np.random.default_rng(C.SEED)
    k = sum(g in z.index for g in gimp)
    null_or = []
    obs = assoc[(assoc.program == "GIMP") & (assoc.outcome == "N1")].OR.iloc[0]
    for i in range(C.N_RANDOM_SETS):
        gs = rng.choice(z.index, k, replace=False)
        s = z.loc[gs].mean(0); df["_rnd"] = ((s - s.mean()) / s.std()).reindex(df.index)
        null_or.append(logit_or(df, "N1", "_rnd")["OR"])
    null_or = np.array(null_or)
    emp_p = (np.sum(null_or >= obs) + 1) / (len(null_or) + 1)
    save_json({"GIMP_OR_N1": obs, "null_median_OR": float(np.nanmedian(null_or)), "empirical_p": emp_p,
               "n_random_sets": len(null_or)}, OUT / "random_geneset_null_N1.json")
    log.info(f"随机基因集检验：GIMP OR={obs:.3f}，经验P={emp_p:.4f}")

    tn = tumor_vs_normal(gimp, bnp)
    if tn:
        save_json(tn, OUT / "tumor_vs_normal.json")

    fig, ax = plt.subplots(1, 4, figsize=(14, 3.6))
    for a_, y in zip(ax[:2], ["N1", "ETE"]):
        d = df[[y, "GIMP"]].dropna()
        sns.boxplot(data=d, x=y, y="GIMP", ax=a_, color="lightsteelblue", showfliers=False)
        sns.stripplot(data=d, x=y, y="GIMP", ax=a_, size=2, color="k", alpha=.4)
    for a_, y in zip(ax[2:], ["BRS", "ERK"] if df["BRS"].notna().sum() > 30 else ["TDS", "ERK"]):
        d = df[[y, "GIMP"]].dropna()
        a_.scatter(d[y], d["GIMP"], s=5, alpha=.5); a_.set_xlabel(y); a_.set_ylabel("GIMP")
    fig.tight_layout(); fig.savefig(OUT / "fig_tcga_GIMP.png", dpi=300); plt.close(fig)
    log.info("Step07 完成")


if __name__ == "__main__":
    main()
