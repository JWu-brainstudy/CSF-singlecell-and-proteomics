from _bootstrap import benjamini_hochberg, finish_figure, np, partial_spearman_table, pd, plt, setup, welch_differential, write_table
import seaborn as sns
from scipy.stats import fisher_exact, kruskal


def sample_key(value):
    value = str(value).strip()
    return value[:-2] if value.endswith(".0") else value


def recompute_cime(repo):
    categorical = repo.table(9, "ST09a").copy()
    count_columns = ["a: group and program", "b: group other programs", "c: other groups and program", "d: other groups other programs"]
    for column in count_columns:
        categorical[column] = pd.to_numeric(categorical[column], errors="coerce")
    results = []
    for _, row in categorical.dropna(subset=count_columns).iterrows():
        table = [[row[count_columns[0]], row[count_columns[1]]], [row[count_columns[2]], row[count_columns[3]]]]
        odds_ratio, p_value = fisher_exact(table, alternative="two-sided")
        results.append({"Disease group": row["Disease group"], "Program": row["Program"], "odds_ratio_recomputed": odds_ratio, "p_value_recomputed": p_value})
    fisher = pd.DataFrame(results)
    fisher["q_value_recomputed"] = benjamini_hochberg(fisher["p_value_recomputed"])
    published = categorical[["Disease group", "Program", "Two-sided Fisher P", "BH q (48 tests)"]]
    fisher = fisher.merge(published, on=["Disease group", "Program"], how="left")
    usage = repo.table(9, "ST09c").copy()
    donor_groups = repo.table(7, "ST07a")[["Donor ID", "Clinical group"]].drop_duplicates("Donor ID")
    usage = usage.merge(donor_groups, on="Donor ID", how="left")
    programs = [column for column in usage.columns if str(column).startswith("P") and str(column).endswith(" usage")]
    rows = []
    for program in programs:
        usage[program] = pd.to_numeric(usage[program], errors="coerce")
        samples = [frame[program].dropna().to_numpy() for _, frame in usage.groupby("Clinical group", observed=True)]
        samples = [sample for sample in samples if len(sample)]
        statistic, p_value = kruskal(*samples) if len(samples) > 1 else (np.nan, np.nan)
        rows.append({"Program": program.replace(" usage", ""), "kruskal_H_recomputed": statistic, "p_value_recomputed": p_value})
    continuous = pd.DataFrame(rows)
    continuous["q_value_recomputed"] = benjamini_hochberg(continuous["p_value_recomputed"])
    return fisher, continuous


def recompute_proteomics(repo):
    proteins = repo.source(4)
    annotations = ["ProteinAccessions", "Accession", "PG.Genes", "PG.ProteinDescriptions", "PG.ProteinNames"]
    activity = repo.table(33, "ST33c")
    activity["sample_key"] = activity["samplename"].map(sample_key)
    group_map = activity.set_index("sample_key")["group"].astype(str).to_dict()
    sample_columns = [column for column in proteins.columns if column not in annotations and not str(column).startswith("QC")]
    group_map.update({str(column): "Control" for column in sample_columns if str(column).startswith("C")})
    sample_columns = [column for column in sample_columns if str(column) in group_map]
    accession = proteins["Accession"].fillna(proteins["ProteinAccessions"]).astype(str)
    gene = proteins["PG.Genes"].fillna(accession).astype(str).str.split(";").str[0]
    feature = gene + "|" + accession
    duplicate = feature.groupby(feature).cumcount()
    feature = feature + duplicate.map(lambda value: "" if value == 0 else f"#{value + 1}")
    matrix = proteins[sample_columns].apply(pd.to_numeric, errors="coerce").T
    matrix.columns = feature
    matrix = np.log2(matrix.where(matrix.gt(0)))
    groups = pd.Series({column: group_map[str(column)] for column in sample_columns})
    differential = welch_differential(matrix, groups, reference="Control")
    annotation = pd.DataFrame({"feature": feature, "gene": gene, "accession": accession})
    return differential.merge(annotation, on="feature", how="left")


def recompute_bcr_module_correlations(repo):
    cohort = repo.table(36, "ST36a").copy()
    modules = repo.table(36, "ST36b").copy()
    cohort["sample_key"] = cohort["sample_ID"].map(sample_key)
    modules["sample_key"] = modules["sample_ID"].map(sample_key)
    data = cohort.merge(modules.drop(columns=["group", "sample_ID"], errors="ignore"), on="sample_key", how="inner")
    module_columns = [column for column in modules.columns if str(column).startswith("M")]
    metric_columns = [
        "expanded_fraction", "unique_clonotypes", "shannon_clonality", "simpson_concentration",
        "max_clone_fraction", "secretory_fraction", "plasmablast_fraction", "plasma_IGHG4_fraction",
    ]
    metric_columns = [column for column in metric_columns if column in data]
    for column in metric_columns + module_columns:
        data[column] = pd.to_numeric(data[column], errors="coerce")
    pairs = [(metric, module) for metric in metric_columns for module in module_columns]
    unadjusted = partial_spearman_table(data, pairs)
    unadjusted["model"] = "Spearman"
    adjusted = partial_spearman_table(data, pairs, covariates=("group",))
    adjusted["model"] = "Partial Spearman: disease group"
    return data, pd.concat([unadjusted, adjusted], ignore_index=True)


def recompute_multiplicity(repo):
    cd4 = repo.table(15, "ST15d").copy()
    cd4["p_val"] = pd.to_numeric(cd4["p_val"], errors="coerce")
    cd4["bonferroni_within_subtype"] = cd4.groupby("Subtype", observed=True)["p_val"].transform(lambda values: np.minimum(values * values.notna().sum(), 1))
    cd4["bonferroni_all_rows"] = np.minimum(cd4["p_val"] * cd4["p_val"].notna().sum(), 1)
    myeloid = repo.table(29, "ST29b").copy()
    myeloid["P value"] = pd.to_numeric(myeloid["P value"], errors="coerce")
    myeloid["FDR_recomputed_within_group"] = myeloid.groupby("Group", observed=True)["P value"].transform(benjamini_hochberg)
    return cd4, myeloid


def plot_overview(differential, correlations, out):
    counts = differential.assign(significant=differential["q_value"].lt(0.05)).groupby("comparison", observed=True)["significant"].sum().reset_index()
    adjusted = correlations[correlations["model"].str.startswith("Partial")]
    matrix = adjusted.pivot(index="measure_1", columns="measure_2", values="rho")
    fig, axes = plt.subplots(1, 2, gridspec_kw={"width_ratios": [1, 1.5]})
    sns.barplot(data=counts, x="comparison", y="significant", hue="comparison", palette="Set2", legend=False, ax=axes[0])
    axes[0].tick_params(axis="x", rotation=45)
    axes[0].set(xlabel="Group versus control", ylabel="Proteins at BH q < 0.05")
    sns.heatmap(matrix, cmap="vlag", center=0, annot=True, fmt=".2f", cbar_kws={"label": "Group-adjusted rho"}, ax=axes[1])
    axes[1].set(xlabel="Proteomic NMF module", ylabel="BCR metric")
    finish_figure(fig, out / "statistical_reanalysis_overview", 13, 6)


def main():
    _, repo, out = setup("Recompute supplementary-table statistical analyses", "20_statistical_reanalysis")
    cime_fisher, cime_kw = recompute_cime(repo)
    proteomics = recompute_proteomics(repo)
    matched, correlations = recompute_bcr_module_correlations(repo)
    cd4, myeloid = recompute_multiplicity(repo)
    write_table(cime_fisher, out / "ST09a_fisher_recomputed.tsv")
    write_table(cime_kw, out / "ST09b_kruskal_recomputed.tsv")
    write_table(proteomics, out / "ST31_proteomics_welch_differential_recomputed.tsv")
    write_table(proteomics[proteomics["q_value"].lt(0.05)], out / "ST31_proteomics_significant_recomputed.tsv")
    write_table(matched, out / "ST36_matched_patient_analysis_matrix.tsv")
    write_table(correlations, out / "ST36_bcr_nmf_correlations_recomputed.tsv")
    write_table(cd4, out / "ST15d_bonferroni_recomputed.tsv")
    write_table(myeloid, out / "ST29b_fdr_recomputed.tsv")
    plot_overview(proteomics, correlations, out)


if __name__ == "__main__":
    main()
