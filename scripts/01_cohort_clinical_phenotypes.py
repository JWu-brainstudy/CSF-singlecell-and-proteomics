from _bootstrap import PALETTE, finish_figure, np, pd, plt, setup, write_table
import seaborn as sns


GROUPS = ["HIV ", "HIV _NS", "HIV _TBM", "HIV _CM", "HIV _TE", "HIV _CMV"]


def numeric_prefix(value):
    return pd.to_numeric(str(value).replace(",", "").split(" ")[0].split("(")[0], errors="coerce")


def plot_fig1b(laboratory, symptoms, out):
    available = [column for column in GROUPS if column in laboratory.columns]
    lab = laboratory.set_index("Variable")[available].map(numeric_prefix)
    lab = lab.loc[lab.notna().sum(axis=1).ge(3)]
    lab = lab.sub(lab.mean(axis=1), axis=0).div(lab.std(axis=1).replace(0, np.nan), axis=0)
    symptom_groups = [column for column in symptoms.columns if str(column).startswith("HIV")]
    symptom = symptoms.set_index("Symptom")[symptom_groups].map(numeric_prefix)
    symptom = symptom.loc[symptom.notna().sum(axis=1).ge(3)]
    symptom = symptom.div(symptom.max(axis=1).replace(0, np.nan), axis=0)
    fig, axes = plt.subplots(1, 2, gridspec_kw={"width_ratios": [1.2, 1]})
    sns.heatmap(lab, cmap="vlag", center=0, yticklabels=True, cbar_kws={"label": "Row z score"}, ax=axes[0])
    sns.heatmap(symptom, cmap="YlOrRd", vmin=0, vmax=1, yticklabels=True, cbar_kws={"label": "Relative frequency"}, ax=axes[1])
    axes[0].set(xlabel="Clinical group", ylabel="Clinical laboratory variable", title="Laboratory profile")
    axes[1].set(xlabel="Clinical group", ylabel="Neurological symptom", title="Symptom profile")
    for ax in axes:
        ax.tick_params(axis="x", rotation=45, labelsize=6)
        ax.tick_params(axis="y", labelsize=5)
    finish_figure(fig, out / "Fig1b_clinical_phenotype_overview", 12, 8)


def plot_cohort_flow(cohort, out):
    numeric = [column for column in cohort.columns if column not in ["Infection tier", "Disease group"]]
    long = cohort.melt(id_vars=["Disease group"], value_vars=numeric, var_name="Dataset", value_name="Participants")
    long["Participants"] = pd.to_numeric(long["Participants"], errors="coerce")
    fig, ax = plt.subplots()
    sns.barplot(data=long, x="Disease group", y="Participants", hue="Dataset", palette="Blues", ax=ax)
    ax.tick_params(axis="x", rotation=40)
    ax.legend(frameon=False, fontsize=6, bbox_to_anchor=(1.02, 1), loc="upper left")
    ax.set(xlabel="Clinical group", ylabel="Participants")
    finish_figure(fig, out / "EDFig1a_dataset_availability", 8, 4.8)


def main():
    _, repo, out = setup("Cohort and clinical phenotype summaries", "01_cohort")
    cohort, laboratory, symptoms, metadata = repo.table(1), repo.table(2), repo.table(3), repo.table(4)
    group = "Disease group"
    counts = metadata.groupby(group, observed=True).size().rename("participants").reset_index()
    write_table(counts, out / "cohort_counts.tsv")
    for number, table in enumerate([cohort, laboratory, symptoms, metadata], start=1):
        write_table(table, out / f"ST{number:02d}_analysis_table.tsv")
    plot_fig1b(laboratory, symptoms, out)
    plot_cohort_flow(cohort, out)
    fig, ax = plt.subplots()
    ax.bar(counts[group].astype(str), counts["participants"], color=PALETTE[: len(counts)])
    ax.set(xlabel="Clinical group", ylabel="Participants")
    ax.tick_params(axis="x", rotation=35)
    finish_figure(fig, out / "EDFig1b_cohort_composition")


if __name__ == "__main__":
    main()
