from _bootstrap import finish_figure, pd, plt, setup, write_table
import seaborn as sns


def plot_model_path(data, label_column, selected, path):
    frame = data[data[label_column].eq(selected)].copy()
    frame["Correlation coefficient (rho)"] = pd.to_numeric(frame["Correlation coefficient (rho)"], errors="coerce")
    frame["BH-adjusted q value"] = pd.to_numeric(frame["BH-adjusted q value"], errors="coerce")
    frame = frame.dropna(subset=["Correlation coefficient (rho)", "Statistical model"])
    frame["Statistical model"] = frame["Statistical model"].astype(str)
    fig, ax = plt.subplots()
    ax.plot(frame["Statistical model"], frame["Correlation coefficient (rho)"], marker="o", color="#4C78A8")
    for i, row in frame.reset_index(drop=True).iterrows():
        ax.text(i, row["Correlation coefficient (rho)"] + 0.02, f"q={row['BH-adjusted q value']:.3g}", ha="center", fontsize=6)
    ax.axhline(0, color="#808080", linewidth=0.7)
    ax.tick_params(axis="x", rotation=25)
    ax.set(xlabel="Adjustment model", ylabel="Partial Spearman rho", title=selected)
    finish_figure(fig, path, 7, 4.5)


def main():
    _, repo, out = setup("BCR-proteome integration", "17_bcr_proteome")
    cohort = repo.table(36, "ST36a")
    models = repo.table(36, "ST36d")
    program = repo.table(35, "ST35b")
    v_gene = repo.table(35, "ST35c")
    write_table(cohort, out / "matched_bcr_proteome_cohort.tsv")
    write_table(models, out / "bcr_proteome_adjusted_models.tsv")
    combined = pd.concat([program.assign(analysis="BCR program"), v_gene.assign(analysis="IGHV usage")], ignore_index=True)
    combined["Correlation coefficient (rho)"] = pd.to_numeric(combined["Correlation coefficient (rho)"], errors="coerce")
    combined["BH-adjusted q value"] = pd.to_numeric(combined["BH-adjusted q value"], errors="coerce")
    write_table(combined, out / "bcr_proteome_associations.tsv")
    combined = combined.sort_values("Correlation coefficient (rho)")
    labels = combined.iloc[:, 0].astype(str) + " | " + combined["Statistical model"].astype(str)
    colors = combined["BH-adjusted q value"].lt(0.05).map({True: "#E45756", False: "#A0A0A0"})
    fig, ax = plt.subplots()
    ax.scatter(combined["Correlation coefficient (rho)"], range(len(combined)), c=colors, s=24)
    ax.axvline(0, color="black", linewidth=0.7)
    ax.set_yticks(range(len(combined)), labels, fontsize=6)
    ax.set(xlabel="Spearman rho", ylabel="Association model")
    finish_figure(fig, out / "bcr_proteome_correlations", 8, 9)
    plot_model_path(program, "Protein program", "Antibody output", out / "EDFig8a_bcr_expansion_antibody_output")
    plot_model_path(v_gene, "V gene", "IGHV6-1", out / "Fig7e_ighv6_1_model_summary")
    model_frame = models[models["BCR_metric"].eq("expanded_fraction")].copy()
    metric = "rho_full_adjusted" if "rho_full_adjusted" in model_frame else "rho_unadjusted"
    model_frame[metric] = pd.to_numeric(model_frame[metric], errors="coerce")
    fig, ax = plt.subplots()
    sns.barplot(data=model_frame, x=metric, y="Module", hue="Module", palette="vlag", legend=False, ax=ax)
    ax.axvline(0, color="#606060", linewidth=0.7)
    ax.set(xlabel="Adjusted Spearman rho", ylabel="Proteomic module")
    finish_figure(fig, out / "Fig7f_bcr_expansion_nmf_module_associations", 6.5, 5)


if __name__ == "__main__":
    main()
