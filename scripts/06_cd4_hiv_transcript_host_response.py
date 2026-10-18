from _bootstrap import finish_figure, np, pd, plt, setup, write_table
import seaborn as sns


def plot_fig4b(counts, out):
    values = pd.to_numeric(counts["HIV_positive_cells"], errors="coerce").fillna(0)
    fig, ax = plt.subplots()
    ax.pie(values, labels=counts["Subtype"].astype(str), autopct=lambda p: f"{p:.1f}%" if p >= 4 else "", startangle=90, colors=sns.color_palette("Set2", len(values)), wedgeprops={"linewidth": 0.7, "edgecolor": "white"}, textprops={"fontsize": 7})
    ax.set_title("HIV RNA-positive CD4 T cells")
    finish_figure(fig, out / "Fig4b_hiv_positive_cd4_state_composition", 5.5, 5)


def plot_fig4c(genes, out):
    states = genes.groupby("Subtype").size().sort_values(ascending=False).index
    rows = int(np.ceil(len(states) / 3))
    fig, axes = plt.subplots(rows, 3, squeeze=False)
    for ax, state in zip(axes.flat, states):
        frame = genes[genes["Subtype"].eq(state)]
        colors = frame["up_down"].astype(str).str.lower().map({"up": "#D1495B", "down": "#4C78A8"}).fillna("#A0A0A0")
        ax.scatter(frame["avg_log2FC"], frame["minus_log10_p"], c=colors, s=9, alpha=0.75)
        label = frame.nlargest(3, "minus_log10_p")
        for row in label.itertuples(index=False):
            ax.text(row.avg_log2FC, row.minus_log10_p, str(row.Gene), fontsize=5)
        ax.axvline(0, color="#808080", linewidth=0.6)
        ax.set_title(str(state), fontsize=8)
        ax.set(xlabel="Average log2 fold change", ylabel="−log10 adjusted P")
    for ax in axes.flat[len(states):]:
        ax.axis("off")
    finish_figure(fig, out / "Fig4c_hiv_positive_cd4_differential_expression", 10, max(5, rows * 3.2))


def plot_edfig5_hiv_counts(counts, out):
    data = counts.copy()
    data["HIV_positive_cells"] = pd.to_numeric(data["HIV_positive_cells"], errors="coerce")
    data["Percent_of_512"] = pd.to_numeric(data["Percent_of_512"], errors="coerce")
    fig, axes = plt.subplots(1, 2)
    sns.barplot(data=data, y="Subtype", x="HIV_positive_cells", hue="Subtype", legend=False, palette="Reds", ax=axes[0])
    sns.barplot(data=data, y="Subtype", x="Percent_of_512", hue="Subtype", legend=False, palette="Reds", ax=axes[1])
    axes[0].set(xlabel="HIV-positive cells", ylabel="CD4 state")
    axes[1].set(xlabel="HIV-positive cells (%)", ylabel="")
    finish_figure(fig, out / "EDFig5e_hiv_positive_cd4_summary", 8, 4.5)


def main():
    _, repo, out = setup("HIV-positive CD4 T-cell response", "06_cd4_hiv_response")
    counts = repo.table(15, "ST15a")
    genes = repo.table(15, "ST15d")
    write_table(counts, out / "hiv_positive_cd4_counts.tsv")
    write_table(genes, out / "hiv_positive_cd4_differential_genes.tsv")
    genes["avg_log2FC"] = pd.to_numeric(genes["avg_log2FC"], errors="coerce")
    genes["Bonferroni_adjusted_P"] = pd.to_numeric(genes["Bonferroni_adjusted_P"], errors="coerce")
    genes["minus_log10_p"] = -np.log10(genes["Bonferroni_adjusted_P"].clip(lower=np.finfo(float).tiny))
    plot_fig4b(counts, out)
    plot_fig4c(genes, out)
    plot_edfig5_hiv_counts(counts, out)


if __name__ == "__main__":
    main()
