from _bootstrap import finish_figure, np, pd, plt, setup, write_table
from upsetplot import UpSet, from_contents
import warnings
import seaborn as sns


def plot_shared_sets(differential, direction, path):
    selected = differential[differential["Direction"].eq(direction)]
    contents = {
        comparison: set(frame["Representative gene"].dropna().astype(str))
        for comparison, frame in selected.groupby("Comparison", observed=True)
    }
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=FutureWarning, module="upsetplot")
        data = from_contents(contents)
        figure = plt.figure(figsize=(10, 5.5))
        UpSet(data, subset_size="count", show_counts=True, sort_by="cardinality", min_subset_size=2).plot(fig=figure)
    figure.savefig(path.with_suffix(".pdf"), bbox_inches="tight")
    figure.savefig(path.with_suffix(".png"), dpi=300, bbox_inches="tight")
    plt.close(figure)


def plot_fig7c(enrichment_up, enrichment_down, out):
    frames = []
    for label, table in [("Core up", enrichment_up), ("Core down", enrichment_down)]:
        table = table.copy()
        table["Minus log10 FDR"] = pd.to_numeric(table["Minus log10 FDR"], errors="coerce")
        table["Fold enrichment"] = pd.to_numeric(table["Fold enrichment"], errors="coerce")
        frames.append(table.nlargest(12, "Minus log10 FDR").assign(Direction=label))
    data = pd.concat(frames, ignore_index=True)
    fig, axes = plt.subplots(1, 2)
    for ax, (direction, frame) in zip(axes, data.groupby("Direction", observed=True)):
        frame = frame.sort_values("Minus log10 FDR")
        sns.scatterplot(data=frame, x="Minus log10 FDR", y="Term", size="Fold enrichment", hue="Source", sizes=(30, 220), palette="Set2", ax=ax)
        ax.set(title=direction, xlabel="−log10 FDR", ylabel="Enriched term")
        ax.legend(frameon=False, fontsize=5, bbox_to_anchor=(1.02, 1), loc="upper left")
    finish_figure(fig, out / "Fig7c_shared_core_enrichment", 13, 7)


def main():
    _, repo, out = setup("Differential proteomics and shared core", "15_proteomics_differential")
    differential = repo.table(31)
    core_up = repo.table(32, "ST32a")
    core_down = repo.table(32, "ST32b")
    enrichment_up = repo.table(32, "ST32c")
    enrichment_down = repo.table(32, "ST32d")
    differential["Fold change"] = pd.to_numeric(differential["Fold change"], errors="coerce")
    differential["Adjusted P"] = pd.to_numeric(differential["Adjusted P"], errors="coerce")
    differential["signed_log2_fold_change"] = np.sign(differential["Fold change"]) * np.log2(differential["Fold change"].abs())
    differential["minus_log10_q"] = -np.log10(differential["Adjusted P"].clip(lower=np.finfo(float).tiny))
    write_table(differential, out / "differential_proteins.tsv")
    write_table(pd.concat([core_up.assign(core_direction="up"), core_down.assign(core_direction="down")]), out / "shared_proteomic_core.tsv")
    write_table(pd.concat([enrichment_up.assign(core_direction="up"), enrichment_down.assign(core_direction="down")]), out / "shared_core_enrichment.tsv")
    plot_shared_sets(differential, "Increased", out / "Fig7b_shared_upregulated_proteins")
    plot_shared_sets(differential, "Decreased", out / "Fig7b_shared_downregulated_proteins")
    plot_fig7c(enrichment_up, enrichment_down, out)
    comparisons = list(differential["Comparison"].dropna().unique())
    fig, axes = plt.subplots(2, 3)
    for ax, comparison in zip(axes.flat, comparisons):
        frame = differential[differential["Comparison"] == comparison]
        colors = frame["Direction"].map({"Increased": "#E45756", "Decreased": "#4C78A8"}).fillna("#A0A0A0")
        ax.scatter(frame["signed_log2_fold_change"], frame["minus_log10_q"], c=colors, s=8, alpha=0.7)
        ax.set_title(str(comparison), fontsize=8)
        ax.set(xlabel="log2 fold change", ylabel="−log10 q")
    for ax in axes.flat[len(comparisons):]:
        ax.axis("off")
    finish_figure(fig, out / "differential_proteome", 10, 6.5)


if __name__ == "__main__":
    main()
