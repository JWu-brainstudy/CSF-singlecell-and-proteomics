from _bootstrap import PALETTE, donor_composition, finish_figure, kruskal_by_feature, pd, plt, setup, write_table
import seaborn as sns


GROUP_ORDER = ["HIV", "HIV_TBM", "HIV_CMV", "HIV_NS", "HIV_TE", "HIV_CM"]


def plot_fig1d(cells, out):
    summary = cells[cells["Disease group"].isin(GROUP_ORDER)].groupby(["Disease group", "Type"], observed=True).size().rename("cells").reset_index()
    summary["proportion"] = summary["cells"] / summary.groupby("Disease group")["cells"].transform("sum")
    pivot = summary.pivot(index="Disease group", columns="Type", values="proportion").reindex(GROUP_ORDER).fillna(0)
    order = pivot.sum().sort_values(ascending=False).index
    fig, ax = plt.subplots()
    pivot[order].plot.bar(stacked=True, ax=ax, colormap="tab20", width=0.82)
    ax.set(xlabel="Clinical group", ylabel="Cell fraction", ylim=(0, 1))
    ax.legend(title="Major cell type", bbox_to_anchor=(1.02, 1), loc="upper left", frameon=False, fontsize=6)
    finish_figure(fig, out / "Fig1d_major_cell_composition", 8.5, 5)


def plot_edfig2a(repo, out):
    markers = repo.table(5)
    for column in ["avg_log2FC", "pct.1", "p_val_adj", "marker_rank"]:
        markers[column] = pd.to_numeric(markers[column], errors="coerce")
    selected = markers.sort_values(["cell_type", "marker_rank"]).groupby("cell_type", observed=True).head(3)
    selected["gene_cell"] = selected["gene_symbol"].astype(str)
    fig, ax = plt.subplots()
    sns.scatterplot(data=selected, x="gene_cell", y="cell_type", size="pct.1", hue="avg_log2FC", sizes=(20, 180), palette="coolwarm", edgecolor="black", linewidth=0.25, ax=ax)
    ax.set(xlabel="Top marker genes", ylabel="Major cell type")
    ax.tick_params(axis="x", rotation=60)
    ax.legend(bbox_to_anchor=(1.02, 1), loc="upper left", frameon=False, fontsize=6)
    finish_figure(fig, out / "EDFig2a_major_cell_marker_dotplot", 11, 5.5)


def plot_edfig2b(composition, tests, out):
    order = composition.groupby("Type")["cells"].sum().sort_values(ascending=False).index
    data = composition[composition["Disease group"].isin(GROUP_ORDER)].copy()
    data["percentage"] = 100 * data["proportion"]
    grid = sns.catplot(
        data=data, x="Disease group", y="percentage", col="Type", col_order=order,
        kind="box", col_wrap=5, order=GROUP_ORDER, sharey=False, height=2.3,
        hue="Disease group", hue_order=GROUP_ORDER, legend=False,
        palette=PALETTE, linewidth=0.8, fliersize=0,
    )
    for ax, cell_type in zip(grid.axes.flat, order):
        sns.stripplot(
            data=data[data["Type"].eq(cell_type)], x="Disease group", y="percentage",
            order=GROUP_ORDER, color="#303030", size=2.2, alpha=0.65, ax=ax,
        )
        row = tests[tests["Type"].eq(cell_type)]
        if len(row):
            ax.text(0.02, 0.95, f"KW P={row.iloc[0]['p_value']:.3g}\nBH q={row.iloc[0]['q_value']:.3g}", transform=ax.transAxes, va="top", fontsize=6)
        ax.tick_params(axis="x", rotation=55, labelsize=5)
        ax.set_title(cell_type, fontsize=8)
    grid.set_axis_labels("", "Cells (% of donor CSF cells)")
    grid.figure.savefig(out / "EDFig2b_major_cell_composition.pdf", bbox_inches="tight")
    grid.figure.savefig(out / "EDFig2b_major_cell_composition.png", dpi=300, bbox_inches="tight")
    plt.close(grid.figure)


def main():
    _, repo, out = setup("Major-cell composition and marker analysis", "02_major_cell_composition")
    cells = repo.source(1, usecols=["orig.ident", "Disease group", "Type"])
    composition = donor_composition(cells, "orig.ident", "Disease group", "Type")
    tests = kruskal_by_feature(composition, "Type", "Disease group", "proportion")
    write_table(composition, out / "donor_major_cell_composition.tsv")
    write_table(tests, out / "major_cell_global_tests.tsv")
    plot_fig1d(cells, out)
    plot_edfig2a(repo, out)
    plot_edfig2b(composition, tests, out)


if __name__ == "__main__":
    main()
