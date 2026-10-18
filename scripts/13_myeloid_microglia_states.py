from _bootstrap import donor_composition, finish_figure, heatmap, pd, plt, setup, write_table
import seaborn as sns


def main():
    _, repo, out = setup("Myeloid and microglia-like states", "13_myeloid_states")
    cells = repo.source(1, usecols=["orig.ident", "Disease group", "Type", "SubType_260501"])
    cells = cells[cells["Type"].isin(["Monocyte", "Macrophage"])]
    composition = donor_composition(cells, "orig.ident", "Disease group", "SubType_260501")
    markers = repo.table(27)
    markers["avg_log2FC"] = pd.to_numeric(markers["avg_log2FC"], errors="coerce")
    top = markers.sort_values(["cell_state", "avg_log2FC"], ascending=[True, False]).groupby("cell_state", observed=True).head(5)
    matrix = top.pivot_table(index="cell_state", columns="gene_symbol", values="avg_log2FC", aggfunc="max", fill_value=0)
    write_table(composition, out / "myeloid_donor_composition.tsv")
    write_table(top, out / "myeloid_top_markers.tsv")
    top["pct.1"] = pd.to_numeric(top["pct.1"], errors="coerce")
    fig, ax = plt.subplots()
    sns.scatterplot(data=top, x="gene_symbol", y="cell_state", size="pct.1", hue="avg_log2FC", sizes=(20, 220), palette="vlag", edgecolor="#404040", linewidth=0.3, ax=ax)
    ax.set(xlabel="Marker gene", ylabel="Myeloid state")
    ax.tick_params(axis="x", rotation=55)
    ax.legend(frameon=False, fontsize=6, bbox_to_anchor=(1.02, 1), loc="upper left")
    finish_figure(fig, out / "Fig6e_myeloid_marker_dotplot", 11, 4.8)
    composition["percentage"] = 100 * composition["proportion"]
    states = composition.groupby("SubType_260501")["cells"].sum().nlargest(8).index
    grid = sns.catplot(data=composition[composition["SubType_260501"].isin(states)], x="Disease group", y="percentage", col="SubType_260501", col_wrap=4, kind="box", sharey=False, height=2.6, fliersize=0)
    for ax in grid.axes.flat:
        ax.tick_params(axis="x", rotation=55, labelsize=5)
    grid.figure.savefig(out / "EDFig7a_myeloid_state_composition.pdf", bbox_inches="tight")
    grid.figure.savefig(out / "EDFig7a_myeloid_state_composition.png", dpi=300, bbox_inches="tight")
    plt.close(grid.figure)


if __name__ == "__main__":
    main()
