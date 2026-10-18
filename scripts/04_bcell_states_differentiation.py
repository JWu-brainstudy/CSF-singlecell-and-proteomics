from _bootstrap import PALETTE, donor_composition, finish_figure, kruskal_by_feature, mann_whitney_by_feature, pd, plt, setup, write_table
import seaborn as sns


GROUP_ORDER = ["HIV", "HIV_TBM", "HIV_CMV", "HIV_NS", "HIV_TE", "HIV_CM"]


def plot_fig3b(trajectory, out):
    value = "Median_CytoTRACE_relative_score"
    data = trajectory.sort_values(value)
    fig, ax = plt.subplots()
    sns.barplot(data=data, x=value, y="Fine cell subtype", hue="Fine cell subtype", palette="viridis", legend=False, ax=ax)
    ax.set(xlabel="Median relative CytoTRACE score", ylabel="B-cell state")
    finish_figure(fig, out / "Fig3b_bcell_cytotrace_order", 7, 5.5)


def plot_fig3c(composition, contrasts, out):
    effect = "HIV_NS median proportion (%)"
    contrasts[effect] = pd.to_numeric(contrasts[effect], errors="coerce")
    selected = contrasts.nlargest(4, effect)["Fine cell subtype"].astype(str)
    data = composition[composition["SubType_260501"].isin(selected)].copy()
    data["percentage"] = 100 * data["proportion"]
    grid = sns.catplot(
        data=data, x="Disease group", y="percentage", col="SubType_260501", col_order=selected,
        order=GROUP_ORDER, hue="Disease group", hue_order=GROUP_ORDER, palette=PALETTE,
        kind="box", sharey=False, height=3.1, aspect=0.9, legend=False, fliersize=0,
    )
    for ax, state in zip(grid.axes.flat, selected):
        sns.stripplot(data=data[data["SubType_260501"].eq(state)], x="Disease group", y="percentage", order=GROUP_ORDER, color="#303030", size=2.5, alpha=0.65, ax=ax)
        row = contrasts[contrasts["Fine cell subtype"].astype(str).eq(state)]
        if len(row):
            ax.text(0.03, 0.97, f"MW P={pd.to_numeric(row.iloc[0]['Two-sided Mann-Whitney P'], errors='coerce'):.3g}\nBH q={pd.to_numeric(row.iloc[0]['BH q (11 tests)'], errors='coerce'):.3g}", transform=ax.transAxes, va="top", fontsize=6)
        ax.tick_params(axis="x", rotation=55, labelsize=5)
        ax.set_title(state, fontsize=8)
    grid.set_axis_labels("", "Cells per donor (%)")
    grid.figure.savefig(out / "Fig3c_bcell_state_composition.pdf", bbox_inches="tight")
    grid.figure.savefig(out / "Fig3c_bcell_state_composition.png", dpi=300, bbox_inches="tight")
    plt.close(grid.figure)
    all_states = contrasts["Fine cell subtype"].astype(str).tolist()
    full = composition[composition["SubType_260501"].isin(all_states)].copy()
    full["percentage"] = 100 * full["proportion"]
    grid = sns.catplot(data=full, x="Disease group", y="percentage", col="SubType_260501", col_order=all_states, col_wrap=4, order=GROUP_ORDER, kind="box", sharey=False, height=2.6, fliersize=0)
    for ax in grid.axes.flat:
        ax.tick_params(axis="x", rotation=55, labelsize=5)
    grid.figure.savefig(out / "EDFig3b_bcell_state_composition.pdf", bbox_inches="tight")
    grid.figure.savefig(out / "EDFig3b_bcell_state_composition.png", dpi=300, bbox_inches="tight")
    plt.close(grid.figure)


def plot_edfig3a(markers, out):
    for column in ["avg_log2FC", "pct.1", "marker_rank"]:
        markers[column] = pd.to_numeric(markers[column], errors="coerce")
    data = markers.sort_values(["cell_type", "marker_rank"]).groupby("cell_type", observed=True).head(4)
    fig, ax = plt.subplots()
    sns.scatterplot(data=data, x="gene_symbol", y="cell_type", size="pct.1", hue="avg_log2FC", sizes=(20, 220), palette="vlag", edgecolor="#404040", linewidth=0.3, ax=ax)
    ax.tick_params(axis="x", rotation=55)
    ax.set(xlabel="B/plasma marker gene", ylabel="B/plasma state")
    ax.legend(frameon=False, fontsize=6, bbox_to_anchor=(1.02, 1), loc="upper left")
    finish_figure(fig, out / "EDFig3a_bcell_marker_dotplot", 11, 5)


def main():
    _, repo, out = setup("B-cell states and differentiation", "04_bcell_states")
    trajectory = repo.table(11, "ST11a")
    global_tests = repo.table(11, "ST11b")
    contrasts = repo.table(11, "ST11c")
    write_table(global_tests, out / "bcell_global_tests.tsv")
    write_table(contrasts, out / "bcell_hiv_ns_contrasts.tsv")
    value = "Median_CytoTRACE_relative_score"
    trajectory[value] = pd.to_numeric(trajectory[value], errors="coerce")
    cells = repo.source(1, usecols=["orig.ident", "Disease group", "Type", "SubType_260501"])
    cells = cells[cells["Type"].isin(["B", "Bcell", "B_cell", "Plasma", "Plasmablast"])]
    composition = donor_composition(cells, "orig.ident", "Disease group", "SubType_260501")
    global_recomputed = kruskal_by_feature(composition, "SubType_260501", "Disease group", "proportion")
    binary = composition.assign(contrast_group=composition["Disease group"].map(lambda value: "HIV_NS" if value == "HIV_NS" else "Other groups"))
    contrast_recomputed = mann_whitney_by_feature(binary, "SubType_260501", "contrast_group", "proportion", "Other groups")
    write_table(composition, out / "bcell_donor_composition.tsv")
    write_table(global_recomputed, out / "ST11b_bcell_kruskal_recomputed.tsv")
    write_table(contrast_recomputed, out / "ST11c_bcell_hiv_ns_mann_whitney_recomputed.tsv")
    plot_fig3b(trajectory, out)
    plot_fig3c(composition, contrasts, out)
    plot_edfig3a(repo.table(10), out)


if __name__ == "__main__":
    main()
