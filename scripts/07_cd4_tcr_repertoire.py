from _bootstrap import PALETTE, donor_composition, finish_figure, kruskal_by_feature, mann_whitney_by_feature, pd, plt, setup, write_table
import networkx as nx
import seaborn as sns


def stacked(data, index, column, value, path, ylabel, cmap):
    frame = data.copy()
    frame[value] = pd.to_numeric(frame[value], errors="coerce")
    pivot = frame.pivot_table(index=index, columns=column, values=value, aggfunc="sum", fill_value=0)
    fig, ax = plt.subplots()
    pivot.plot.bar(stacked=True, colormap=cmap, width=0.82, ax=ax)
    ax.set(xlabel="Clinical group", ylabel=ylabel)
    ax.tick_params(axis="x", rotation=40)
    ax.legend(frameon=False, fontsize=6, bbox_to_anchor=(1.02, 1), loc="upper left")
    finish_figure(fig, path, 8, 5)


def plot_fig4g(edges, out):
    graph = nx.DiGraph()
    for row in edges.itertuples(index=False):
        graph.add_edge(str(row.source), str(row.target), weight=float(pd.to_numeric(row.n, errors="coerce") or 0))
    layers = {}
    for node in graph.nodes:
        if node.startswith("CD4"):
            layers[node] = 0
        elif node.startswith("TRBV"):
            layers[node] = 1
        elif node.startswith("TRAV"):
            layers[node] = 2
        else:
            layers[node] = 3
    nx.set_node_attributes(graph, layers, "layer")
    position = nx.multipartite_layout(graph, subset_key="layer", align="vertical")
    weights = [graph[u][v]["weight"] for u, v in graph.edges]
    scale = max(weights, default=1)
    fig, ax = plt.subplots()
    nx.draw_networkx_edges(graph, position, width=[0.3 + 4 * w / scale for w in weights], edge_color="#7F8C8D", alpha=0.35, arrows=False, ax=ax)
    nx.draw_networkx_nodes(graph, position, node_size=120, node_color="#72B7B2", edgecolors="white", ax=ax)
    nx.draw_networkx_labels(graph, position, font_size=5, ax=ax)
    ax.axis("off")
    finish_figure(fig, out / "Fig4g_cd4_tcr_alluvial_network", 10, 6)


def plot_edfig5a(markers, out):
    for column in ["avg_log2FC", "pct.1", "marker_rank"]:
        markers[column] = pd.to_numeric(markers[column], errors="coerce")
    data = markers.sort_values(["cell_type", "marker_rank"]).groupby("cell_type", observed=True).head(4)
    fig, ax = plt.subplots()
    sns.scatterplot(data=data, x="gene_symbol", y="cell_type", size="pct.1", hue="avg_log2FC", sizes=(20, 220), palette="vlag", edgecolor="#404040", linewidth=0.3, ax=ax)
    ax.tick_params(axis="x", rotation=55)
    ax.set(xlabel="CD4 marker gene", ylabel="CD4 state")
    ax.legend(frameon=False, fontsize=6, bbox_to_anchor=(1.02, 1), loc="upper left")
    finish_figure(fig, out / "EDFig5a_cd4_marker_dotplot", 9, 4.5)


def plot_edfig5b(repo, out):
    cells = repo.source(1, usecols=["orig.ident", "Disease group", "Type", "SubType_260501"])
    cells = cells[cells["Type"].eq("CD4T")]
    data = donor_composition(cells, "orig.ident", "Disease group", "SubType_260501")
    tests = kruskal_by_feature(data, "SubType_260501", "Disease group", "proportion")
    write_table(data, out / "cd4_donor_state_composition.tsv")
    write_table(tests, out / "ST14c_cd4_state_kruskal_recomputed.tsv")
    data["percentage"] = 100 * data["proportion"]
    states = data.groupby("SubType_260501")["cells"].sum().nlargest(5).index
    order = ["HIV", "HIV_TBM", "HIV_CMV", "HIV_NS", "HIV_TE", "HIV_CM"]
    grid = sns.catplot(data=data[data["SubType_260501"].isin(states)], x="Disease group", y="percentage", col="SubType_260501", col_order=states, kind="box", order=order, hue="Disease group", hue_order=order, palette=dict(zip(order, PALETTE)), legend=False, sharey=False, height=2.8, fliersize=0)
    for ax in grid.axes.flat:
        ax.tick_params(axis="x", rotation=55, labelsize=5)
    grid.figure.savefig(out / "EDFig5b_cd4_state_composition.pdf", bbox_inches="tight")
    grid.figure.savefig(out / "EDFig5b_cd4_state_composition.png", dpi=300, bbox_inches="tight")
    plt.close(grid.figure)


def main():
    _, repo, out = setup("CD4 T-cell receptor repertoire", "07_cd4_tcr")
    clone = repo.table(16, "ST16a")
    donor = repo.table(16, "ST16b")
    subtype = repo.table(17, "ST17a")
    expanded_subtype = repo.table(17, "ST17b")
    receptor = repo.source(3, usecols=["SubType_260501", "IR_VDJ_1_v_call", "IR_VJ_1_v_call"])
    receptor = receptor[receptor["SubType_260501"].astype(str).str.startswith("CD4")].dropna()
    first = receptor.groupby(["SubType_260501", "IR_VDJ_1_v_call"], observed=True).size().rename("n").reset_index().rename(columns={"SubType_260501": "source", "IR_VDJ_1_v_call": "target"})
    second = receptor.groupby(["IR_VDJ_1_v_call", "IR_VJ_1_v_call"], observed=True).size().rename("n").reset_index().rename(columns={"IR_VDJ_1_v_call": "source", "IR_VJ_1_v_call": "target"})
    edges = pd.concat([first, second], ignore_index=True)
    edges = edges.sort_values("n", ascending=False).groupby(edges["source"].map(lambda value: 0 if str(value).startswith("CD4") else 1), group_keys=False).head(80)
    write_table(donor, out / "cd4_tcr_donor_expansion.tsv")
    write_table(subtype, out / "cd4_tcr_subtype_composition.tsv")
    write_table(clone, out / "cd4_tcr_clone_size_distribution.tsv")
    write_table(expanded_subtype, out / "cd4_expanded_tcr_subtypes.tsv")
    write_table(edges, out / "cd4_tcr_alluvial_edges.tsv")
    donor["expanded_fraction"] = pd.to_numeric(donor["Expanded proportion"], errors="coerce")
    donor["analysis"] = "CD4 TCR expansion"
    write_table(kruskal_by_feature(donor, "analysis", "Clinical group", "expanded_fraction"), out / "ST16b_cd4_tcr_expansion_kruskal_recomputed.tsv")
    write_table(mann_whitney_by_feature(donor, "analysis", "Clinical group", "expanded_fraction", "HIV"), out / "ST16b_cd4_tcr_expansion_mann_whitney_recomputed.tsv")
    groups = list(donor["Clinical group"].dropna().unique())
    fig, ax = plt.subplots()
    ax.boxplot([donor.loc[donor["Clinical group"] == group, "expanded_fraction"] for group in groups], tick_labels=groups, showfliers=False)
    ax.set(xlabel="Clinical group", ylabel="Expanded CD4 TCR fraction")
    ax.tick_params(axis="x", rotation=35)
    finish_figure(fig, out / "EDFig5c_cd4_tcr_donor_expansion")
    stacked(clone, "Clinical group", "Clone-size category", "Proportion", out / "Fig4e_cd4_tcr_clone_size", "Fraction of CD4 TCR cells", "Blues")
    stacked(subtype, "Clinical group", "CD4 T-cell subtype", "Proportion", out / "Fig4f_cd4_tcr_subtype_distribution", "Fraction of CD4 TCR cells", "tab20")
    plot_fig4g(edges, out)
    plot_edfig5a(repo.table(14, "ST14a"), out)
    plot_edfig5b(repo, out)


if __name__ == "__main__":
    main()
