from _bootstrap import finish_figure, kruskal_by_feature, mann_whitney_by_feature, np, pd, plt, setup, write_table
import networkx as nx


def stacked(data, index, column, value, path, ylabel, cmap):
    data[value] = pd.to_numeric(data[value], errors="coerce")
    pivot = data.pivot_table(index=index, columns=column, values=value, aggfunc="sum", fill_value=0)
    fig, ax = plt.subplots()
    pivot.plot.bar(stacked=True, ax=ax, colormap=cmap, width=0.82)
    ax.set(xlabel="Clinical group", ylabel=ylabel)
    ax.tick_params(axis="x", rotation=40)
    ax.legend(frameon=False, fontsize=6, bbox_to_anchor=(1.02, 1), loc="upper left")
    finish_figure(fig, path, 8, 5)


def plot_fig5e(edges, out):
    edges = edges.copy()
    edges["n"] = pd.to_numeric(edges["n"], errors="coerce").fillna(0)
    graph = nx.DiGraph()
    for row in edges.nlargest(180, "n").itertuples(index=False):
        graph.add_edge(str(row.source), str(row.target), weight=row.n)
    position = nx.spring_layout(graph, seed=2026, k=1.4)
    weights = [graph[u][v]["weight"] for u, v in graph.edges]
    scale = max(weights, default=1)
    fig, ax = plt.subplots()
    nx.draw_networkx_edges(graph, position, width=[0.2 + 4 * w / scale for w in weights], edge_color="#808080", alpha=0.3, arrows=False, ax=ax)
    nx.draw_networkx_nodes(graph, position, node_size=80, node_color="#FF9D4D", edgecolors="white", ax=ax)
    nx.draw_networkx_labels(graph, position, font_size=4.5, ax=ax)
    ax.axis("off")
    finish_figure(fig, out / "Fig5e_cd8_tcr_alluvial_network", 10, 7)


def main():
    _, repo, out = setup("CD8 TCR expansion and Sankey inputs", "10_cd8_tcr_overview")
    expansion = repo.table(21, "ST21b")
    subtype = repo.table(21, "ST21d")
    donor = repo.table(21, "ST21c")
    edges = repo.table(22, "ST22b")
    write_table(expansion, out / "cd8_tcr_clone_size_distribution.tsv")
    write_table(edges, out / "cd8_tcr_sankey_edges.tsv")
    write_table(subtype, out / "cd8_tcr_subtype_distribution.tsv")
    write_table(donor, out / "cd8_tcr_donor_expansion.tsv")
    stacked(expansion, "Disease.group", "clone_size_category", "pct_within_group", out / "Fig5c_cd8_tcr_clone_size", "Fraction of CD8 TCR cells (%)", "Oranges")
    stacked(subtype, "Disease.group", "seurat_clusters", "pct_within_group", out / "Fig5d_cd8_tcr_subtype_distribution", "Fraction of CD8 TCR cells (%)", "tab20")
    plot_fig5e(edges, out)
    donor["expanded_pct"] = pd.to_numeric(donor["expanded_pct"], errors="coerce")
    donor["analysis"] = "CD8 TCR expansion"
    write_table(kruskal_by_feature(donor, "analysis", "Disease.group", "expanded_pct"), out / "ST21c_cd8_tcr_expansion_kruskal_recomputed.tsv")
    write_table(mann_whitney_by_feature(donor, "analysis", "Disease.group", "expanded_pct", "HIV"), out / "ST21c_cd8_tcr_expansion_mann_whitney_recomputed.tsv")
    groups = donor["Disease.group"].dropna().unique()
    fig, ax = plt.subplots()
    ax.boxplot([donor.loc[donor["Disease.group"].eq(group), "expanded_pct"] for group in groups], tick_labels=groups, showfliers=False)
    ax.set(xlabel="Clinical group", ylabel="Expanded CD8 TCR cells (%)")
    ax.tick_params(axis="x", rotation=40)
    finish_figure(fig, out / "EDFig6b_cd8_tcr_donor_expansion")


if __name__ == "__main__":
    main()
