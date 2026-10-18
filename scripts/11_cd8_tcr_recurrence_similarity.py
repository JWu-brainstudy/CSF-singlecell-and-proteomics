from _bootstrap import PALETTE, finish_figure, pd, plt, setup, write_table
import networkx as nx
import seaborn as sns


GROUP_ORDER = ["HIV", "HIV_TBM", "HIV_CMV", "HIV_NS", "HIV_TE", "HIV_CM"]


def plot_fig5f(clonotypes, out):
    data = clonotypes.head(20).copy()
    data["total_cells"] = pd.to_numeric(data["total_cells"], errors="coerce")
    data["overall_sample_prevalence"] = pd.to_numeric(data["overall_sample_prevalence"], errors="coerce")
    data = data.sort_values("total_cells")
    label = data["clonotype"].astype(str) + " | " + data["TRBV"].astype(str) + " + " + data["TRAV"].astype(str)
    fig, ax = plt.subplots()
    points = ax.scatter(data["total_cells"], range(len(data)), s=80 + 500 * data["overall_sample_prevalence"].fillna(0), c=data["overall_sample_prevalence"], cmap="YlOrRd", edgecolor="white")
    ax.set(yticks=range(len(data)), yticklabels=label, xlabel="Total expanded cells", ylabel="Clonotype")
    fig.colorbar(points, ax=ax, label="Sample prevalence")
    finish_figure(fig, out / "Fig5f_top_expanded_cd8_tcr_clonotypes", 8, 7)


def plot_fig5g(highlighted, out):
    data = highlighted.copy()
    data["sample_prevalence"] = pd.to_numeric(data["sample_prevalence"], errors="coerce")
    data["expanded_fraction_within_group"] = pd.to_numeric(data["expanded_fraction_within_group"], errors="coerce")
    data["expanded_cells"] = pd.to_numeric(data["expanded_cells"], errors="coerce").fillna(0)
    data["group"] = data["group"].astype(str)
    top = data.groupby("clonotype")["expanded_cells"].sum().nlargest(5).index
    data = data[data["clonotype"].isin(top)]
    fig, ax = plt.subplots()
    sns.scatterplot(
        data=data, x="group", y="clonotype", size="sample_prevalence",
        hue="expanded_fraction_within_group", sizes=(20, 260), palette="Reds",
        edgecolor="#606060", linewidth=0.4, ax=ax,
    )
    ax.set(xlabel="Clinical group", ylabel="Selected clonotype")
    ax.tick_params(axis="x", rotation=45)
    ax.legend(bbox_to_anchor=(1.02, 1), loc="upper left", frameon=False, fontsize=6)
    finish_figure(fig, out / "Fig5g_pathogen_biased_clonotypes", 8, 5.2)


def plot_network(frame, path, title, seed):
    graph = nx.Graph()
    for row in frame.itertuples(index=False):
        source = getattr(row, "source_cell_subtype")
        target = getattr(row, "target_cell_subtype")
        weight = float(getattr(row, "jaccard_similarity"))
        if weight > 0:
            graph.add_edge(source, target, weight=weight)
    if graph.number_of_nodes() == 0:
        return
    position = nx.spring_layout(graph, seed=seed, weight="weight", k=1.1)
    weights = [graph[u][v]["weight"] for u, v in graph.edges]
    fig, ax = plt.subplots()
    nx.draw_networkx_edges(graph, position, width=[0.5 + 8 * value for value in weights], edge_color=weights, edge_cmap=plt.cm.Reds, alpha=0.65, ax=ax)
    nx.draw_networkx_nodes(graph, position, node_color=PALETTE[1], node_size=180, edgecolors="white", linewidths=0.8, ax=ax)
    nx.draw_networkx_labels(graph, position, font_size=6, ax=ax)
    ax.set_title(title)
    ax.axis("off")
    finish_figure(fig, path, 7, 5.5)


def plot_fig5h(edges, out, seed):
    for index, (anchor, frame) in enumerate(edges.groupby("anchor_clonotype", observed=True), start=1):
        plot_network(frame, out / f"Fig5h_{anchor}_tcr_similarity_network", f"{anchor} TCR similarity network", seed + index)


def main():
    args, repo, out = setup("CD8 TCR recurrence and similarity networks", "11_cd8_tcr_similarity")
    clonotypes = repo.table(23, "ST23a")
    highlighted = repo.table(24, "ST24b")
    edges = repo.table(25, "ST25b")
    edges["jaccard_similarity"] = pd.to_numeric(edges["jaccard_similarity"], errors="coerce").fillna(0)
    write_table(clonotypes, out / "top_cd8_clonotypes.tsv")
    write_table(highlighted, out / "highlighted_cd8_clonotypes.tsv")
    write_table(edges, out / "cd8_tcr_similarity_edges.tsv")
    plot_fig5f(clonotypes, out)
    plot_fig5g(highlighted, out)
    plot_fig5h(edges, out, args.seed)


if __name__ == "__main__":
    main()
