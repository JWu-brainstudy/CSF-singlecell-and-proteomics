from _bootstrap import finish_figure, heatmap, pd, plt, setup, write_table
import networkx as nx
import seaborn as sns


def plot_differential_network(data, path, title):
    data = data.copy()
    data["Difference"] = pd.to_numeric(data["Difference"], errors="coerce").fillna(0)
    graph = nx.DiGraph()
    for row in data.itertuples(index=False):
        graph.add_edge(row.Sender, row.Receiver, weight=abs(row.Difference), signed=row.Difference)
    position = nx.circular_layout(graph)
    edges = list(graph.edges)
    widths = [0.5 + 5 * graph[u][v]["weight"] / max(1e-9, max(abs(data["Difference"]))) for u, v in edges]
    colors = ["#D1495B" if graph[u][v]["signed"] > 0 else "#5B9BD5" for u, v in edges]
    fig, ax = plt.subplots()
    nx.draw_networkx_nodes(graph, position, node_size=430, node_color="#E5E5E5", edgecolors="#505050", ax=ax)
    nx.draw_networkx_labels(graph, position, font_size=7, ax=ax)
    nx.draw_networkx_edges(graph, position, width=widths, edge_color=colors, alpha=0.75, arrows=True, arrowsize=10, connectionstyle="arc3,rad=0.08", ax=ax)
    ax.set_title(title)
    ax.axis("off")
    finish_figure(fig, path, 6.5, 5.5)


def main():
    _, repo, out = setup("Microglia-like differential communication", "14_microglia_cellchat")
    interactions = repo.table(28)
    interactions["prob"] = pd.to_numeric(interactions["prob"], errors="coerce")
    pathway = interactions.groupby(["relation", "pathway_name"], observed=True)["prob"].sum().rename("interaction_strength").reset_index()
    mhc2 = repo.table(30, "ST30a")
    mhc1 = repo.table(30, "ST30b")
    write_table(pathway, out / "microglia_outgoing_pathways.tsv")
    write_table(mhc2, out / "mhc_ii_differential_edges.tsv")
    write_table(mhc1, out / "mhc_i_differential_edges.tsv")
    receiver = interactions.groupby(["relation", "target"], observed=True)["prob"].sum().rename("incoming_strength").reset_index()
    selected = receiver.groupby("relation", observed=True).apply(lambda frame: frame.nlargest(6, "incoming_strength"), include_groups=False).reset_index()
    fig, axes = plt.subplots(1, max(1, selected["relation"].nunique()), squeeze=False)
    for ax, (relation, frame) in zip(axes.flat, selected.groupby("relation", observed=True)):
        ax.pie(frame["incoming_strength"], labels=frame["target"].astype(str), startangle=90, colors=sns.color_palette("Set3", len(frame)), wedgeprops={"width": 0.45, "edgecolor": "white"}, textprops={"fontsize": 5.5})
        ax.set_title(str(relation), fontsize=8)
    finish_figure(fig, out / "Fig6f_microglia_receiver_landscape", max(7, selected["relation"].nunique() * 3.2), 4)
    top = pathway.groupby("pathway_name")["interaction_strength"].sum().nlargest(15).index
    plot_data = pathway[pathway["pathway_name"].isin(top)].pivot_table(index="pathway_name", columns="relation", values="interaction_strength", fill_value=0)
    fig, ax = plt.subplots()
    image = heatmap(ax, plot_data.to_numpy(), plot_data.columns, plot_data.index, cmap="viridis")
    fig.colorbar(image, ax=ax, label="Summed communication probability")
    ax.set(xlabel="Comparison", ylabel="Signaling pathway")
    finish_figure(fig, out / "Fig6g_microglia_outgoing_communication", 8, 6)
    for label, data, prefix in [("MHC-II", mhc2, "Fig6h"), ("MHC-I", mhc1, "Fig6i")]:
        for index, (comparison, frame) in enumerate(data.groupby("Comparison group", observed=True), start=1):
            plot_differential_network(frame, out / f"{prefix}_{index}_differential_network", f"{label}: {comparison}")


if __name__ == "__main__":
    main()
