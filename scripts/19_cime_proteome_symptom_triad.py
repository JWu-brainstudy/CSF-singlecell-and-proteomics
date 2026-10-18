from _bootstrap import finish_figure, heatmap, pd, plt, setup, write_table
import networkx as nx
import numpy as np


def plot_fig7i(screen, out):
    data = screen.copy()
    data["q_BH_720_by_weakest_edge"] = pd.to_numeric(data["q_BH_720_by_weakest_edge"], errors="coerce")
    data["minimum_absolute_rho"] = pd.to_numeric(data["minimum_absolute_rho"], errors="coerce")
    data = data.sort_values(["q_BH_720_by_weakest_edge", "minimum_absolute_rho"], ascending=[True, False]).head(30)
    graph = nx.Graph()
    for row in data.itertuples(index=False):
        symptom = f"Symptom:{row.domain_code}"
        cime = f"CIME:{row.CIME_program}"
        protein = f"Protein:{row.proteomic_module}"
        graph.add_edge(symptom, cime, weight=abs(float(row.rho_CIME_symptom)), signed=float(row.rho_CIME_symptom))
        graph.add_edge(symptom, protein, weight=abs(float(row.rho_protein_symptom)), signed=float(row.rho_protein_symptom))
        graph.add_edge(cime, protein, weight=abs(float(row.rho_CIME_protein)), signed=float(row.rho_CIME_protein))
    nodes = sorted(graph.nodes)
    angles = np.linspace(0, 2 * np.pi, len(nodes), endpoint=False)
    position = {node: np.array([np.cos(angle), np.sin(angle)]) for node, angle in zip(nodes, angles)}
    colors = ["#7B4EA3" if node.startswith("Symptom") else ("#2C9BC7" if node.startswith("CIME") else "#E77C5C") for node in nodes]
    widths = [0.4 + 4 * graph[u][v]["weight"] for u, v in graph.edges]
    edge_colors = ["#D1495B" if graph[u][v]["signed"] > 0 else "#4C78A8" for u, v in graph.edges]
    fig, ax = plt.subplots()
    nx.draw_networkx_edges(graph, position, width=widths, edge_color=edge_colors, alpha=0.35, connectionstyle="arc3,rad=0.15", ax=ax)
    nx.draw_networkx_nodes(graph, position, node_color=colors, node_size=240, edgecolors="white", ax=ax)
    nx.draw_networkx_labels(graph, position, labels={node: node.split(":", 1)[1] for node in nodes}, font_size=6, ax=ax)
    ax.axis("off")
    finish_figure(fig, out / "Fig7i_cime_proteome_symptom_network", 8, 8)


def plot_edfig9_forest(pooled, out):
    data = pooled.copy()
    for column in ["odds_ratio", "ci_low", "ci_high", "p"]:
        data[column] = pd.to_numeric(data[column], errors="coerce")
    data = data.dropna(subset=["odds_ratio", "ci_low", "ci_high"])
    labels = data["predictor"].astype(str) + " | " + data["model"].astype(str)
    y = np.arange(len(data))
    fig, ax = plt.subplots()
    ax.errorbar(data["odds_ratio"], y, xerr=[data["odds_ratio"] - data["ci_low"], data["ci_high"] - data["odds_ratio"]], fmt="o", color="#2F6690", ecolor="#7A9EBC", capsize=2)
    ax.axvline(1, color="#606060", linestyle="--", linewidth=0.8)
    ax.set_xscale("log")
    ax.set_yticks(y, labels, fontsize=6)
    ax.set(xlabel="Odds ratio per rank SD (95% CI)", ylabel="")
    finish_figure(fig, out / "EDFig9e_symptom_item_models", 9, max(4, 0.3 * len(data)))


def main():
    _, repo, out = setup("CIME-proteome-symptom triad", "19_triad")
    screen = repo.table(39, "ST39b")
    primary = repo.table(39, "ST39c")
    pooled = repo.table(39, "ST39d")
    write_table(screen, out / "triad_exploratory_screen.tsv")
    write_table(primary, out / "triad_primary_estimates.tsv")
    write_table(pooled, out / "triad_pooled_item_models.tsv")
    for column in ["rho_CIME_symptom", "rho_protein_symptom", "rho_CIME_protein"]:
        screen[column] = pd.to_numeric(screen[column], errors="coerce")
    score = screen[["rho_CIME_symptom", "rho_protein_symptom", "rho_CIME_protein"]].abs().min(axis=1)
    top = screen.assign(triad_strength=score).nlargest(25, "triad_strength")
    labels = top["domain_code"].astype(str) + " | " + top["CIME_program"].astype(str) + " | " + top["proteomic_module"].astype(str)
    matrix = top[["rho_CIME_symptom", "rho_protein_symptom", "rho_CIME_protein"]].to_numpy()
    fig, ax = plt.subplots()
    image = heatmap(ax, matrix, ["CIME–symptom", "Protein–symptom", "CIME–protein"], labels, center=0)
    fig.colorbar(image, ax=ax, label="Spearman rho")
    ax.set(ylabel="Triad")
    finish_figure(fig, out / "cime_proteome_symptom_triad", 7.5, 8)
    plot_fig7i(screen, out)
    plot_edfig9_forest(pooled, out)


if __name__ == "__main__":
    main()
