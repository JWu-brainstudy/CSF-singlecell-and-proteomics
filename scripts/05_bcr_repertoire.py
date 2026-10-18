from _bootstrap import PALETTE, finish_figure, kruskal_by_feature, mann_whitney_by_feature, np, pd, plt, setup, write_table
import seaborn as sns


GROUP_ORDER = ["HIV", "HIV_TBM", "HIV_CMV", "HIV_NS", "HIV_TE", "HIV_CM"]


def stacked_panel(data, index, column, value, path, ylabel, colormap):
    frame = data.copy()
    frame[value] = pd.to_numeric(frame[value], errors="coerce")
    pivot = frame.pivot_table(index=index, columns=column, values=value, aggfunc="sum", fill_value=0).reindex(GROUP_ORDER).fillna(0)
    fig, ax = plt.subplots()
    pivot.plot.bar(stacked=True, ax=ax, colormap=colormap, width=0.82)
    ax.set(xlabel="Clinical group", ylabel=ylabel)
    ax.tick_params(axis="x", rotation=40)
    ax.legend(title=column, frameon=False, fontsize=6, bbox_to_anchor=(1.02, 1), loc="upper left")
    finish_figure(fig, path, 8, 5)


def plot_fig3g(sankey, out):
    data = sankey.copy()
    data["Cell count"] = pd.to_numeric(data["Cell count"], errors="coerce").fillna(0)
    dimensions = ["B/plasma subtype", "IGHV", "Light-chain V", "Clinical group"]
    top = {column: data.groupby(column)["Cell count"].sum().nlargest(10).index for column in dimensions[1:3]}
    data = data[data["IGHV"].isin(top["IGHV"]) & data["Light-chain V"].isin(top["Light-chain V"])]
    positions = {column: {label: i for i, label in enumerate(data.groupby(column)["Cell count"].sum().sort_values().index)} for column in dimensions}
    fig, ax = plt.subplots()
    maximum = max(data["Cell count"].max(), 1)
    for _, row in data.iterrows():
        labels = [row[column] for column in dimensions]
        y = [positions[column][label] for column, label in zip(dimensions, labels)]
        ax.plot(range(4), y, color=PALETTE[GROUP_ORDER.index(labels[-1])] if labels[-1] in GROUP_ORDER else "#999999", alpha=0.15 + 0.7 * row["Cell count"] / maximum, linewidth=0.4 + 3 * row["Cell count"] / maximum)
    ax.set_xticks(range(4), dimensions)
    ax.set_ylabel("Category rank")
    ax.set_title("B-cell subtype–V gene–clinical group flow")
    finish_figure(fig, out / "Fig3g_bcr_alluvial", 9, 5.5)


def plot_fig3h(pairs, out):
    data = pairs.copy()
    for column in ["Within-group frequency", "Donors represented", "Pair cells in group"]:
        data[column] = pd.to_numeric(data[column], errors="coerce")
    top = data.groupby("IGHV-light-chain V pair")["Pair cells in group"].sum().nlargest(20).index
    data = data[data["IGHV-light-chain V pair"].isin(top)]
    fig, ax = plt.subplots()
    sns.scatterplot(data=data, x="Clinical group", y="IGHV-light-chain V pair", size="Donors represented", hue="Within-group frequency", sizes=(20, 240), palette="viridis", edgecolor="#505050", linewidth=0.3, ax=ax)
    ax.tick_params(axis="x", rotation=45)
    ax.legend(frameon=False, fontsize=6, bbox_to_anchor=(1.02, 1), loc="upper left")
    ax.set(xlabel="Clinical group", ylabel="IGHV–light-chain V pair")
    finish_figure(fig, out / "Fig3h_bcr_v_pair_recurrence", 9, 7)


def main():
    _, repo, out = setup("B-cell receptor repertoire analysis", "05_bcr_repertoire")
    clone_summary = repo.table(12, "ST12b")
    expanded_subtype = repo.table(12, "ST12c")
    isotype = repo.table(12, "ST12d")
    sankey = repo.table(13, "ST13a")
    pairs = repo.table(13, "ST13b")
    metadata = repo.source(1, usecols=["CellName", "orig.ident", "Disease group"])
    receptor = repo.source(2, usecols=["cellnames", "cc_aa_alignment_size"])
    donor = receptor.merge(metadata, left_on="cellnames", right_on="CellName", how="left")
    donor["expanded"] = pd.to_numeric(donor["cc_aa_alignment_size"], errors="coerce").gt(2)
    donor = donor.groupby(["orig.ident", "Disease group"], observed=True).agg(bcr_cells=("cellnames", "size"), expanded_cells=("expanded", "sum")).reset_index()
    donor["expanded_fraction"] = donor["expanded_cells"] / donor["bcr_cells"]
    donor["analysis"] = "BCR expansion"
    global_test = kruskal_by_feature(donor, "analysis", "Disease group", "expanded_fraction")
    contrasts = mann_whitney_by_feature(donor, "analysis", "Disease group", "expanded_fraction", "HIV")
    write_table(clone_summary, out / "bcr_clone_size_distribution.tsv")
    write_table(isotype, out / "bcr_heavy_chain_isotypes.tsv")
    write_table(pairs, out / "bcr_v_gene_pairs.tsv")
    write_table(donor, out / "bcr_donor_expansion_recomputed.tsv")
    write_table(global_test, out / "ST12a_bcr_expansion_kruskal_recomputed.tsv")
    write_table(contrasts, out / "ST12a_bcr_expansion_mann_whitney_recomputed.tsv")
    write_table(expanded_subtype, out / "bcr_expanded_subtypes.tsv")
    write_table(sankey, out / "bcr_alluvial_input.tsv")
    stacked_panel(clone_summary, "Clinical group", "Clone-size category", "Proportion", out / "Fig3d_bcr_clone_size_distribution", "Fraction of BCR cells", "Blues")
    stacked_panel(expanded_subtype, "Clinical group", "B/plasma subtype", "Proportion of group expanded pool", out / "Fig3e_expanded_bcr_subtypes", "Fraction of expanded BCR cells", "tab20")
    stacked_panel(isotype, "Clinical group", "Heavy-chain constant call", "Proportion", out / "Fig3f_bcr_isotype_usage", "Fraction of heavy-chain annotated cells", "Set3")
    plot_fig3g(sankey, out)
    plot_fig3h(pairs, out)


if __name__ == "__main__":
    main()
