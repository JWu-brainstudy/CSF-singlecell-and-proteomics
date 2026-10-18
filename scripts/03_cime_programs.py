from _bootstrap import PALETTE, finish_figure, np, pd, plt, setup, write_table
import seaborn as sns
from matplotlib.colors import ListedColormap
from scipy.cluster.hierarchy import dendrogram, linkage
from sklearn.preprocessing import StandardScaler


GROUP_ORDER = ["HIV", "HIV-TBM", "HIV-CMV", "HIV-NS", "HIV-TE", "HIV-CM"]
PROGRAM_ORDER = ["P1", "P2", "P3", "P4", "P5", "P6", "P7", "P9"]
PROGRAM_LABELS = {
    "P1": "CIME-CD8Eff", "P2": "CIME-gdT", "P3": "CIME-CD4", "P4": "CIME-NK",
    "P5": "CIME-Mye", "P6": "CIME-CD8Act", "P7": "CIME-B/Plasma", "P9": "CIME-Neu",
}


def normalize_group(values):
    return values.astype(str).str.replace("_", "-", regex=False)


def enrichment_table(repo):
    data = repo.table(9, "ST09a").copy()
    numeric = [
        "a: group and program", "b: group other programs", "c: other groups and program",
        "d: other groups other programs", "Group n", "Two-sided Fisher P", "BH q (48 tests)",
    ]
    for column in numeric:
        data[column] = pd.to_numeric(data[column], errors="coerce")
    data["program_total"] = data["a: group and program"] + data["c: other groups and program"]
    data["cohort_total"] = data[[
        "a: group and program", "b: group other programs", "c: other groups and program",
        "d: other groups other programs",
    ]].sum(axis=1)
    data["expected_donors"] = data["Group n"] * data["program_total"] / data["cohort_total"]
    data["observed_expected_ratio"] = data["a: group and program"] / data["expected_donors"]
    return data


def plot_fig2a(repo, out):
    loadings = repo.table(8, "ST08c")
    programs = [column for column in loadings.columns if column.startswith("P")]
    matrix = loadings[programs].apply(pd.to_numeric, errors="coerce").fillna(0)
    tree = linkage(matrix, method="ward", optimal_ordering=True)
    cells = repo.source(1, usecols=["Disease group", "SubType_260501"])
    cells["group"] = normalize_group(cells["Disease group"])
    cells = cells[cells["group"].isin(GROUP_ORDER)]
    counts = pd.crosstab(cells["SubType_260501"], cells["group"])
    counts = counts.reindex(index=loadings.iloc[:, 0].astype(str), columns=GROUP_ORDER, fill_value=0)
    fig = plt.figure()
    grid = fig.add_gridspec(3, 1, height_ratios=[5.0, 1.0, 1.2], hspace=0.08)
    ax_tree = fig.add_subplot(grid[0])
    result = dendrogram(
        tree, labels=loadings.iloc[:, 0].astype(str).to_numpy(), leaf_rotation=90,
        leaf_font_size=5, color_threshold=None, above_threshold_color="#202020", ax=ax_tree,
    )
    order = result["ivl"]
    x = np.arange(len(order)) * 10 + 5
    ordered = counts.reindex(order)
    ax_tree.set(ylabel="Ward distance")
    ax_count = fig.add_subplot(grid[1], sharex=ax_tree)
    ax_count.bar(x, ordered.sum(axis=1), width=8, color="#E87824")
    ax_count.set_ylabel("Cells")
    ax_count.tick_params(axis="x", labelbottom=False)
    ax_prop = fig.add_subplot(grid[2], sharex=ax_tree)
    fractions = ordered.div(ordered.sum(axis=1).replace(0, np.nan), axis=0).fillna(0)
    bottom = np.zeros(len(order))
    for group, color in zip(GROUP_ORDER, PALETTE):
        ax_prop.bar(x, fractions[group], width=8, bottom=bottom, color=color, label=group)
        bottom += fractions[group].to_numpy()
    ax_prop.set(ylabel="Group fraction", ylim=(0, 1))
    ax_prop.tick_params(axis="x", labelbottom=False)
    ax_prop.legend(ncol=6, frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.25))
    finish_figure(fig, out / "Fig2a_hierarchical_cell_state_organization", 12, 8)


def plot_fig2b(repo, out):
    counts = repo.table(7, "ST07a")
    counts["Cell count"] = pd.to_numeric(counts["Cell count"], errors="coerce").fillna(0)
    counts["Clinical group"] = normalize_group(counts["Clinical group"])
    counts = counts[counts["Clinical group"].isin(GROUP_ORDER)]
    group_map = counts.drop_duplicates("Donor ID").set_index("Donor ID")["Clinical group"]
    usage = repo.table(9, "ST09c")
    usage = usage[usage["Included in core categorical analysis"].astype(str).eq("Yes")]
    dominant = usage.set_index("Donor ID")["Dominant CIME"].astype(str)
    donors = [donor for donor in usage["Donor ID"] if donor in group_map.index]
    donors = sorted(donors, key=lambda donor: (
        PROGRAM_ORDER.index(dominant[donor]) if dominant[donor] in PROGRAM_ORDER else 99,
        GROUP_ORDER.index(group_map[donor]), donor,
    ))
    membership = repo.table(8, "ST08b")
    membership["Within-program rank"] = pd.to_numeric(membership["Within-program rank"], errors="coerce")
    membership = membership.sort_values(["Assigned program", "Within-program rank"])
    subtypes = membership["Fine cell subtype"].astype(str).tolist()
    abundance = counts.pivot_table(index="Fine cell subtype", columns="Donor ID", values="Cell count", aggfunc="sum", fill_value=0)
    abundance = abundance.reindex(index=subtypes, columns=donors, fill_value=0)
    proportions = abundance.div(abundance.sum(axis=0).replace(0, np.nan), axis=1).fillna(0)
    score = StandardScaler().fit_transform(proportions.T).T
    score = pd.DataFrame(np.clip(score, -2.5, 2.5), index=proportions.index, columns=proportions.columns)
    program_codes = [PROGRAM_ORDER.index(dominant[donor]) for donor in donors]
    group_codes = [GROUP_ORDER.index(group_map[donor]) for donor in donors]
    fig = plt.figure()
    grid = fig.add_gridspec(3, 1, height_ratios=[0.18, 0.18, 6], hspace=0.02)
    ax_program = fig.add_subplot(grid[0])
    ax_program.imshow([program_codes], aspect="auto", cmap=ListedColormap(sns.color_palette("Set3", len(PROGRAM_ORDER))))
    ax_program.set(xticks=[], yticks=[0], yticklabels=["CIME"])
    ax_group = fig.add_subplot(grid[1])
    ax_group.imshow([group_codes], aspect="auto", cmap=ListedColormap(PALETTE))
    ax_group.set(xticks=[], yticks=[0], yticklabels=["Group"])
    ax_heat = fig.add_subplot(grid[2])
    sns.heatmap(
        score, cmap="vlag", center=0, xticklabels=False, yticklabels=True,
        cbar_kws={"label": "Scaled abundance"}, ax=ax_heat,
    )
    ax_heat.set(xlabel="Donors ordered by dominant CIME", ylabel="Fine cell subtype")
    ax_heat.tick_params(axis="y", labelsize=5)
    finish_figure(fig, out / "Fig2b_cime_program_heatmap", 11, 12)


def plot_fig2c(enrichment, out):
    row_order = ["P7", "P6", "P4", "P3", "P2", "P9", "P5", "P1"]
    ratio = enrichment.pivot(index="Program", columns="Disease group", values="observed_expected_ratio")
    ratio = ratio.reindex(index=row_order, columns=GROUP_ORDER)
    annotations = ratio.map(lambda value: "+++" if value > 2.5 else ("++" if value > 2 else ("+" if value > 1.5 else "")))
    fig, ax = plt.subplots()
    sns.heatmap(
        ratio, cmap="Purples", vmin=0, vmax=max(3, np.nanmax(ratio.to_numpy())),
        annot=annotations, fmt="", linewidths=1, linecolor="white",
        cbar_kws={"label": "Observed / expected"}, ax=ax,
    )
    ax.set(xlabel="Clinical group", ylabel="CIME program")
    ax.set_yticklabels([PROGRAM_LABELS[program] for program in row_order], rotation=0)
    ax.tick_params(axis="x", rotation=45)
    finish_figure(fig, out / "Fig2c_cime_enrichment_matrix", 7.5, 5.8)


def plot_lollipop(data, program, path, title):
    frame = data[data["Program"].eq(program)].set_index("Disease group").reindex(GROUP_ORDER)
    values = frame["observed_expected_ratio"].fillna(0).to_numpy()
    x = np.arange(len(GROUP_ORDER))
    fig, ax = plt.subplots()
    ax.vlines(x, 0, values, color="#B0B0B0", linewidth=2)
    ax.scatter(x, values, color=PALETTE, s=90, zorder=3)
    for position, value in zip(x, values):
        ax.text(position, value + max(values.max() * 0.04, 0.08), f"{value:.2f}", ha="center", va="bottom", fontsize=7)
    ax.axhline(1, color="#606060", linestyle="--", linewidth=0.8)
    ax.set(xticks=x, xticklabels=GROUP_ORDER, ylabel="Observed / expected", title=title)
    ax.tick_params(axis="x", rotation=45)
    finish_figure(fig, path, 5.8, 4.2)


def plot_fig2d_e_and_ed2(enrichment, repo, out):
    plot_lollipop(enrichment, "P4", out / "Fig2d_cime_nk_enrichment", "CIME-NK")
    plot_lollipop(enrichment, "P7", out / "Fig2e_cime_b_plasma_enrichment", "CIME-B/Plasma")
    coverage = repo.table(8, "ST08a")
    coverage["Dominant donors, all 73"] = pd.to_numeric(coverage["Dominant donors, all 73"], errors="coerce")
    fig, ax = plt.subplots()
    colors = ["#B0B0B0" if program in ["P8", "P10"] else "#4C78A8" for program in coverage["Program"]]
    ax.barh(coverage["Program"], coverage["Dominant donors, all 73"], color=colors)
    ax.invert_yaxis()
    ax.set(xlabel="Dominant donors", ylabel="Program", title="Dominant-program donor coverage")
    finish_figure(fig, out / "EDFig2c_dominant_program_coverage", 5.2, 4.6)
    for panel, program in zip("defghi", ["P1", "P2", "P3", "P5", "P6", "P9"]):
        plot_lollipop(enrichment, program, out / f"EDFig2{panel}_{program.lower()}_enrichment", f"{program} enrichment")


def main():
    _, repo, out = setup("CIME program structure and disease enrichment", "03_cime")
    enrichment = enrichment_table(repo)
    write_table(repo.table(8, "ST08c"), out / "cime_loading_matrix.tsv")
    write_table(repo.table(9, "ST09c"), out / "donor_cime_usage.tsv")
    write_table(enrichment, out / "cime_enrichment_statistics.tsv")
    plot_fig2a(repo, out)
    plot_fig2b(repo, out)
    plot_fig2c(enrichment, out)
    plot_fig2d_e_and_ed2(enrichment, repo, out)


if __name__ == "__main__":
    main()
