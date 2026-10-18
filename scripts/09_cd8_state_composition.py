from _bootstrap import PALETTE, donor_composition, finish_figure, kruskal_by_feature, plt, setup, write_table
import seaborn as sns


GROUP_ORDER = ["HIV", "HIV_TBM", "HIV_CMV", "HIV_NS", "HIV_TE", "HIV_CM"]


def panel_grid(composition, tests, states, path, columns=5):
    rows = (len(states) + columns - 1) // columns
    fig, axes = plt.subplots(rows, columns, squeeze=False)
    data = composition.copy()
    data["percentage"] = 100 * data["proportion"]
    for ax, state in zip(axes.flat, states):
        frame = data[data["SubType_260501"].eq(state)]
        sns.boxplot(data=frame, x="Disease group", y="percentage", order=GROUP_ORDER, hue="Disease group", palette=PALETTE, legend=False, fliersize=0, ax=ax)
        sns.stripplot(data=frame, x="Disease group", y="percentage", order=GROUP_ORDER, color="#303030", size=2, alpha=0.65, ax=ax)
        test = tests[tests["SubType_260501"].eq(state)]
        if len(test):
            ax.text(0.02, 0.98, f"KW P={test.iloc[0]['p_value']:.3g}\nBH q={test.iloc[0]['q_value']:.3g}", transform=ax.transAxes, va="top", fontsize=5.5)
        ax.set_title(state, fontsize=7)
        ax.tick_params(axis="x", rotation=60, labelsize=5)
        ax.set(xlabel="", ylabel="Cells (% of donor CSF cells)")
    for ax in axes.flat[len(states):]:
        ax.axis("off")
    finish_figure(fig, path, 11, max(4.2, rows * 3.1))


def main():
    _, repo, out = setup("CD8 T-cell state composition", "09_cd8_composition")
    cells = repo.source(1, usecols=["orig.ident", "Disease group", "Type", "SubType_260501"])
    cells = cells[cells["Type"].eq("CD8T")]
    composition = donor_composition(cells, "orig.ident", "Disease group", "SubType_260501")
    tests = kruskal_by_feature(composition, "SubType_260501", "Disease group", "proportion")
    write_table(composition, out / "cd8_donor_composition.tsv")
    write_table(tests, out / "cd8_global_tests.tsv")
    states = composition.groupby("SubType_260501")["cells"].sum().nlargest(13).index
    selected = tests.sort_values("p_value").head(4)["SubType_260501"]
    panel_grid(composition, tests, list(selected), out / "Fig5b_selected_cd8_state_composition", columns=4)
    panel_grid(composition, tests, list(states), out / "EDFig6a_cd8_state_composition", columns=5)


if __name__ == "__main__":
    main()
