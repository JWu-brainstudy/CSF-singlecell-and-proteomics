from _bootstrap import finish_figure, pd, plt, setup, spearman_table, write_table
import seaborn as sns
from scipy.stats import spearmanr


def scatter_panel(data, x, y, path, title):
    frame = data[[x, y]].apply(pd.to_numeric, errors="coerce").dropna()
    rho, p_value = spearmanr(frame[x], frame[y])
    fig, ax = plt.subplots()
    sns.regplot(data=frame, x=x, y=y, scatter_kws={"s": 24, "alpha": 0.75}, line_kws={"linewidth": 1.4}, ci=95, ax=ax)
    ax.set_title(title)
    ax.text(0.03, 0.97, f"Spearman ρ={rho:.2f}\nP={p_value:.3g}; n={len(frame)}", transform=ax.transAxes, va="top")
    finish_figure(fig, path, 5.5, 4.5)


def main():
    _, repo, out = setup("CD4 cross-compartment associations", "08_cd4_cross_compartment")
    data = repo.table(18, "ST18c")
    blood = [column for column in data.columns if str(column).startswith("Blood")]
    csf = [column for column in data.columns if str(column).startswith("CSF")]
    pairs = [(left, right) for left in blood for right in csf]
    correlations = spearman_table(data, pairs)
    write_table(correlations, out / "cross_compartment_spearman.tsv")
    scatter_panel(data, blood[0], csf[0], out / "Fig4h_csf_cd4_vs_blood_cd4", "CSF CD4 abundance vs blood CD4")
    scatter_panel(data, blood[1], csf[1], out / "Fig4i_csf_vs_blood_cd4_cd8_ratio", "CSF vs blood CD4/CD8 ratio")
    matrix = correlations.pivot(index="measure_2", columns="measure_1", values="rho").reindex(index=csf, columns=blood)
    qvalues = correlations.pivot(index="measure_2", columns="measure_1", values="q_value").reindex(index=csf, columns=blood)
    annotations = matrix.copy().astype(object)
    for row in matrix.index:
        for column in matrix.columns:
            value = matrix.loc[row, column]
            annotations.loc[row, column] = "" if pd.isna(value) else f"{value:.2f}{'*' if qvalues.loc[row, column] < 0.05 else ''}"
    fig, ax = plt.subplots()
    sns.heatmap(matrix, cmap="vlag", center=0, vmin=-0.6, vmax=0.6, annot=annotations, fmt="", linewidths=1, cbar_kws={"label": "Spearman ρ"}, ax=ax)
    ax.set(xlabel="Peripheral blood measure", ylabel="CSF CD4 measure")
    finish_figure(fig, out / "Fig4j_cd4_subtype_concordance", 7.5, 5.5)


if __name__ == "__main__":
    main()
