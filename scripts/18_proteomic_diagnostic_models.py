from _bootstrap import finish_figure, np, pd, plt, setup, write_table
from sklearn.calibration import calibration_curve
from sklearn.metrics import auc, brier_score_loss, roc_curve
import seaborn as sns


def plot_prediction_distributions(data, prefix, out):
    data = data.copy()
    data["Mean out-of-fold probability"] = pd.to_numeric(data["Mean out-of-fold probability"], errors="coerce")
    comparisons = data["Comparison"].dropna().unique()
    fig, axes = plt.subplots(1, len(comparisons), squeeze=False)
    for ax, comparison in zip(axes.flat, comparisons):
        frame = data[data["Comparison"].eq(comparison)]
        order = list(frame.groupby("Observed group")["Mean out-of-fold probability"].median().sort_values(ascending=False).index)
        values = [frame.loc[frame["Observed group"].eq(group), "Mean out-of-fold probability"].dropna().to_numpy() for group in order]
        ax.boxplot(values, showfliers=False, patch_artist=True, boxprops={"facecolor": "#8EC9DC"})
        ax.set_xticks(range(1, len(order) + 1), order)
        random = np.random.default_rng(2026)
        for position, sample in enumerate(values, start=1):
            ax.scatter(position + random.normal(0, 0.045, len(sample)), sample, color="#303030", s=8, alpha=0.6)
        threshold = pd.to_numeric(frame["Internal Youden threshold"], errors="coerce").dropna()
        if len(threshold):
            ax.axhline(threshold.iloc[0], color="#606060", linestyle="--", linewidth=0.8)
        ax.set_title(str(comparison), fontsize=8)
        ax.tick_params(axis="x", rotation=50, labelsize=5)
        ax.set(xlabel="", ylabel="Mean out-of-fold probability")
    finish_figure(fig, out / f"EDFig8_{prefix}_prediction_distributions", max(8, 3.2 * len(comparisons)), 4.8)


def plot_coefficients(repo, number, prefix, out):
    terms = repo.table(number, f"ST{number}d").copy()
    terms["Coefficient"] = pd.to_numeric(terms["Coefficient"], errors="coerce")
    terms = terms[terms["Model term"].astype(str).ne("Intercept")].dropna(subset=["Coefficient"])
    fig, axes = plt.subplots(1, max(1, terms["Comparison"].nunique()), squeeze=False)
    for ax, (comparison, frame) in zip(axes.flat, terms.groupby("Comparison", observed=True)):
        frame = frame.sort_values("Coefficient")
        ax.hlines(range(len(frame)), 0, frame["Coefficient"], color="#A0A0A0")
        ax.scatter(frame["Coefficient"], range(len(frame)), color="#2F6690")
        ax.set_yticks(range(len(frame)), frame["Model term"].astype(str), fontsize=6)
        ax.axvline(0, color="#606060", linewidth=0.7)
        ax.set(title=str(comparison), xlabel="Standardized coefficient")
    finish_figure(fig, out / f"EDFig8_{prefix}_model_coefficients", max(8, 3.2 * terms["Comparison"].nunique()), 4.8)


def model_roc(repo, number):
    data = repo.table(number, f"ST{number}c")
    rows = []
    curves = []
    for comparison, frame in data.groupby("Comparison", observed=True):
        frame = frame.copy()
        frame["Observed outcome (1=case)"] = pd.to_numeric(frame["Observed outcome (1=case)"], errors="coerce")
        frame["Mean out-of-fold probability"] = pd.to_numeric(frame["Mean out-of-fold probability"], errors="coerce")
        frame = frame.dropna(subset=["Observed outcome (1=case)", "Mean out-of-fold probability"])
        y = frame["Observed outcome (1=case)"].astype(int)
        score = frame["Mean out-of-fold probability"].astype(float)
        if len(frame) == 0 or y.nunique() < 2:
            continue
        fpr, tpr, _ = roc_curve(y, score)
        value = auc(fpr, tpr)
        rows.append({"comparison": comparison, "n": len(frame), "auc": value, "brier_score": brier_score_loss(y, score)})
        curves.append((comparison, fpr, tpr, value, y.to_numpy(), score.to_numpy()))
    return data, rows, curves


def main():
    _, repo, out = setup("Proteomic diagnostic models", "18_diagnostic_models")
    tbm, tbm_rows, tbm_curves = model_roc(repo, 37)
    cm, cm_rows, cm_curves = model_roc(repo, 38)
    metrics = pd.DataFrame(tbm_rows + cm_rows)
    write_table(tbm, out / "tbm_out_of_fold_predictions.tsv")
    write_table(cm, out / "cm_out_of_fold_predictions.tsv")
    write_table(metrics, out / "diagnostic_auc.tsv")
    plot_prediction_distributions(tbm, "tbm", out)
    plot_prediction_distributions(cm, "cm", out)
    plot_coefficients(repo, 37, "tbm", out)
    plot_coefficients(repo, 38, "cm", out)
    fig, axes = plt.subplots(1, 2)
    for ax, title, curves in zip(axes, ["TBM", "CM"], [tbm_curves, cm_curves]):
        for comparison, fpr, tpr, value, _, _ in curves:
            ax.plot(fpr, tpr, label=f"{comparison}; AUC={value:.2f}")
        ax.plot([0, 1], [0, 1], color="#A0A0A0", linestyle="--")
        ax.set(title=title, xlabel="1 − specificity", ylabel="Sensitivity")
        ax.legend(frameon=False, fontsize=6)
    finish_figure(fig, out / "Fig7h_diagnostic_roc_curves", 9, 4.5)
    fig, axes = plt.subplots(1, 2)
    for ax, title, curves in zip(axes, ["TBM", "CM"], [tbm_curves, cm_curves]):
        for comparison, _, _, _, y, score in curves:
            observed, predicted = calibration_curve(y, score, n_bins=6, strategy="quantile")
            ax.plot(predicted, observed, marker="o", label=comparison)
        ax.plot([0, 1], [0, 1], color="#A0A0A0", linestyle="--")
        ax.set(title=title, xlabel="Mean predicted probability", ylabel="Observed event fraction")
        ax.legend(frameon=False, fontsize=6)
    finish_figure(fig, out / "EDFig8_diagnostic_calibration", 9, 4.5)


if __name__ == "__main__":
    main()
