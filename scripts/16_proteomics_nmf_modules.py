from _bootstrap import finish_figure, heatmap, pd, plt, setup, write_table
from sklearn.decomposition import NMF
from sklearn.impute import SimpleImputer
from sklearn.preprocessing import MinMaxScaler
import seaborn as sns


def refit_source_nmf(repo, seed):
    proteins = repo.source(4)
    annotation = ["ProteinAccessions", "Accession", "PG.Genes", "PG.ProteinDescriptions", "PG.ProteinNames"]
    samples = [column for column in proteins.columns if column not in annotation and not str(column).startswith("QC")]
    matrix = proteins[samples].apply(pd.to_numeric, errors="coerce").T
    accessions = proteins["Accession"].astype(str)
    matrix.columns = accessions
    matrix = matrix.loc[:, matrix.notna().any(axis=0)]
    matrix = pd.DataFrame(SimpleImputer(strategy="median").fit_transform(matrix), index=samples, columns=matrix.columns)
    matrix = pd.DataFrame(MinMaxScaler().fit_transform(matrix), index=matrix.index, columns=matrix.columns)
    model = NMF(n_components=8, init="nndsvda", solver="cd", max_iter=2000, random_state=seed)
    scores = pd.DataFrame(model.fit_transform(matrix), index=samples, columns=[f"M{i}" for i in range(1, 9)]).reset_index(names="sample")
    loadings = pd.DataFrame(model.components_.T, index=matrix.columns, columns=[f"M{i}" for i in range(1, 9)]).reset_index(names="Accession")
    return scores, loadings, model.reconstruction_err_


def plot_fig7d(activity, definitions, modules, out):
    ordered = activity.sort_values("group")
    fig, axes = plt.subplots(1, 2, gridspec_kw={"width_ratios": [5, 1.2]})
    sns.heatmap(ordered[modules], cmap="magma", yticklabels=False, cbar_kws={"label": "NMF activity"}, ax=axes[0])
    axes[0].set(xlabel="Proteomic module", ylabel="Samples ordered by clinical group")
    size = definitions.set_index("Module").reindex(modules)
    size["Module_size"] = pd.to_numeric(size["Module_size"], errors="coerce")
    axes[1].barh(modules, size["Module_size"], color="#4C78A8")
    axes[1].invert_yaxis()
    axes[1].set(xlabel="Proteins", ylabel="")
    finish_figure(fig, out / "Fig7d_proteomic_nmf_activity_and_size", 10, 7)


def plot_fig7g(associations, out):
    for metric, suffix in [("HIV_related_Cohens_d", "hiv_related"), ("Vs_HIV_Cohens_d", "versus_hiv")]:
        data = associations.copy()
        data[metric] = pd.to_numeric(data[metric], errors="coerce")
        matrix = data.pivot_table(index="Module", columns="Group", values=metric, aggfunc="first")
        fig, ax = plt.subplots()
        sns.heatmap(matrix, cmap="vlag", center=0, annot=True, fmt=".2f", linewidths=0.7, cbar_kws={"label": "Cohen's d"}, ax=ax)
        ax.set(xlabel="Clinical group", ylabel="Proteomic module")
        finish_figure(fig, out / f"Fig7g_module_effects_{suffix}", 8, 5.5)


def main():
    args, repo, out = setup("Proteomic NMF module activity", "16_proteomics_nmf")
    definitions = repo.table(33, "ST33a")
    loadings = repo.table(33, "ST33b")
    activity = repo.table(33, "ST33c")
    associations = repo.table(34, "ST34a")
    refit_scores, refit_loadings, reconstruction_error = refit_source_nmf(repo, args.seed)
    modules = [column for column in activity.columns if str(column).startswith("M")]
    for column in modules:
        activity[column] = pd.to_numeric(activity[column], errors="coerce")
    group_activity = activity.groupby("group", observed=True)[modules].median()
    write_table(definitions, out / "nmf_module_definitions.tsv")
    write_table(loadings, out / "nmf_protein_loadings.tsv")
    write_table(associations, out / "nmf_disease_associations.tsv")
    write_table(refit_scores, out / "source_data4_nmf_sample_scores.tsv")
    write_table(refit_loadings, out / "source_data4_nmf_protein_loadings.tsv")
    write_table(pd.DataFrame([{"components": 8, "reconstruction_error": reconstruction_error, "random_seed": args.seed}]), out / "source_data4_nmf_fit.tsv")
    plot_fig7d(activity, definitions, modules, out)
    plot_fig7g(associations, out)
    fig, ax = plt.subplots()
    image = heatmap(ax, group_activity.to_numpy(), modules, group_activity.index, cmap="viridis")
    fig.colorbar(image, ax=ax, label="Median module activity")
    ax.set(xlabel="Proteomic module", ylabel="Clinical group")
    finish_figure(fig, out / "EDFig8_proteomic_module_group_medians")


if __name__ == "__main__":
    main()
