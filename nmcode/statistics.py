import numpy as np
import pandas as pd
import warnings
from scipy.stats import kruskal, mannwhitneyu, pearsonr, spearmanr, ttest_ind
from statsmodels.stats.multitest import multipletests


def benjamini_hochberg(values):
    series = pd.Series(values, dtype=float)
    valid = series.notna()
    result = pd.Series(np.nan, index=series.index)
    if valid.any():
        result.loc[valid] = multipletests(series.loc[valid], method="fdr_bh")[1]
    return result.to_numpy()


def donor_composition(data, donor, group, feature):
    counts = data.groupby([donor, group, feature], observed=True).size().rename("cells").reset_index()
    donors = data[[donor, group]].drop_duplicates().assign(_key=1)
    features = data[[feature]].drop_duplicates().assign(_key=1)
    result = donors.merge(features, on="_key").drop(columns="_key")
    result = result.merge(counts, on=[donor, group, feature], how="left")
    result["cells"] = result["cells"].fillna(0).astype(int)
    totals = data.groupby([donor, group], observed=True).size().rename("total_cells").reset_index()
    result = result.merge(totals, on=[donor, group], how="left")
    result["proportion"] = result["cells"] / result["total_cells"]
    return result


def kruskal_by_feature(data, feature, group, value):
    rows = []
    for label, frame in data.groupby(feature, observed=True):
        samples = [part[value].dropna().to_numpy() for _, part in frame.groupby(group, observed=True)]
        samples = [sample for sample in samples if len(sample)]
        statistic, p_value = kruskal(*samples) if len(samples) > 1 else (np.nan, np.nan)
        rows.append({feature: label, "statistic": statistic, "p_value": p_value})
    result = pd.DataFrame(rows)
    result["q_value"] = benjamini_hochberg(result["p_value"])
    return result


def mann_whitney_by_feature(data, feature, group, value, reference):
    rows = []
    for label, frame in data.groupby(feature, observed=True):
        baseline = frame.loc[frame[group] == reference, value].dropna()
        for comparison, current in frame.loc[frame[group] != reference].groupby(group, observed=True):
            observed = current[value].dropna()
            statistic, p_value = mannwhitneyu(observed, baseline, alternative="two-sided") if len(observed) and len(baseline) else (np.nan, np.nan)
            rows.append({feature: label, "comparison": comparison, "reference": reference, "statistic": statistic, "p_value": p_value})
    result = pd.DataFrame(rows)
    result["q_value"] = benjamini_hochberg(result["p_value"])
    return result


def spearman_table(data, pairs):
    rows = []
    for left, right in pairs:
        frame = data[[left, right]].apply(pd.to_numeric, errors="coerce").dropna()
        rho, p_value = spearmanr(frame[left], frame[right]) if len(frame) > 2 else (np.nan, np.nan)
        rows.append({"measure_1": left, "measure_2": right, "n": len(frame), "rho": rho, "p_value": p_value})
    result = pd.DataFrame(rows)
    result["q_value"] = benjamini_hochberg(result["p_value"])
    return result


def welch_differential(matrix, groups, reference="Control"):
    matrix = matrix.apply(pd.to_numeric, errors="coerce")
    groups = pd.Series(groups, index=matrix.index).astype(str)
    rows = []
    baseline = matrix.loc[groups.eq(reference)]
    for comparison in groups[~groups.eq(reference)].dropna().unique():
        current = matrix.loc[groups.eq(comparison)]
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            statistic, p_value = ttest_ind(current, baseline, axis=0, equal_var=False, nan_policy="omit")
        mean_current = current.mean(axis=0)
        mean_reference = baseline.mean(axis=0)
        n_current = current.notna().sum(axis=0)
        n_reference = baseline.notna().sum(axis=0)
        variance = ((n_current - 1) * current.var(axis=0) + (n_reference - 1) * baseline.var(axis=0)) / (n_current + n_reference - 2)
        effect = (mean_current - mean_reference) / np.sqrt(variance.replace(0, np.nan))
        frame = pd.DataFrame({
            "feature": matrix.columns,
            "comparison": comparison,
            "reference": reference,
            "n_comparison": n_current.to_numpy(),
            "n_reference": n_reference.to_numpy(),
            "mean_difference": (mean_current - mean_reference).to_numpy(),
            "cohens_d": effect.to_numpy(),
            "welch_t": statistic,
            "p_value": p_value,
        })
        frame["q_value"] = benjamini_hochberg(frame["p_value"])
        rows.append(frame)
    return pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()


def partial_spearman_table(data, pairs, covariates=()):
    rows = []
    for left, right in pairs:
        columns = [left, right, *covariates]
        frame = data[columns].dropna().copy()
        if len(frame) < 4:
            rows.append({"measure_1": left, "measure_2": right, "n": len(frame), "rho": np.nan, "p_value": np.nan})
            continue
        ranked = frame[[left, right]].rank(method="average")
        if covariates:
            design = pd.get_dummies(frame[list(covariates)], drop_first=True, dtype=float)
            design = np.column_stack([np.ones(len(frame)), design.to_numpy(dtype=float)])
            residual_left = ranked[left].to_numpy() - design @ np.linalg.lstsq(design, ranked[left].to_numpy(), rcond=None)[0]
            residual_right = ranked[right].to_numpy() - design @ np.linalg.lstsq(design, ranked[right].to_numpy(), rcond=None)[0]
            rho, p_value = pearsonr(residual_left, residual_right)
        else:
            rho, p_value = spearmanr(ranked[left], ranked[right])
        rows.append({"measure_1": left, "measure_2": right, "n": len(frame), "rho": rho, "p_value": p_value})
    result = pd.DataFrame(rows)
    result["q_value"] = benjamini_hochberg(result["p_value"])
    return result
