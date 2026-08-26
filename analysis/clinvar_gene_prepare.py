#!/usr/bin/env python3
"""ClinVar gene-level source tables and sequence-confounder controls.

This script uses the ClinVar cohort to create gene-level summaries for
downstream codon-contribution analyses and to test whether CaLM signal
persists after controlling explicit sequence features.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy.stats import rankdata
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold


INPUT = Path("Results/Revision/len1022_model_control/len1022_model_control_score_table_complete_cases.csv")
GENE_FASTA_DIR = Path("/Users/cassie/Desktop/Gene")
OUT_DIR = Path("Results/Revision/len1022_core_clinvar")

BASES = ("A", "C", "G", "T")
TRANSITIONS = {("A", "G"), ("G", "A"), ("C", "T"), ("T", "C")}
CODON_TABLE = {
    "TTT": "F",
    "TTC": "F",
    "TTA": "L",
    "TTG": "L",
    "TCT": "S",
    "TCC": "S",
    "TCA": "S",
    "TCG": "S",
    "TAT": "Y",
    "TAC": "Y",
    "TAA": "*",
    "TAG": "*",
    "TGT": "C",
    "TGC": "C",
    "TGA": "*",
    "TGG": "W",
    "CTT": "L",
    "CTC": "L",
    "CTA": "L",
    "CTG": "L",
    "CCT": "P",
    "CCC": "P",
    "CCA": "P",
    "CCG": "P",
    "CAT": "H",
    "CAC": "H",
    "CAA": "Q",
    "CAG": "Q",
    "CGT": "R",
    "CGC": "R",
    "CGA": "R",
    "CGG": "R",
    "ATT": "I",
    "ATC": "I",
    "ATA": "I",
    "ATG": "M",
    "ACT": "T",
    "ACC": "T",
    "ACA": "T",
    "ACG": "T",
    "AAT": "N",
    "AAC": "N",
    "AAA": "K",
    "AAG": "K",
    "AGT": "S",
    "AGC": "S",
    "AGA": "R",
    "AGG": "R",
    "GTT": "V",
    "GTC": "V",
    "GTA": "V",
    "GTG": "V",
    "GCT": "A",
    "GCC": "A",
    "GCA": "A",
    "GCG": "A",
    "GAT": "D",
    "GAC": "D",
    "GAA": "E",
    "GAG": "E",
    "GGT": "G",
    "GGC": "G",
    "GGA": "G",
    "GGG": "G",
}


def auc_rank_sum(y: np.ndarray, score: np.ndarray) -> float:
    y_bool = y.astype(bool)
    n_pos = int(y_bool.sum())
    n_neg = int(len(y_bool) - n_pos)
    if n_pos == 0 or n_neg == 0:
        raise ValueError("AUROC requires at least one positive and one negative label")
    ranks = rankdata(score, method="average")
    pos_rank_sum = float(ranks[y_bool].sum())
    return (pos_rank_sum - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)


def best_weight(y: np.ndarray, a: np.ndarray, b: np.ndarray) -> tuple[float, float]:
    best_auc = -np.inf
    best_w = np.nan
    for w in np.linspace(0, 1, 101):
        score = (1 - w) * a + w * b
        auc = auc_rank_sum(y, score)
        if auc > best_auc:
            best_auc = auc
            best_w = w
    return float(best_w), float(best_auc)


def nested_pair_auc(
    y: np.ndarray,
    a: np.ndarray,
    b: np.ndarray,
    *,
    n_splits: int,
    seed: int,
) -> tuple[float, float, float, float, int, list[dict[str, float | int]]]:
    """Select pair weights on inner training folds and evaluate held-out variants."""
    splitter = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=seed)
    heldout_score = np.full(len(y), np.nan, dtype=float)
    weights: list[float] = []
    fold_rows: list[dict[str, float | int]] = []

    for fold, (train_idx, test_idx) in enumerate(splitter.split(a, y), start=1):
        weight, train_auc = best_weight(y[train_idx], a[train_idx], b[train_idx])
        test_score = (1 - weight) * a[test_idx] + weight * b[test_idx]
        heldout_score[test_idx] = test_score
        weights.append(weight)
        fold_rows.append(
            {
                "fold": fold,
                "weight_second": weight,
                "train_auc_at_selected_weight": train_auc,
                "test_auc_at_selected_weight": auc_rank_sum(y[test_idx], test_score),
                "n_train": int(len(train_idx)),
                "n_test": int(len(test_idx)),
                "n_test_pathogenic": int(y[test_idx].sum()),
                "n_test_benign": int(len(test_idx) - y[test_idx].sum()),
            }
        )

    if np.isnan(heldout_score).any():
        raise RuntimeError("Nested CV did not score every held-out variant")
    return (
        float(np.mean(weights)),
        float(np.std(weights, ddof=1)) if len(weights) > 1 else 0.0,
        float(auc_rank_sum(y, heldout_score)),
        float(np.mean([row["train_auc_at_selected_weight"] for row in fold_rows])),
        int(len(fold_rows)),
        fold_rows,
    )


def per_gene_nested_summary(
    df: pd.DataFrame,
    min_pos: int = 5,
    min_neg: int = 5,
    max_splits: int = 5,
    seed: int = 16,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    rows: list[dict[str, float | int | str]] = []
    fold_rows: list[dict[str, float | int | str]] = []
    pair_specs = [
        ("ESM-2 150M + 650M", "esm2_150m_score", "esm2_650m_score"),
        ("ESM-2 650M + ESM-1b 650M", "esm2_650m_score", "esm1b_650m_score"),
        ("ESM-2 650M + CaLM", "esm2_650m_score", "calm_score"),
        ("ESM-2 150M + CaLM", "esm2_150m_score", "calm_score"),
        ("ESM-1b 650M + CaLM", "esm1b_650m_score", "calm_score"),
    ]

    for gene, gene_df in df.groupby("Gene_gene", sort=True):
        y = gene_df["label"].to_numpy()
        n_pos = int(y.sum())
        n_neg = int(len(y) - n_pos)
        if n_pos < min_pos or n_neg < min_neg:
            continue

        s150 = gene_df["esm2_150m_score"].to_numpy(dtype=float)
        s650 = gene_df["esm2_650m_score"].to_numpy(dtype=float)
        s1b = gene_df["esm1b_650m_score"].to_numpy(dtype=float)
        calm = gene_df["calm_score"].to_numpy(dtype=float)
        score_map = {
            "esm2_150m_score": s150,
            "esm2_650m_score": s650,
            "esm1b_650m_score": s1b,
            "calm_score": calm,
        }

        row: dict[str, float | int | str] = {
            "gene": gene,
            "n": int(len(gene_df)),
            "n_pathogenic": n_pos,
            "n_benign": n_neg,
            "protein_length": int(gene_df["protein_length"].iloc[0]),
            "nested_cv_splits": int(min(max_splits, n_pos, n_neg)),
            "ESM-2 150M AUROC": roc_auc_score(y, s150),
            "ESM-2 650M AUROC": roc_auc_score(y, s650),
            "ESM-1b 650M AUROC": roc_auc_score(y, s1b),
            "CaLM AUROC": roc_auc_score(y, calm),
        }

        for offset, (label, first_col, second_col) in enumerate(pair_specs):
            weight, weight_sd, auc, train_auc, n_folds, pair_folds = nested_pair_auc(
                y,
                score_map[first_col],
                score_map[second_col],
                n_splits=int(row["nested_cv_splits"]),
                seed=seed + 1009 * offset + len(rows),
            )
            row[f"{label} AUROC"] = auc
            row[f"{label} weight_second"] = weight
            row[f"{label} weight_second_sd"] = weight_sd
            row[f"{label} mean_train_auc_at_selected_weights"] = train_auc
            row[f"{label} n_outer_folds"] = n_folds
            for fold_row in pair_folds:
                fold_rows.append({"gene": gene, "model": label, **fold_row})

        row["650M_minus_150M"] = row["ESM-2 650M AUROC"] - row["ESM-2 150M AUROC"]
        row["150M650M_minus_650M"] = row["ESM-2 150M + 650M AUROC"] - row["ESM-2 650M AUROC"]
        row["650M1b_minus_650M"] = row["ESM-2 650M + ESM-1b 650M AUROC"] - row["ESM-2 650M AUROC"]
        row["650MCaLM_minus_650M"] = row["ESM-2 650M + CaLM AUROC"] - row["ESM-2 650M AUROC"]
        row["1bCaLM_minus_1b"] = row["ESM-1b 650M + CaLM AUROC"] - row["ESM-1b 650M AUROC"]
        row["650MCaLM_minus_150M650M"] = row["ESM-2 650M + CaLM AUROC"] - row["ESM-2 150M + 650M AUROC"]
        row["650MCaLM_minus_650M1b"] = row["ESM-2 650M + CaLM AUROC"] - row["ESM-2 650M + ESM-1b 650M AUROC"]
        rows.append(row)

    return pd.DataFrame(rows), pd.DataFrame(fold_rows)


def fit_regression(df: pd.DataFrame, x_col: str, y_col: str, label: str, weighted: bool) -> dict[str, float | str | int]:
    reg_df = df.dropna(subset=[x_col, y_col, "n_pathogenic", "n_benign"]).copy()
    x = reg_df[x_col].to_numpy(dtype=float)
    y = reg_df[y_col].to_numpy(dtype=float)
    X = sm.add_constant(x)
    if weighted:
        weights = (reg_df["n_pathogenic"] * reg_df["n_benign"]) / (reg_df["n_pathogenic"] + reg_df["n_benign"])
        model = sm.WLS(y, X, weights=weights.to_numpy(dtype=float)).fit(cov_type="HC3")
        weight_label = "n_pathogenic*n_benign/(n_pathogenic+n_benign)"
        method = "WLS_HC3"
    else:
        model = sm.OLS(y, X).fit(cov_type="HC3")
        weight_label = ""
        method = "OLS_HC3"
    ci_low, ci_high = [float(v) for v in model.conf_int(alpha=0.05)[1]]
    return {
        "analysis": label,
        "method": method,
        "weighted": int(weighted),
        "n_genes": int(len(reg_df)),
        "x_column": x_col,
        "y_column": y_col,
        "weight_definition": weight_label,
        "intercept": float(model.params[0]),
        "slope": float(model.params[1]),
        "slope_95ci_low": ci_low,
        "slope_95ci_high": ci_high,
        "p_value": float(model.pvalues[1]),
        "r_squared": float(model.rsquared),
    }


def read_fasta(path: Path) -> str | None:
    if not path.exists():
        return None
    parts: list[str] = []
    with path.open() as handle:
        for line in handle:
            line = line.strip()
            if line and not line.startswith(">"):
                parts.append(line.upper().replace("U", "T"))
    return "".join(parts) or None


def gene_composition(gene: str) -> dict[str, float]:
    seq = read_fasta(GENE_FASTA_DIR / f"{gene}.fasta")
    if seq is None:
        return {"gene_gc": np.nan, "gene_cpg_density": np.nan}
    clean = "".join(base for base in seq if base in BASES)
    if not clean:
        return {"gene_gc": np.nan, "gene_cpg_density": np.nan}
    return {
        "gene_gc": (clean.count("G") + clean.count("C")) / len(clean),
        "gene_cpg_density": clean.count("CG") / max(len(clean) - 1, 1),
    }


def mutation_features(ref: str, mut: str) -> dict[str, float]:
    ref = str(ref).upper().replace("U", "T")
    mut = str(mut).upper().replace("U", "T")
    out = {
        "valid_snv_codon": 0.0,
        "ref_codon_gc": np.nan,
        "mut_codon_gc": np.nan,
        "delta_codon_gc": np.nan,
        "ref_cpg": np.nan,
        "mut_cpg": np.nan,
        "delta_cpg": np.nan,
        "is_transition": np.nan,
        "is_transversion": np.nan,
        "mutated_codon_position": np.nan,
        "local_degeneracy": np.nan,
    }
    if len(ref) != 3 or len(mut) != 3 or any(base not in BASES for base in ref + mut):
        return out
    diffs = [idx for idx, (a, b) in enumerate(zip(ref, mut)) if a != b]
    if len(diffs) != 1:
        return out
    pos = diffs[0]
    ref_aa = CODON_TABLE.get(ref)
    if ref_aa is None or ref_aa == "*":
        return out

    syn_count = 0
    for base in BASES:
        candidate = ref[:pos] + base + ref[pos + 1 :]
        if CODON_TABLE.get(candidate) == ref_aa:
            syn_count += 1

    ref_gc = sum(base in {"G", "C"} for base in ref)
    mut_gc = sum(base in {"G", "C"} for base in mut)
    ref_cpg = int("CG" in ref)
    mut_cpg = int("CG" in mut)
    transition = int((ref[pos], mut[pos]) in TRANSITIONS)
    out.update(
        {
            "valid_snv_codon": 1.0,
            "ref_codon_gc": ref_gc / 3,
            "mut_codon_gc": mut_gc / 3,
            "delta_codon_gc": (mut_gc - ref_gc) / 3,
            "ref_cpg": float(ref_cpg),
            "mut_cpg": float(mut_cpg),
            "delta_cpg": float(mut_cpg - ref_cpg),
            "is_transition": float(transition),
            "is_transversion": float(1 - transition),
            "mutated_codon_position": float(pos + 1),
            "local_degeneracy": float(syn_count),
        }
    )
    return out


def zscore(series: pd.Series) -> pd.Series:
    std = series.std(ddof=0)
    if std == 0 or np.isnan(std):
        return series * np.nan
    return (series - series.mean()) / std


def design_matrix(df: pd.DataFrame, columns: list[str], add_constant: bool = True) -> pd.DataFrame:
    x = df[columns].copy()
    for col in ["local_degeneracy", "mutated_codon_position"]:
        if col in x:
            dummies = pd.get_dummies(x[col].astype("Int64").astype(str), prefix=col, drop_first=True)
            x = pd.concat([x.drop(columns=[col]), dummies], axis=1)
    x = x.astype(float)
    if add_constant:
        x = sm.add_constant(x, has_constant="add")
    return x


def fit_glm_clustered(
    df: pd.DataFrame,
    predictors: list[str],
    analysis: str,
    focal_predictor: str | None,
) -> pd.Series:
    model_df = df.dropna(subset=["label", "Gene_gene", *predictors]).copy()
    y = model_df["label"].astype(float)
    x = design_matrix(model_df, predictors)
    model = sm.GLM(y, x, family=sm.families.Binomial()).fit(
        cov_type="cluster",
        cov_kwds={"groups": model_df["Gene_gene"]},
        maxiter=200,
    )
    pred = pd.Series(model.predict(x), index=model_df.index)
    out = {
        "analysis": analysis,
        "n_variants": int(len(model_df)),
        "n_genes": int(model_df["Gene_gene"].nunique()),
        "focal_predictor": focal_predictor or "",
        "coef": np.nan,
        "coef_95ci_low": np.nan,
        "coef_95ci_high": np.nan,
        "odds_ratio": np.nan,
        "or_95ci_low": np.nan,
        "or_95ci_high": np.nan,
        "p_value": np.nan,
        "auroc": float(roc_auc_score(y, pred)),
        "pseudo_r2_mcfadden": float(1 - model.llf / model.llnull),
    }
    if focal_predictor is not None:
        ci_low, ci_high = model.conf_int().loc[focal_predictor].tolist()
        coef = float(model.params[focal_predictor])
        out.update(
            {
                "coef": coef,
                "coef_95ci_low": float(ci_low),
                "coef_95ci_high": float(ci_high),
                "odds_ratio": float(np.exp(coef)),
                "or_95ci_low": float(np.exp(ci_low)),
                "or_95ci_high": float(np.exp(ci_high)),
                "p_value": float(model.pvalues[focal_predictor]),
            }
        )
    return pd.Series(out)


def sequence_confounder_control(df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    feature_df = pd.DataFrame([mutation_features(ref, mut) for ref, mut in zip(df["Ref_gene"], df["Mut_gene"])])
    valid = pd.concat([df.reset_index(drop=True), feature_df], axis=1)
    gene_comp = pd.DataFrame.from_dict(
        {gene: gene_composition(gene) for gene in sorted(valid["Gene_gene"].dropna().unique())},
        orient="index",
    )
    gene_comp.index.name = "Gene_gene"
    valid = valid.merge(gene_comp.reset_index(), on="Gene_gene", how="left")
    valid = valid[valid["valid_snv_codon"] == 1].copy()
    valid = valid.dropna(subset=["calm_score", "esm2_650m_score", "gene_gc", "gene_cpg_density"])
    if valid.empty:
        summary = pd.DataFrame(
            [
                {
                    "analysis": "sequence_confounder_control_skipped",
                    "n_variants": 0,
                    "n_genes": 0,
                    "focal_predictor": "",
                    "coef": np.nan,
                    "coef_95ci_low": np.nan,
                    "coef_95ci_high": np.nan,
                    "odds_ratio": np.nan,
                    "or_95ci_low": np.nan,
                    "or_95ci_high": np.nan,
                    "p_value": np.nan,
                    "auroc": np.nan,
                    "pseudo_r2_mcfadden": np.nan,
                    "delta_auroc_vs_base": np.nan,
                    "skip_reason": "no variants with available gene FASTA composition covariates",
                }
            ]
        )
        residual_meta = pd.DataFrame(
            [
                {
                    "analysis": "calm_residualization_skipped",
                    "outcome": "z_calm_score",
                    "n_variants": 0,
                    "n_genes": 0,
                    "r_squared": np.nan,
                    "covariates": "",
                    "skip_reason": "no variants with available gene FASTA composition covariates",
                }
            ]
        )
        return summary, residual_meta, valid

    valid["z_calm_score"] = zscore(valid["calm_score"])
    valid["z_esm2_650m_score"] = zscore(valid["esm2_650m_score"])
    valid["z_gene_gc"] = zscore(valid["gene_gc"])
    valid["z_gene_cpg_density"] = zscore(valid["gene_cpg_density"])

    confounders = [
        "z_esm2_650m_score",
        "ref_codon_gc",
        "delta_codon_gc",
        "ref_cpg",
        "delta_cpg",
        "is_transition",
        "local_degeneracy",
        "mutated_codon_position",
        "z_gene_gc",
        "z_gene_cpg_density",
    ]
    residual_x = design_matrix(valid, confounders)
    residual_model = sm.OLS(valid["z_calm_score"], residual_x).fit()
    valid["calm_incremental_residual"] = residual_model.resid
    valid["z_calm_incremental_residual"] = zscore(valid["calm_incremental_residual"])

    models = [
        fit_glm_clustered(valid, confounders, "base_esm2_plus_sequence_controls", None),
        fit_glm_clustered(valid, confounders + ["z_calm_score"], "base_plus_calm", "z_calm_score"),
        fit_glm_clustered(
            valid,
            ["z_calm_incremental_residual"],
            "calm_incremental_residual_only",
            "z_calm_incremental_residual",
        ),
    ]
    summary = pd.DataFrame(models)
    base_auroc = float(summary.loc[summary["analysis"] == "base_esm2_plus_sequence_controls", "auroc"].iloc[0])
    summary["delta_auroc_vs_base"] = summary["auroc"] - base_auroc

    residual_meta = pd.DataFrame(
        [
            {
                "analysis": "calm_residualization",
                "outcome": "z_calm_score",
                "n_variants": int(len(valid)),
                "n_genes": int(valid["Gene_gene"].nunique()),
                "r_squared": float(residual_model.rsquared),
                "covariates": ";".join(confounders),
            }
        ]
    )

    keep_cols = [
        "label",
        "Gene_gene",
        "variant_id",
        "calm_score",
        "esm2_650m_score",
        "valid_snv_codon",
        "ref_codon_gc",
        "mut_codon_gc",
        "delta_codon_gc",
        "ref_cpg",
        "mut_cpg",
        "delta_cpg",
        "gene_gc",
        "gene_cpg_density",
        "is_transition",
        "is_transversion",
        "mutated_codon_position",
        "local_degeneracy",
        "z_calm_score",
        "z_esm2_650m_score",
        "z_calm_incremental_residual",
    ]
    return summary, residual_meta, valid[keep_cols]


def summarize_core_analyses(new_per_gene: pd.DataFrame, new_regs: pd.DataFrame, seq_summary: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, float | str | int]] = []
    rows.append(
        {
            "analysis": "per_gene_panel",
            "cohort": "len1022_complete",
            "n_genes": int(len(new_per_gene)),
            "mean_delta_delta": float(new_per_gene["650MCaLM_minus_150M650M"].mean()),
            "median_delta_delta": float(new_per_gene["650MCaLM_minus_150M650M"].median()),
            "fraction_delta_delta_positive": float((new_per_gene["650MCaLM_minus_150M650M"] > 0).mean()),
        }
    )
    for _, row in new_regs.iterrows():
        rows.append(
            {
                "analysis": row["analysis"],
                "cohort": "len1022_complete",
                "n_genes": int(row["n_genes"]),
                "slope": float(row["slope"]),
                "slope_95ci_low": float(row["slope_95ci_low"]),
                "slope_95ci_high": float(row["slope_95ci_high"]),
                "p_value": float(row["p_value"]),
                "r_squared": float(row["r_squared"]),
            }
        )
    for _, row in seq_summary.iterrows():
        if row["focal_predictor"]:
            rows.append(
                {
                    "analysis": row["analysis"],
                    "cohort": "len1022_complete",
                    "n_genes": int(row["n_genes"]),
                    "n_variants": int(row["n_variants"]),
                    "odds_ratio": float(row["odds_ratio"]),
                    "or_95ci_low": float(row["or_95ci_low"]),
                    "or_95ci_high": float(row["or_95ci_high"]),
                    "p_value": float(row["p_value"]),
                    "delta_auroc_vs_base": float(row["delta_auroc_vs_base"]),
                }
            )
    return pd.DataFrame(rows)


def add_gene_level_derived_columns(per_gene: pd.DataFrame) -> pd.DataFrame:
    out = per_gene.copy()
    out["codon_dependency_index"] = out["ESM-2 650M + CaLM weight_second"]
    out["delta_delta_auroc"] = out["650MCaLM_minus_150M650M"]
    out["independent_calm_minus_650m"] = out["CaLM AUROC"] - out["ESM-2 650M AUROC"]
    out["inverse_variance_weight"] = (
        out["n_pathogenic"] * out["n_benign"] / (out["n_pathogenic"] + out["n_benign"])
    )
    return out


def gene_delta_delta_regressions(per_gene: pd.DataFrame) -> pd.DataFrame:
    return pd.DataFrame(
        [
            fit_regression(
                per_gene,
                "codon_dependency_index",
                "delta_delta_auroc",
                "codon_dependency_predicts_cross_modal_gain",
                False,
            ),
            fit_regression(
                per_gene,
                "codon_dependency_index",
                "delta_delta_auroc",
                "weighted_codon_dependency_predicts_cross_modal_gain",
                True,
            ),
            fit_regression(
                per_gene,
                "independent_calm_minus_650m",
                "delta_delta_auroc",
                "independent_calm_signal_predicts_cross_modal_gain",
                False,
            ),
            fit_regression(
                per_gene,
                "independent_calm_minus_650m",
                "delta_delta_auroc",
                "weighted_independent_calm_signal_predicts_cross_modal_gain",
                True,
            ),
            fit_regression(
                per_gene,
                "CaLM AUROC",
                "delta_delta_auroc",
                "calm_standalone_predicts_cross_modal_gain",
                False,
            ),
            fit_regression(
                per_gene,
                "CaLM AUROC",
                "delta_delta_auroc",
                "weighted_calm_standalone_predicts_cross_modal_gain",
                True,
            ),
        ]
    )


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(INPUT)
    df = df.drop_duplicates(subset=["variant_id", "Gene_gene"]).copy()

    cohort_summary = pd.DataFrame(
        [
            {
                "cohort": "len1022_complete_cases",
                "n_variants": int(len(df)),
                "n_genes": int(df["Gene_gene"].nunique()),
                "n_pathogenic": int((df["label"] == 1).sum()),
                "n_benign": int((df["label"] == 0).sum()),
                "min_protein_length": int(df["protein_length"].min()),
                "max_protein_length": int(df["protein_length"].max()),
            }
        ]
    )
    cohort_summary.to_csv(OUT_DIR / "cohort_summary.csv", index=False)

    nested_per_gene, nested_folds = per_gene_nested_summary(df)
    nested_per_gene = add_gene_level_derived_columns(nested_per_gene)
    nested_per_gene.to_csv(OUT_DIR / "per_gene_nested_ensemble_summary.csv", index=False)
    nested_folds.to_csv(OUT_DIR / "per_gene_nested_cv_fold_metrics.csv", index=False)
    nested_regs = gene_delta_delta_regressions(nested_per_gene)
    nested_regs.to_csv(OUT_DIR / "per_gene_nested_delta_delta_regressions.csv", index=False)

    seq_summary, residual_meta, seq_features = sequence_confounder_control(df)
    seq_summary.to_csv(OUT_DIR / "calm_sequence_confounder_regression_summary.csv", index=False)
    residual_meta.to_csv(OUT_DIR / "calm_sequence_confounder_residualization_summary.csv", index=False)
    seq_features.to_csv(OUT_DIR / "calm_sequence_confounder_variant_features.csv", index=False)

    comparison = summarize_core_analyses(nested_per_gene, nested_regs, seq_summary)
    comparison.to_csv(OUT_DIR / "core_analysis_conclusion_check.csv", index=False)

    print("Cohort")
    print(cohort_summary.to_string(index=False))
    print("\nNested CV per-gene panel")
    print(
        nested_per_gene[
            [
                "gene",
                "n",
                "n_pathogenic",
                "n_benign",
                "nested_cv_splits",
                "ESM-2 650M AUROC",
                "CaLM AUROC",
                "ESM-2 150M + 650M AUROC",
                "ESM-2 650M + CaLM AUROC",
                "delta_delta_auroc",
                "codon_dependency_index",
                "independent_calm_minus_650m",
            ]
        ]
        .describe(include="all")
        .to_string()
    )
    print("\nNested CV regressions")
    print(nested_regs.to_string(index=False, float_format=lambda v: f"{v:.6g}"))
    print("\nSequence confounder control")
    print(seq_summary.to_string(index=False, float_format=lambda v: f"{v:.6g}"))
    print("\nResidualization")
    print(residual_meta.to_string(index=False, float_format=lambda v: f"{v:.6g}"))
    print(f"\nWrote outputs to {OUT_DIR}")


if __name__ == "__main__":
    main()
