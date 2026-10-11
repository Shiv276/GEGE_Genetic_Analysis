"""GEGE3888: compare GMMAT GRM-only, experiment, PCA, and combined models.

Run after gwas_infection_all_models.R has written four gwas_*.tsv files.
Uses the SAME SNPs across models for fair diagnostics and plots.
"""
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import chi2
from statsmodels.stats.multitest import multipletests

from config import PHENOTYPE_DIR

GWAS_DIR = PHENOTYPE_DIR.parent / "gwas"

MODELS = {
    "GRM only": "grm_only",
    "GRM + experiment": "grm_experiment",
    "GRM + 3 PCs": "grm_pc3",
    "GRM + experiment + 3 PCs": "grm_experiment_pc3",
}


def load_results():
    """Load four GMMAT scans and ensure the very same SNPs are compared."""
    data = {}
    for label, name in MODELS.items():
        filename = GWAS_DIR / f"gwas_{name}.tsv"
        df = pd.read_csv(filename, sep="\t")

        needed = {"SNP", "CHR", "POS", "PVAL"}
        if not needed.issubset(df.columns):
            raise ValueError(f"{filename.name} lacks {needed - set(df.columns)}")
        if df["SNP"].duplicated().any():
            raise ValueError(f"Duplicate SNP IDs in {filename.name}")
        if df["PVAL"].isna().any() or not df["PVAL"].between(0, 1, inclusive="right").all():
            raise ValueError(f"Missing or invalid p-values in {filename.name}")

        df["CHR"] = pd.to_numeric(df["CHR"], errors="raise").astype(int)
        df["POS"] = pd.to_numeric(df["POS"], errors="raise").astype(int)
        data[label] = df.set_index("SNP", drop=False)

    first = next(iter(data.values()))
    snp_ids = first.index
    for label, df in data.items():
        if len(df) != len(first) or not df.index.isin(snp_ids).all():
            raise ValueError(f"Different tested SNP sets for {label}; verify GWAS inputs first")
        df = df.loc[snp_ids].copy()  # SNP order identical across models
        if not (df["CHR"].to_numpy() == first["CHR"].to_numpy()).all():
            raise ValueError(f"Chromosome mismatch in {label}")
        if not (df["POS"].to_numpy() == first["POS"].to_numpy()).all():
            raise ValueError(f"Position mismatch in {label}")
        data[label] = df.reset_index(drop=True)

    print(f"Confirmed: identical set of {len(first):,} SNPs in all four models.")
    return data


def add_adjusted_pvalues(df):
    """Calculate model-specific Benjamini-Hochberg and Bonferroni corrections."""
    df = df.copy()
    p = df["PVAL"].to_numpy(dtype=float)
    df["FDR_BH"] = multipletests(p, method="fdr_bh")[1]
    df["Bonferroni"] = np.minimum(p * len(p), 1.0)
    df["minus_log10_p"] = -np.log10(p)
    return df


def lambda_gc(p):
    """Ratio of median observed chi-square statistic to the null median (1 df)."""
    return np.median(chi2.isf(p, df=1)) / chi2.ppf(0.5, df=1)


def make_summary(data):
    """Report model calibration and SNP-level multiple testing results."""
    rows = []
    for label, df in data.items():
        p = df["PVAL"].to_numpy()
        rows.append({
            "model": label,
            "snps": len(df),
            "lambda_gc": lambda_gc(p),
            "min_p": p.min(),
            "nominal_p_below_0.05": int((p < 0.05).sum()),
            "fdr_below_0.05": int((df["FDR_BH"] < 0.05).sum()),
            "bonferroni_below_0.05": int((df["Bonferroni"] < 0.05).sum()),
        })
    summary = pd.DataFrame(rows)
    summary.to_csv(GWAS_DIR / "gwas_model_comparison.csv", index=False)
    print("\nMODEL COMPARISON\n", summary.to_string(index=False))
    return summary


def plot_qq_comparison(data):
    """Overlay null-expected versus observed p-value quantiles for four models."""
    fig, ax = plt.subplots(figsize=(7, 6))
    max_observed = 0
    for label, df in data.items():
        p = np.sort(df["PVAL"].to_numpy())[::-1]  # large to small
        n = len(p)
        expected = -np.log10((np.arange(n, 0, -1) - 0.5) / n)
        observed = -np.log10(p)
        max_observed = max(max_observed, observed.max())
        ax.plot(expected, observed, linewidth=1.5,
                label=f"{label} (λ={lambda_gc(p):.3f})")

    limit = 1.05 * max(max_observed, expected.max())
    ax.plot([0, limit], [0, limit], linestyle="--", color="grey", label="Null expectation")
    ax.set(xlim=(0, limit), ylim=(0, limit),
           xlabel="Expected −log10(p)", ylabel="Observed −log10(p)",
           title="Bd infection GWAS — QQ comparison")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(GWAS_DIR / "gwas_qq_comparison.png", dpi=250)
    fig.savefig(GWAS_DIR / "gwas_qq_comparison.pdf")
    return fig


def get_genomic_positions(df):
    """Calculate chromosome-offset x coordinates from physical SNP positions."""
    loci = df[["CHR", "POS"]].copy()
    chr_lengths = loci.groupby("CHR")["POS"].max().sort_index()
    offsets = chr_lengths.cumsum().shift(fill_value=0)
    loci["x"] = loci["POS"] + loci["CHR"].map(offsets)
    ticks = loci.groupby("CHR")["x"].agg(["min", "max"]).mean(axis=1)
    return loci, ticks


def plot_manhattan_comparison(data):
    """Four Manhattan panels on identical chromosome and significance scales."""
    first = next(iter(data.values()))
    loci, ticks = get_genomic_positions(first)
    n = len(first)
    threshold = -np.log10(0.05 / n)
    highest = max(threshold, *(df["minus_log10_p"].max() for df in data.values()))

    fig, axes = plt.subplots(2, 2, figsize=(15, 9), sharex=True, sharey=True)
    for ax, (label, df) in zip(axes.flat, data.items()):
        for chrom in sorted(loci["CHR"].unique()):
            subset = loci["CHR"] == chrom
            ax.scatter(loci.loc[subset, "x"], df.loc[subset, "minus_log10_p"],
                       s=3, alpha=0.65, rasterized=True)
        ax.axhline(threshold, color="red", linestyle="--", linewidth=1)
        ax.set_title(label)
        ax.set_ylim(0, highest * 1.07)
        ax.set_xticks(ticks.values)
        ax.set_xticklabels(ticks.index.astype(str))
        ax.set_xlabel("Chromosome")
        ax.set_ylabel("−log10(p)")

    fig.suptitle("Bd infection GWAS — four-model comparison")
    fig.tight_layout()
    fig.savefig(GWAS_DIR / "gwas_manhattan_comparison.png", dpi=250)
    fig.savefig(GWAS_DIR / "gwas_manhattan_comparison.pdf")
    return fig


def save_top_snps(data, n=10):
    """Compare p-values of leading SNPs across all four association models."""
    # Include the top N from EACH model so changed rankings are visible.
    top_ids = set()
    for df in data.values():
        top_ids.update(df.nsmallest(n, "PVAL")["SNP"])

    first = next(iter(data.values()))
    comparison = first.loc[first["SNP"].isin(top_ids),
                           ["SNP", "CHR", "POS"]].copy().set_index("SNP")
    for label, df in data.items():
        comparison[label] = df.set_index("SNP").loc[comparison.index, "PVAL"]
    comparison["minimum_p_across_models"] = comparison[list(data)].min(axis=1)
    comparison = comparison.sort_values("minimum_p_across_models")
    comparison.to_csv(GWAS_DIR / "gwas_top_snps_across_models.csv")
    print("\nTOP ASSOCIATIONS ACROSS MODELS\n", comparison.head(15).to_string())


def main():
    GWAS_DIR.mkdir(parents=True, exist_ok=True)
    data = {label: add_adjusted_pvalues(df)
            for label, df in load_results().items()}
    for label, name in MODELS.items():
        data[label].to_csv(GWAS_DIR / f"gwas_{name}_processed.csv", index=False)

    make_summary(data)
    save_top_snps(data)
    plot_qq_comparison(data)
    plot_manhattan_comparison(data)
    plt.show()


if __name__ == "__main__":
    main()
    import pandas as pd
    from config import GWAS_DIR

    summary = pd.read_csv(
        GWAS_DIR / "gwas_model_comparison.csv"
    )

    print(summary.to_string(index=False))
