
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from scipy.stats import chi2
from statsmodels.stats.multitest import multipletests

from config import GWAS_DIR


def load_gwas_results():
    """Load GMMAT results and calculate multiple-testing corrections."""

    df = pd.read_csv(GWAS_DIR / "gwas_grm_only.tsv", sep="\t")

    # Retain valid association p-values.
    df = df.loc[
        df["PVAL"].notna()
        & (df["PVAL"] > 0)
        & (df["PVAL"] <= 1)
    ].copy()

    # Benjamini-Hochberg false discovery rate correction.
    df["FDR"] = multipletests(
        df["PVAL"], method="fdr_bh"
    )[1]

    # Bonferroni-adjusted p-values.
    df["Bonferroni"] = np.minimum(
        df["PVAL"] * len(df), 1
    )

    df["neglog10_p"] = -np.log10(df["PVAL"])

    return df


def plot_qq(df):
    """Compare observed association p-values with the null expectation."""

    p = np.sort(df["PVAL"].to_numpy())
    n = len(p)

    expected = -np.log10(
        (np.arange(n, 0, -1) - 0.5) / n
    )

    observed = -np.log10(p[::-1])

    # Genomic inflation factor for one-degree-of-freedom tests.
    lambda_gc = (
        np.median(chi2.isf(p, df=1))
        / chi2.ppf(0.5, df=1)
    )

    fig, ax = plt.subplots(figsize=(6, 6))

    ax.scatter(expected, observed, s=8, alpha=0.5)

    limit = max(expected.max(), observed.max()) * 1.05

    ax.plot([0, limit], [0, limit], "--", color="grey")

    ax.set(
        xlabel="Expected -log10(p)",
        ylabel="Observed -log10(p)",
        title=f"GRM-only GWAS (λGC = {lambda_gc:.3f})"
    )

    ax.set_xlim(0, limit)
    ax.set_ylim(0, limit)

    plt.tight_layout()
    plt.show()

    return lambda_gc


def plot_manhattan(df):
    """Plot genome-wide SNP association results by chromosome."""

    df = df.copy()
    df["CHR"] = pd.to_numeric(df["CHR"])
    df["POS"] = pd.to_numeric(df["POS"])

    df = df.sort_values(["CHR", "POS"])
    df["plot_position"] = np.arange(len(df))

    fig, ax = plt.subplots(figsize=(12, 5))

    for chromosome, group in df.groupby("CHR"):

        ax.scatter(
            group["plot_position"],
            group["neglog10_p"],
            s=5,
            alpha=0.7
        )

    # Genome-wide Bonferroni significance threshold.
    threshold = -np.log10(0.05 / len(df))

    ax.axhline(
        threshold,
        linestyle="--",
        color="red",
        label="Bonferroni threshold"
    )

    ticks = df.groupby("CHR")["plot_position"].median()

    ax.set_xticks(ticks.values)
    ax.set_xticklabels(ticks.index.astype(str))

    ax.set(
        xlabel="Chromosome",
        ylabel="-log10(p)",
        title="Bd infection GWAS — GRM only"
    )

    ax.legend()
    plt.tight_layout()
    plt.show()


def main():

    df = load_gwas_results()

    print(f"SNPs tested: {len(df)}")
    print(f"Nominal p < 0.05: {(df['PVAL'] < 0.05).sum()}")
    print(f"FDR < 0.05: {(df['FDR'] < 0.05).sum()}")
    print(f"Bonferroni < 0.05: {(df['Bonferroni'] < 0.05).sum()}")

    print("\nTop 10 SNPs:")
    print(
        df.nsmallest(10, "PVAL")[
            ["CHR", "SNP", "POS", "PVAL", "FDR"]
        ]
    )

    lambda_gc = plot_qq(df)
    print(f"\nGenomic inflation factor: {lambda_gc:.3f}")

    plot_manhattan(df)

    df.to_csv(
        GWAS_DIR / "gwas_grm_only_processed.csv",
        index=False
    )


if __name__ == "__main__":
    main()
