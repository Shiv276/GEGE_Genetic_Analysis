"""GEGE3888 | Compare exploratory survival / early-Bd GWAS results.

Run after R/GMMAT. For each analysis: summary, BH FDR, Bonferroni,
QQ comparison and two-panel Manhattan figure. Analysis-specific SNP sets
may differ; models WITHIN each analysis must test the same SNPs.
"""

from pathlib import Path
import argparse
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import chi2
from statsmodels.stats.multitest import multipletests

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import config

ROOT = Path(config.PHENOTYPE_DIR).parent / "gwas" / "ai_validation_work"
LABELS = {
    "grm_only": "GRM only",
    "grm_experiment": "GRM + experiment",
    "grm_pc3": "GRM + 3 PCs",
    "grm_experiment_pc3": "GRM + experiment + 3 PCs",
}


def load_scans(folder):
    """Read all available GMMAT models, requiring matching SNPs within trait."""
    data = {}
    for name in LABELS:
        file = folder / f"gwas_{name}.tsv"
        if not file.is_file():
            continue
        df = pd.read_csv(file, sep="\t")
        if not {"SNP", "CHR", "POS", "PVAL"}.issubset(df.columns):
            raise ValueError(f"Missing SNP/CHR/POS/PVAL columns in {file}")
        if df["SNP"].duplicated().any():
            raise ValueError(f"Duplicate SNP names in {file}")
        df["CHR"] = pd.to_numeric(df["CHR"], errors="raise").astype(int)
        df["POS"] = pd.to_numeric(df["POS"], errors="raise").astype(int)
        df["PVAL"] = pd.to_numeric(df["PVAL"], errors="raise")
        if not np.isfinite(df["PVAL"]).all() or not df["PVAL"].between(0, 1, inclusive="right").all():
            raise ValueError(f"Invalid or missing p-values in {file}")
        data[name] = df.set_index("SNP", drop=False)
    if not data:
        raise FileNotFoundError(f"No GMMAT output files found in {folder}")
    if len(data) < 2:
        raise ValueError(f"Need >=2 models to compare in {folder}")

    # Join by marker ID, never assume different output tables have same order.
    first = next(iter(data.values()))
    snps = first.index
    for name, df in data.items():
        if len(df) != len(first) or set(df.index) != set(snps):
            raise ValueError(f"Different SNP sets between models in {folder}: {name}")
        df = df.loc[snps].copy()
        if not np.array_equal(df["CHR"].values, first["CHR"].values) or not np.array_equal(df["POS"].values, first["POS"].values):
            raise ValueError(f"Different SNP chromosome/position in {folder}: {name}")
        data[name] = df.reset_index(drop=True)
    return data


def add_corrections(df):
    """Correct multiple SNP tests within a single trait/model GWAS scan."""
    df = df.copy()
    p = df["PVAL"].to_numpy()
    df["FDR_BH"] = multipletests(p, method="fdr_bh")[1]
    df["BONF"] = np.minimum(p * len(p), 1)
    df["minus_log10_p"] = -np.log10(p)
    return df


def lambda_gc(p):
    """Median observed one-df chi-square statistic / theoretical null median."""
    return float(np.median(chi2.isf(p, df=1)) / chi2.ppf(0.5, df=1))


def qq_plot(name, models, folder):
    """Overlay null-expected vs observed p-value quantiles across models."""
    fig, ax = plt.subplots(figsize=(7, 6))
    lim = 0
    for key, df in models.items():
        p = np.sort(df["PVAL"].to_numpy())
        n = len(p)
        expected = -np.log10((np.arange(1, n + 1) - 0.5) / n)
        observed = -np.log10(p)
        ax.plot(expected, observed, linewidth=1.5,
                label=f"{LABELS[key]} (λ={lambda_gc(p):.3f})")
        lim = max(lim, float(expected.max()), float(observed.max()))
    lim *= 1.05
    ax.plot([0, lim], [0, lim], "--", color="grey", label="Null expectation")
    ax.set(xlim=(0, lim), ylim=(0, lim),
           xlabel="Expected −log10(p)", ylabel="Observed −log10(p)",
           title=f"{name.replace('_',' ').title()} — QQ comparison")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(folder / "qq_models.png", dpi=250)
    fig.savefig(folder / "qq_models.pdf")
    plt.close(fig)


def chromosome_positions(df):
    """Construct chromosome-respecting cumulative genomic x-coordinates."""
    spans = df.groupby("CHR")["POS"].max().sort_index()
    # Gap between chromosomes is for readability, not a genomic distance.
    gap = int(spans.max() * 0.02)
    offsets = {}
    offset = 0
    for chrom, chrom_size in spans.items():
        offsets[chrom] = offset
        offset += int(chrom_size) + gap
    x = df["POS"].to_numpy() + df["CHR"].map(offsets).to_numpy()
    ticks = [(offsets[c] + spans.loc[c] / 2, str(c)) for c in spans.index]
    return x, ticks


def manhattan_comparison(name, models, folder):
    """Display genome-wide peaks on comparable axes for all models of a trait."""
    cols = 2
    rows = int(np.ceil(len(models) / cols))
    fig, axs = plt.subplots(rows, cols, figsize=(14, 4.5 * rows),
                            sharex=True, sharey=True, squeeze=False)
    axs_flat = axs.ravel()
    y_max = max(df["minus_log10_p"].max() for df in models.values())
    thresholds = [-np.log10(0.05 / len(df)) for df in models.values()]
    y_lim = max(float(y_max), max(thresholds)) + 0.5

    for ax, (key, df) in zip(axs_flat, models.items()):
        x, ticks = chromosome_positions(df)
        for chrom, sub in df.groupby("CHR", sort=True):
            ix = sub.index.to_numpy()
            ax.scatter(x[ix], df.loc[ix, "minus_log10_p"], s=4, alpha=0.65)
        ax.axhline(-np.log10(0.05 / len(df)), linestyle="--", color="grey",
                   label="Bonferroni 0.05")
        ax.set(title=LABELS[key], ylim=(0, y_lim), ylabel="−log10(p)",
               xlabel="Chromosome")
        ax.set_xticks([tick for tick, _ in ticks])
        ax.set_xticklabels([lab for _, lab in ticks])
        ax.legend(loc="upper right", fontsize=8)

    for ax in axs_flat[len(models):]:
        ax.set_visible(False)
    fig.suptitle(f"{name.replace('_',' ').title()} — SNP associations")
    fig.tight_layout()
    fig.savefig(folder / "manhattan_models.png", dpi=250)
    fig.savefig(folder / "manhattan_models.pdf")
    plt.close(fig)


def top_snps(models, folder, n_top=15):
    """Track the same leading SNPs across adjusted models for signal stability."""
    primary = "grm_experiment" if "grm_experiment" in models else next(iter(models))
    chosen = models[primary].nsmallest(n_top, "PVAL")[["SNP", "CHR", "POS"]].copy()
    for key, df in models.items():
        chosen = chosen.merge(
            df[["SNP", "PVAL", "FDR_BH"]].rename(
                columns={"PVAL": f"p_{key}", "FDR_BH": f"fdr_{key}"}),
            on="SNP", how="left", validate="one_to_one"
        )
    chosen.to_csv(folder / "top_snps_across_models.csv", index=False)


def phenotype_diagnostic(folder, trait):
    """Plot outcome distribution; not evidence of SNP association by itself."""
    df = pd.read_csv(folder / "phenotype.tsv", sep="\t")
    fig, ax = plt.subplots(figsize=(7, 4))
    if trait == "early_log10_max":
        ax.hist(df[trait].to_numpy(), bins=30, edgecolor="white")
        ax.set(xlabel="log10(maximum weeks 1–6 zoospore equivalents + 1)",
               ylabel="Frogs", title="Early positive burden distribution")
    else:
        counts = df[trait].value_counts().reindex([0, 1], fill_value=0)
        ax.bar(["0", "1"], counts.values)
        ax.set(xlabel=f"{trait} (0/1)", ylabel="Frogs",
               title="Binary phenotype counts")
        for x, y in enumerate(counts.values):
            ax.text(x, y, str(y), ha="center", va="bottom")
    fig.tight_layout()
    fig.savefig(folder / "phenotype_distribution.png", dpi=250)
    plt.close(fig)


def main():
    argp = argparse.ArgumentParser()
    argp.add_argument("--show", action="store_true",
                      help="Open a preview of the first generated QQ image separately if desired")
    args = argp.parse_args()
    manifest = pd.read_csv(ROOT / "analysis_manifest.tsv", sep="\t")
    summary = []
    for _, info in manifest.iterrows():
        name, trait = info["analysis"], info["trait"]
        folder = ROOT / name
        models = {key: add_corrections(df) for key, df in load_scans(folder).items()}
        for key, df in models.items():
            df.to_csv(folder / f"processed_{key}.csv", index=False)
            p = df["PVAL"].to_numpy()
            summary.append({
                "analysis": name, "trait": trait,
                "model": key, "n_frogs": int(info["n_frogs"]),
                "n_snps": len(df), "lambda_gc": lambda_gc(p),
                "min_p": float(p.min()),
                "nominal_p_lt_0.05": int((p < 0.05).sum()),
                "FDR_BH_lt_0.05": int((df["FDR_BH"] < 0.05).sum()),
                "BONF_lt_0.05": int((df["BONF"] < 0.05).sum()),
            })
        qq_plot(name, models, folder)
        manhattan_comparison(name, models, folder)
        top_snps(models, folder)
        phenotype_diagnostic(folder, trait)
        print(f"\nCompleted figures and corrected p-values for {name}")

    result_summary = pd.DataFrame(summary)
    result_summary.to_csv(ROOT / "all_gwas_results_summary.csv", index=False)
    print("\nGWAS SUMMARY\n", result_summary.to_string(index=False))
    print("\nAll files saved to", ROOT)
    if args.show:
        from PIL import Image
        fig = plt.figure(figsize=(8, 6))
        plt.imshow(Image.open(ROOT / manifest.iloc[0]["analysis"] / "qq_models.png"))
        plt.axis("off")
        plt.show()


if __name__ == "__main__":
    main()
