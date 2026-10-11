"""GEGE3888 | Exploratory phenotype GWAS preparation.

Inputs (from existing project pipeline):
  data/processed/phenotypes/matched_phenotypes.csv
  data/processed/gwas/gwas_input_with_pcs.tsv
  data/processed/gwas/grm_892.csv and grm_892_ids.tsv
  data/processed/gwas/infection_qc.{bed,bim,fam}

Outputs are isolated in data/processed/gwas/ai_validation_work/.
No original source or completed GWAS file is modified.
"""

from pathlib import Path
import argparse
import subprocess
import sys

import numpy as np
import pandas as pd

# These scripts live one folder below the existing Scripts/config.py.
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import config

PHENOTYPE_DIR = Path(config.PHENOTYPE_DIR)
GWAS_ROOT = PHENOTYPE_DIR.parent / "gwas"
OUTPUT_ROOT = GWAS_ROOT / "ai_validation_work"
PLINK2 = str(getattr(config, "PLINK2", "plink2"))
ITS_COPIES_PER_ZOOSPORE = 86  # Wastepoint isolate; source R script
WEEKS = [f"Bd_wk{i}" for i in range(1, 7)]
BASE_COLUMNS = ["FID", "IID", "experiment", "geno_zoo", "PC1", "PC2", "PC3"]


def load_and_verify_inputs():
    """Match raw repeated-Bd records to the exact earlier 892-frog GWAS order."""
    base = pd.read_csv(
        GWAS_ROOT / "gwas_input_with_pcs.tsv", sep="\t",
        dtype={"FID": str, "IID": str},
    )
    raw = pd.read_csv(
        PHENOTYPE_DIR / "matched_phenotypes.csv", dtype={"IID": str},
    )
    ids = pd.read_csv(
        GWAS_ROOT / "grm_892_ids.tsv", sep="\t", dtype={"IID": str, "FID": str}
    )
    K = np.loadtxt(GWAS_ROOT / "grm_892.csv", delimiter=",")

    if len(base) != 892 or base["IID"].duplicated().any():
        raise ValueError("The original GWAS input must contain 892 unique frogs.")
    if not base["IID"].equals(ids["IID"]):
        raise ValueError("The 892-frog phenotype and GRM ID orders do not match.")
    if K.shape != (892, 892) or not np.isfinite(K).all():
        raise ValueError("GRM must be a finite, 892 x 892 numerical matrix.")
    if not np.allclose(K, K.T):
        raise ValueError("Original GRM is not numerically symmetric.")
    if raw["IID"].duplicated().any():
        raise ValueError("Duplicate frog IDs in the matched phenotype table.")

    # Raw source columns are deliberately retained for reproducible derivation.
    required = ["IID", "group", "experiment", "geno_clutch", "geno_zoo",
                "strain", "dose", "euthanised", "days_survived"] + WEEKS
    missing = set(required) - set(raw.columns)
    if missing:
        raise ValueError(f"Columns missing from matched_phenotypes.csv: {sorted(missing)}")
    if not set(BASE_COLUMNS).issubset(base.columns):
        raise ValueError(f"Missing GWAS input columns: {sorted(set(BASE_COLUMNS)-set(base.columns))}")

    raw = raw[required].rename(columns={
        "experiment": "source_experiment",
        "geno_zoo": "source_zoo",
        "geno_clutch": "source_clutch",
    })
    data = base[BASE_COLUMNS].merge(raw, on="IID", how="left", validate="one_to_one", indicator=True)
    if not (data["_merge"] == "both").all():
        raise ValueError("One or more of the original 892 genotyped frogs lack phenotype records.")
    if not data["group"].str.strip().str.lower().eq("treatment").all():
        raise ValueError("An unexposed control unexpectedly entered the analysis.")
    # Compare experiment numbers, allowing differences such as 1 vs 1.0.
    experiment = pd.to_numeric(data["experiment"], errors="raise")
    source_experiment = pd.to_numeric(data["source_experiment"], errors="raise")

    if not experiment.eq(source_experiment).all():
        raise ValueError("Experiment identifiers disagree between sources.")
    if not (data["FID"].astype(str) == data["source_clutch"].astype(str)).all():
        raise ValueError("Clutch / PLINK family IDs disagree between sources.")
    if not (data["geno_zoo"] == data["source_zoo"]).all():
        raise ValueError("Zoo assignments disagree between sources.")
    if not data["strain"].eq("Wastepoint").all():
        raise ValueError("The 86 ITS copies/zoospore conversion is isolate-specific.")
    if not pd.to_numeric(data["dose"], errors="raise").eq(2000).all():
        raise ValueError("Unexpected Bd exposure dose.")

    # Check against the genotype samples actually used in the earlier GWAS.
    fam = pd.read_csv(
        GWAS_ROOT / "infection_qc.fam", sep=r"\s+", header=None,
        dtype={0: str, 1: str},
    )
    if set(zip(fam[0], fam[1])) != set(zip(data["FID"], data["IID"])):
        raise ValueError("Existing infection_qc PLINK samples differ from the 892-frog cohort.")
    for pc in ("PC1", "PC2", "PC3"):
        data[pc] = pd.to_numeric(data[pc], errors="raise")
    if not np.isfinite(data[["PC1", "PC2", "PC3"]].to_numpy()).all():
        raise ValueError("Some required PCs are missing or nonfinite.")

    return data.drop(columns=["_merge"]), K


def make_phenotypes(data):
    """Construct binary survival, early detection, and conditional early burden."""
    data = data.copy()
    event = pd.to_numeric(data["euthanised"], errors="raise")
    days = pd.to_numeric(data["days_survived"], errors="raise")
    if not event.isin([0, 1]).all() or not days.between(1, 98).all():
        raise ValueError("Unexpected endpoint indicator or survival duration.")
    if not days.loc[event == 0].eq(98).all():
        raise ValueError("Non-event frogs must be observed to the day-98 endpoint.")
    data["survived_98"] = (event == 0).astype(int)
    data["euthanasia_event"] = event.astype(int)
    data["survival_days"] = days

    early = data[WEEKS].apply(pd.to_numeric, errors="raise")
    if (early < 0).to_numpy().any():
        raise ValueError("Negative Bd ITS measurements require review.")
    data["complete_weeks_1_6"] = early.notna().all(axis=1)
    # Complete-case restriction guarantees identical six-week observation windows.
    complete = data.loc[data["complete_weeks_1_6"]].copy()
    early_complete = early.loc[complete.index]
    complete["early_detected"] = (early_complete > 0).any(axis=1).astype(int)
    complete["early_max_ITS"] = early_complete.max(axis=1)
    complete["early_max_zoospore"] = complete["early_max_ITS"] / ITS_COPIES_PER_ZOOSPORE
    complete["early_log10_max"] = np.log10(complete["early_max_zoospore"] + 1)
    positive = complete.loc[complete["early_detected"] == 1].copy()

    # Expected counts were independently checked against the uploaded workbook.
    expected = {
        "survival98": (data, 892, "survived_98", 359),
        "early_detection": (complete, 880, "early_detected", 566),
        "early_positive_burden": (positive, 566, None, None),
    }
    for name, (frame, n, col, n_cases) in expected.items():
        if len(frame) != n:
            raise ValueError(f"{name}: expected {n} frogs, found {len(frame)}")
        if col is not None and int(frame[col].sum()) != n_cases:
            raise ValueError(f"{name}: expected {n_cases} positive outcomes")
    if not (positive["early_log10_max"] > 0).all():
        raise ValueError("Positive-burden subset contains nonpositive loads.")
    return {"survival98": data, "early_detection": complete,
            "early_positive_burden": positive}


def export_task(name, frame, original_grm, base_iids, trait, family, skip_plink):
    """Save ID-ordered phenotype, GRM and filtered PLINK genotypes for one GWAS."""
    dest = OUTPUT_ROOT / name
    dest.mkdir(parents=True, exist_ok=True)
    if frame["IID"].duplicated().any() or frame[trait].isna().any():
        raise ValueError(f"{name}: duplicate IDs or missing trait values.")

    columns = BASE_COLUMNS + [trait, "euthanasia_event", "survival_days"]
    if name != "survival98":
        columns += ["early_detected", "early_max_zoospore", "early_log10_max"]
    # Avoid duplicate column labels (e.g. when trait itself is early_detected).
    columns = list(dict.fromkeys(columns))
    frame[columns].to_csv(dest / "phenotype.tsv", sep="\t", index=False)
    frame[["FID", "IID"]].to_csv(dest / "keep.txt", sep="\t", index=False, header=False)
    frame[["FID", "IID"]].to_csv(dest / "grm_ids.tsv", sep="\t", index=False)

    # Select precisely these individuals, preserving phenotype row order.
    pos = pd.Index(base_iids).get_indexer(frame["IID"])
    if (pos < 0).any():
        raise ValueError(f"{name}: sample not found in original GRM.")
    subK = original_grm[np.ix_(pos, pos)]
    np.savetxt(dest / "grm.csv", subK, delimiter=",")
    if subK.shape != (len(frame), len(frame)) or not np.allclose(subK, subK.T):
        raise ValueError(f"{name}: invalid subset GRM.")

    if not skip_plink:
        # Starting from previously QC'd 39,369-SNP data. Recalculate sample-
        # specific missingness and MAC, especially in the 566-frog subset.
        source = GWAS_ROOT / "infection_qc"
        target = dest / "snps_qc"
        cmd = [PLINK2, "--bfile", str(source), "--keep", str(dest / "keep.txt"),
               "--geno", "0.05", "--mac", "20", "--make-bed", "--out", str(target)]
        print("Running:", " ".join(cmd), flush=True)
        subprocess.run(cmd, check=True)
        fam = pd.read_csv(str(target) + ".fam", sep=r"\s+", header=None,
                          dtype={0: str, 1: str})
        bim = pd.read_csv(str(target) + ".bim", sep=r"\s+", header=None)
        if set(zip(fam[0], fam[1])) != set(zip(frame["FID"], frame["IID"])):
            raise ValueError(f"{name}: PLINK --keep returned unexpected samples.")
        if len(bim) == 0:
            raise ValueError(f"{name}: no SNPs survived QC.")
        print(f"{name}: {len(frame)} frogs; {len(bim)} QC-retained SNPs")
    else:
        print(f"{name}: {len(frame)} frogs (PLINK step skipped)")

    return {"analysis": name, "trait": trait, "family": family,
            "n_frogs": len(frame), "n_families": frame["FID"].nunique()}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--skip-plink", action="store_true",
                        help="Only check and export phenotypes/GRMs; no SNP QC")
    args = parser.parse_args()
    OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)
    data, K = load_and_verify_inputs()
    tasks = make_phenotypes(data)
    definitions = {
        "survival98": ("survived_98", "binomial"),
        "early_detection": ("early_detected", "binomial"),
        "early_positive_burden": ("early_log10_max", "gaussian"),
    }

    # Audit early missingness: partial six-week records must not silently be zeroed.
    audit = data[["IID", "experiment", "euthanised", "days_survived"]].copy()
    early = data[WEEKS].apply(pd.to_numeric, errors="raise")
    audit["n_early_measured"] = early.notna().sum(axis=1)
    audit["complete_weeks_1_6"] = early.notna().all(axis=1)
    audit.to_csv(OUTPUT_ROOT / "early_followup_completeness_audit.csv", index=False)

    # Summary of phenotypes: useful for methods and reporting.
    summary_rows = []
    for name, frame in tasks.items():
        trait, family = definitions[name]
        summary_rows.append(export_task(
            name, frame, K, data["IID"], trait, family, args.skip_plink
        ))
        by_exp = frame.groupby("experiment")[trait].agg(
            n="size", mean="mean", median="median"
        )
        by_exp.to_csv(OUTPUT_ROOT / name / "phenotype_by_experiment.csv")
    summary = pd.DataFrame(summary_rows)
    summary.to_csv(OUTPUT_ROOT / "analysis_manifest.tsv", sep="\t", index=False)
    print("\nPREPARED ANALYSES\n", summary.to_string(index=False))
    print("Results directory:", OUTPUT_ROOT)


if __name__ == "__main__":
    main()
