
from pathlib import Path
import subprocess
import pandas as pd

from config import PLINK_PREFIX, PHENOTYPE_DIR, GWAS_DIR
from LD_pipeline import run_plink, PLINK2


"""
This prepares the genotype data used for GWAS which will be carried out by GMMAT. Filters SNPs that have low frequencies
"""




def main():
    """
    Prepare QC-filtered genotype data for the Bd infection GWAS.

    Retain the 892 phenotype-matched exposed frogs, filter
    SNP missingness and minor allele counts, and export
    PLINK BED/BIM/FAM files for GMMAT.
    """

    keep_file = GWAS_DIR / "gwas_keep.txt"
    output_prefix = GWAS_DIR / "infection_qc"


    #Run PLINK2.
    run_plink([
        PLINK2,
        "--bfile", str(PLINK_PREFIX),
        "--keep", str(keep_file),
        "--geno", "0.05",
        "--mac", "20",
        "--make-bed",
        "--out", str(output_prefix)
    ])


    # Verify the samples and count the retained variants.
    fam = pd.read_csv(f"{output_prefix}.fam", sep=r"\s+", header=None, dtype=str)
    bim = pd.read_csv(f"{output_prefix}.bim", sep=r"\s+", header=None)
    keep = pd.read_csv(keep_file, sep=r"\s+", header=None, dtype=str)


    print(f"Retained frogs: {len(fam)}")
    print(f"Retained SNPs: {len(bim)}")
    print(f"GWAS folder: {GWAS_DIR.resolve()}")


if __name__ == "__main__":
    main()
