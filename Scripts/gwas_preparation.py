from pathlib import Path

import numpy as np
import pandas as pd

from config import PHENOTYPE_DIR, PLINK_PREFIX, GWAS_STRUCTURE_DIR, GWAS_DIR




def main():
    """
    Make the .fam PLINK files pre-GWAS be ready for GWAS by removing the individuals that have uncharacterised phenotypes (keep 892).
    Does the same with GRM.
    """
    
    #Load the disease phenotypes already constructed.
    pheno = pd.read_csv(PHENOTYPE_DIR / "bd_phenotypes_preliminary.csv")

    #Load sample identifiers from the actual genotype file.
    fam = pd.read_csv(f"{PLINK_PREFIX}.fam", sep=r"\s+", header=None, names=["FID", "IID", "father", "mother", "sex", "phenotype"], dtype={"FID": str, "IID": str})

    #Match the phenotype records to the real PLINK sample identifiers.
    samples = fam[["FID", "IID"]].merge(pheno, on="IID", how="inner", validate="one_to_one")
    

    #Selecting only the columns needed for the initial association model
    gwas_input = samples[["FID", "IID", "infected_2plus", "experiment", "geno_zoo"]].copy()
    gwas_input.to_csv(GWAS_DIR / "gwas_input.tsv", sep="\t", index=False)


    #PLINK --keep needs FID and IID, without a header. These are the 892 files that will be kept in GWAS analysis.
    samples[["FID", "IID"]].to_csv(GWAS_DIR / "gwas_keep.txt", sep="\t", index=False, header=False)


    #Load the existing 910-frog GRM.
    grm_file = GWAS_STRUCTURE_DIR / "genomic_relationship_matrix.rel"
    grm_id_file = GWAS_STRUCTURE_DIR / "genomic_relationship_matrix.rel.id"

    grm = np.loadtxt(grm_file)
    grm_ids = pd.read_csv(grm_id_file, sep=r"\s+")

    #Find each selected frog's row/column position in the original GRM.
    positions = pd.Index(grm_ids.iloc[:, 1].astype(str)).get_indexer(samples["IID"])

    #Subset rows AND columns in exactly the phenotype sample order.
    selected_grm = grm[np.ix_(positions, positions)]

    #Save modified and re-ordered GRM + IDs
    np.savetxt(GWAS_DIR / "grm_892.csv", selected_grm, delimiter=",")
    samples[["FID", "IID"]].to_csv(GWAS_DIR / "grm_892_ids.tsv", sep="\t", index=False)

    
    
    print(f"GWAS samples: {len(samples)}")
    print(f"Infected: {samples['infected_2plus'].sum()}")
    print(f"Not infected: {(samples['infected_2plus'] == 0).sum()}")
    print(f"GRM shape: {selected_grm.shape}")


if __name__ == "__main__":
    main()