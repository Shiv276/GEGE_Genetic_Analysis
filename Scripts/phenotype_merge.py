import pandas as pd
from config import (BD_PHENOTYPE_FILE, METADATA_FILE, PHENOTYPE_DIR)

"""
The phenotype experiments tested 972 frogs. (exp1: 309, exp2: 333, exp3: 330)

So there are 972 rows in total in the phenotype data.

Of those, 918 were exposed to bd and 54 were unexposed controls.



Phenotype data is available for 910 frogs that were selected post-QC by the authors.

Of those, 6 did not match records in the genotype/metadata set. 2 were not found at all. 4 were putatively found under different names indicating a re-genotype run.
Of the 4 possible matches, data from one frog did not match clutch/zoo data across both files.

As a result, all 6 unresolved frogs were disregarded from further analysis as to not confound sample identity results (Although 3 of them may have truly been a match).

This leaves us with 904 frogs; 892 treated; 12 untreated controls.

"""



def load_data():
    """
    Load the three Bd experiments and genotype sample metadata.
    """

    sheets = ["Exp1", "Exp2", "Exp3"] #The supplied phenotype data was separated into three experiments that exposed all 910 frogs.

    #Load each sheet into one concatenated dataframe
    phenotypes = pd.concat([pd.read_excel(BD_PHENOTYPE_FILE, sheet_name=sheet) for sheet in sheets], ignore_index=True)

    genotype_metadata = pd.read_csv(METADATA_FILE, sep=r"\s+")

    return phenotypes, genotype_metadata



def normalise_ids(ids):
    """
    Fixing some formatting errors visually found.

    IDs in the pheno data had leading 0's (i.e. psco089 instead of psco89 in our metadata)
    Also some capitalisation normalisation required.
    """
    
    return ids.astype("string").str.strip().str.lower() .str.replace(r"^psco_0+(\d+)$", r"psco_\1", regex=True)


def match_phenotypes(phenotypes, genotype_metadata):
    """
    Match genotype and phenotype samples using normalised IDs.
    """

    #Store the IDs, clutches and zoos reported for each frog in the phenotyped data
    phenotypes = phenotypes.rename(columns={"ID": "phenotype_ID", "clutch": "pheno_clutch", "zoo": "pheno_zoo"})

    #Store the IDs, clutches and zoos reported for each frog in the genotype / metedata data.
    genotype_metadata = genotype_metadata.rename(columns={"ID": "IID", "clutch": "geno_clutch", "zoo": "geno_zoo"})

    #This is to later check of they match and agree to ensure data is the same between the two.


    #Create a standardised matching key for both sources.
    phenotypes["match_id"] = normalise_ids(phenotypes["phenotype_ID"])
    genotype_metadata["match_id"] = normalise_ids(genotype_metadata["IID"])

    # Ensure each matching key uniquely identifies one record.
    assert phenotypes["match_id"].is_unique
    assert genotype_metadata["match_id"].is_unique

    #Merge the two datasets so that IDs from both can be bought together. Indicator=True means that the data shows where a match was found.
    matched = genotype_metadata.merge(phenotypes, on="match_id", how="left", validate="one_to_one", indicator=True)

    return matched


def main():

    phenotypes, genotype_metadata = load_data()

    merged = match_phenotypes(phenotypes, genotype_metadata)

    #This is where the relevant ID existed in both datasets from the 'Indicator' arg in the merge.
    matched = merged[merged["_merge"] == "both"].copy()

    #This is where the relevant ID only existed in the genotype/metadata set from the 'Indicator' arg in the merge.
    unmatched = merged[merged["_merge"] == "left_only"].copy()


    print(f"Phenotype records: {len(phenotypes)}")
    print(f"Genotyped frogs: {len(genotype_metadata)}")
    print(f"Matched frogs: {len(matched)}")
    print(f"Unmatched genotypes: {len(unmatched)}")

    # Compare biological metadata from the independent sources.
    clutch_agreement = pd.to_numeric(matched["geno_clutch"]) == pd.to_numeric(matched["pheno_clutch"])
    zoo_agreement = matched["geno_zoo"] == matched["pheno_zoo"]


    print(f"Clutch disagreements: {(~clutch_agreement).sum()}")
    print(f"Zoo disagreements: {(~zoo_agreement).sum()}")

    print("\nMatched experimental groups:")
    print(matched["group"].str.lower().value_counts())

    print("\nUnmatched genotype IDs:")
    print(unmatched["IID"].to_string(index=False))


    merged.to_csv(PHENOTYPE_DIR / "genotype_phenotype_match_check.csv", index=False)

    # Retain only confidently matched individuals.
    matched.to_csv(PHENOTYPE_DIR / "matched_phenotypes.csv", index=False)

    # Retain exposed frogs for subsequent disease-resistance analysis.
    exposed = matched[matched["group"].str.lower() == "treatment"].copy()
    exposed.to_csv(PHENOTYPE_DIR / "exposed_phenotypes.csv",index=False)



    print(f"\nMatched frogs retained: {len(matched)}")
    print(f"Bd-exposed frogs retained: {len(exposed)}")

if __name__ == "__main__":
    main()