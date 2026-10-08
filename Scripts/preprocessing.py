import pandas as pd
import matplotlib.pyplot as plt

from config import (
    NE_FILE,
    ROH_FILE,
    METADATA_FILE,
    PROCESSED_DATA,
    GENOME_SIZE_KB
)


#Load Datasets
def load_ne_data():
    """
    Load historical effective population size estimates.
    """
    return pd.read_csv(NE_FILE, sep=r"\s+", header=None, names=["Generations_ago", "Ne"])


def load_roh_data():
    """
    Load PLINK ROH output.
    """
    return pd.read_csv(ROH_FILE, sep=r"\s+")


def load_metadata():
    """
    Load frog metadata containing ID, genotype order, pedigree family (clutch), and zoo.
    """
    return pd.read_csv(METADATA_FILE, sep=r"\s+")




#Individual Genomic Metrics
def calculate_froh(df_roh):
    """
    Calculate  genomic inbreeding metrics (ROH) for each individual frog.
    Takes a raw dataframe containg ROH data as input. This would be one row
    per ROH observation.
    """

    #Assigning each ROH to its respective frog + summary stats.
    froh = df_roh.groupby("IID", as_index=False).agg(
        total_ROH_kb=("KB", "sum"),
        number_ROH=("KB", "size"),
        mean_ROH_kb=("KB", "mean"),
    )
    

    # Proportion of analysed genome contained in ROH
    froh["F_ROH"] = (froh["total_ROH_kb"] / GENOME_SIZE_KB)

    # Using percentages for more interpretable plots.
    froh["F_ROH_pct"] = froh["F_ROH"] * 100

    return froh


#Attach pedigree / family data
def attach_metadata(froh, metadata):
    """
    Match each individual's genomic metrics with its
    pedigree family and zoo information (they call it clutch).

    This allows us to later look at ROH within-family / subpopulation
    """

    #Merging ROH data to frog metadata.
    froh_meta = froh.merge(
        metadata,
        left_on="IID",
        right_on="ID",
        how="left",
        validate="one_to_one"
        #^Should ideally have one ROH-grouped individual per metadata individual
    )


    #Every genotyped frog should ideally map to metadata.
    unmatched = froh_meta["clutch"].isna().sum()
    if unmatched != 0:
        raise ValueError(f"{unmatched} frogs could not be matched to metadata.")


    return froh_meta



#Create a family-level summary for descriptive plots
def summarise_family_froh(froh_meta):
    """
    Summarises genomic inbreeding across pedigree families.
    """

    #Generate summary ROH stats for each family.
    family_froh = froh_meta.groupby("clutch", as_index=False).agg(
        mean_FROH=("F_ROH", "mean"),
        sd_FROH=("F_ROH", "std"),
        n=("F_ROH", "size")
        )
    

    # Standard error around family mean (SE = SD / sqrt(n))
    family_froh["se_FROH"] = (family_froh["sd_FROH"] / family_froh["n"] ** 0.5)

    # Percentage transform for neater plots.
    family_froh["mean_FROH_pct"] = (family_froh["mean_FROH"] * 100)
    family_froh["se_FROH_pct"] = (family_froh["se_FROH"] * 100)

    #Sorting so that the plots are shown in order of increasing ROH across family.
    return family_froh.sort_values("mean_FROH")



#Save data for future scripts

def save_processed_data( df_ne, froh, froh_meta, family_froh):
    """
    Save clean intermediate datasets so later scripts
    do not need to repeat preprocessing.
    """

    df_ne.to_csv(PROCESSED_DATA / "historical_ne.csv",index=False)
    froh.to_csv(PROCESSED_DATA / "individual_froh.csv", index=False)
    froh_meta.to_csv(PROCESSED_DATA / "individual_froh_metadata.csv", index=False)
    family_froh.to_csv(PROCESSED_DATA / "family_froh.csv", index=False)


#Run pipeline
def main():

    #Load
    df_ne = load_ne_data()
    df_roh = load_roh_data()
    metadata = load_metadata()

    #Calculate
    froh = calculate_froh(df_roh)

    #Add family/sample information
    froh_meta = attach_metadata(froh,metadata)

    #Family summaries
    family_froh = summarise_family_froh(froh_meta)

    #Save
    save_processed_data(df_ne, froh, froh_meta, family_froh)


    #Basic QC / useful metrics (docstringed out but may be useful to check for others)
    
    """
    print(f"Individuals analysed: {len(froh)}")

    print(f"Pedigree families: {froh_meta['clutch'].nunique()}")

    print("\nF_ROH summary:")
    print(froh["F_ROH"].describe())
    """



if __name__ == "__main__":
    main()