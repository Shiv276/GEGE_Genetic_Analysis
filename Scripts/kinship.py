import pandas as pd

from config import KINSHIP_DIR, PROCESSED_DATA

"""
Establishes / measures genetic relatedness across the whole captive population. May be useful in breeding decisions (mostly explorative)
"""

PRIMARY_KINSHIP_FILE = (KINSHIP_DIR /"kinship_pruned_5000kb_r2_0.5.kin0")
FROG_FILE = (PROCESSED_DATA /"individual_froh_metadata.csv")


def load_data():
    kinship = pd.read_csv(PRIMARY_KINSHIP_FILE, sep=r"\s+")
    frogs = pd.read_csv(FROG_FILE)
    kinship = kinship.rename(columns={"IID1": "ID1", "IID2": "ID2"})

    return kinship, frogs


def attach_pair_metadata(kinship, frogs):
    """
    Largely follows the logic of attach_family_metadata() in kinship_stringency_checker.py to
    merge frog metadata with KINSHIP data (with some extra parameters).
    """
    #Keep only information that will be useful for interpreting pairwise relationships.
    metadata = frogs[["IID", "clutch", "zoo", "F_ROH", "F_ROH_pct"]].copy()


    #Metadata for frog 1 in each pair.
    frog1 = metadata.rename(columns={"IID": "ID1", "clutch": "clutch1", "zoo": "zoo1", "F_ROH": "F_ROH1", "F_ROH_pct": "F_ROH_pct1"})


    #Metadata for frog 2 in each pair.
    frog2 = metadata.rename(columns={"IID": "ID2", "clutch": "clutch2", "zoo": "zoo2", "F_ROH": "F_ROH2", "F_ROH_pct": "F_ROH_pct2"})

    #Merge them all together as done in kinship_stringency_checker.py
    pairs = kinship.merge(frog1, on="ID1", how="left", validate="many_to_one")
    pairs = pairs.merge(frog2, on="ID2", how="left", validate="many_to_one")


    #Known pedigree relationships.
    pairs["same_family"] = (pairs["clutch1"] ==pairs["clutch2"])

    #Captive zoo/institution relationship.
    pairs["same_zoo"] = (pairs["zoo1"] == pairs["zoo2"])

    return pairs



def calculate_individual_kinship(pairs):
    """
    Summarises pairwise genetic KINSHIP at the individual frog level.

    Done by first separting ID1 and ID2 (since some frogs only appear in one, or in unequal amounts across both).
    KINSHIP info is retained for each individual.

    pd.concat() on these two separated dfs essentially puts all individuals in one column, with their associative KINSHIPs
    in the other column. 
    
    Technically we've lost information on which KINSHIP refers to which pair - however it's not required
    because the function provides summary statistics on each frogs relatedness to the pooled population.
    """
    #Separate the first and second ID rows so that we can view each individual separately
    #and determine its KINSHIP metrics relative to the whole population
    first = pairs[["ID1", "KINSHIP"]].rename(columns={"ID1": "ID"})
    second = pairs[["ID2", "KINSHIP"]].rename(columns={"ID2": "ID"})


    # Stack both orientations so each frog receives all of its pairwise relationships.
    long_kinship = pd.concat([first, second], ignore_index=True)

    individual = long_kinship.groupby("ID", as_index=False).agg(
        mean_kinship=( "KINSHIP", "mean"),
        median_kinship=("KINSHIP", "median"),
        min_kinship=("KINSHIP", "min"),
        max_kinship=("KINSHIP", "max"),
        n_pairs=("KINSHIP", "size")
    )

    return individual



def save_results(pairs, individual):

    pairs.to_csv(KINSHIP_DIR / "kinship_pairs_primary.csv", index=False)
    individual.to_csv(KINSHIP_DIR / "individual_kinship_primary.csv", index=False)


def main():

    kinship, frogs = load_data()
    pairs = attach_pair_metadata(kinship,frogs)
    individual = calculate_individual_kinship(pairs)
    
    save_results(pairs,individual)


if __name__ == "__main__":
    main()