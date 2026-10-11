
import numpy as np
import pandas as pd

from config import PHENOTYPE_DIR


def get_bd_phenotypes(matched):
    """
    Construct individual-level disease phenotypes for Bd-exposed frogs.

    Returns one row per exposed frog, with infection-status, observed pathogen-burden, survival, and QC variables.
    """

    #Exclude unexposed controls from disease-resistance comparisons.
    exposed = matched.loc[matched["group"].str.strip().str.lower() == "treatment"].copy()


    #Load post-exposure weekly Bd measurements. T here were 14 weeks in total where levels were measured.
    week_cols = [f"Bd_wk{i}" for i in range(1, 15)]

    #Converting each load measurment into a float, incase they get interpreted as a str from Excel.
    loads = exposed[week_cols].astype(float)
    


    #Ensuring the number of positive weeks is > 0. This was True for all frogs.
    exposed["n_positive_weeks"] = loads.gt(0).sum(axis=1)
    exposed["n_measured_weeks"] = loads.count(axis=1)

    #Athors define infection as two positive observations from qPCR. This only differs by one case from the alternative, so it was kept.
    exposed["infected_2plus"] = ( exposed["n_positive_weeks"] >= 2).astype(int)

    #Maximum observed Bd burden, converted from ITS copies to zoospore-equivalent units (86 copies/zoospore). This figure was taken directly from the original authors
    exposed["max_Bd_load"] = loads.max(axis=1) / 86

    #Log10 transformation because Bd load can differ substantially between individuals. (+1 accounts for loads of 0)
    exposed["log10_max_Bd_load"] = np.log10(exposed["max_Bd_load"] + 1)

    #Preserve both parts of the time-to-event outcome.
    exposed["survival_days"] = exposed["days_survived"]
    exposed["euthanasia_event"] = exposed["euthanised"]



    #Only keeping the variables needed for phenotype QC and future GWAS.
    columns = [
        "IID", "phenotype_ID",
        "geno_clutch", "geno_zoo",
        "experiment", "age_start", "bc_start",
        "n_measured_weeks", "n_positive_weeks",
        "infected_2plus", "max_Bd_load", 
        "log10_max_Bd_load", "survival_days",
        "euthanasia_event"
    ]

    return exposed[columns].copy()
    

def main():

    matched = pd.read_csv(PHENOTYPE_DIR / "matched_phenotypes.csv")

    phenotypes = get_bd_phenotypes(matched)

    #Summarise outcomes across experiments.
    summary = phenotypes.groupby("experiment").agg(
        frogs=("IID", "size"),
        positive_2plus=("infected_2plus", "sum"),
        euthanasia_events=("euthanasia_event", "sum")
        )

    print("\nPhenotype summary by experiment:")
    print(summary)

    phenotypes.to_csv(PHENOTYPE_DIR / "bd_phenotypes_preliminary.csv", index=False)


if __name__ == "__main__":
    main()

