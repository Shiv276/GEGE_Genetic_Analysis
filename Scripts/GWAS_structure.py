import numpy as np
import pandas as pd

from config import GWAS_STRUCTURE_DIR, PROCESSED_DATA


PCA_FILE = GWAS_STRUCTURE_DIR / "population_structure_pca.eigenvec"

#Contains the 'sqare' GRM data
GRM_FILE = GWAS_STRUCTURE_DIR / "genomic_relationship_matrix.rel"

#Contains the ID labels for the GRM data (ordered)
GRM_ID_FILE = GWAS_STRUCTURE_DIR / "genomic_relationship_matrix.rel.id"

FROG_FILE = PROCESSED_DATA / "individual_froh_metadata.csv"



def load_pca():

    pca = pd.read_csv(PCA_FILE, sep=r"\s+")
    pca = pca.rename(columns={"#FID": "FID"})

    frogs = pd.read_csv(FROG_FILE)

    #Attach useful biological information to each frog's PCs
    pca = pca.merge(frogs[["IID", "clutch", "zoo", "F_ROH"]], on="IID", how="left", validate="one_to_one")

    return pca



def load_grm():

    ids = pd.read_csv(GRM_ID_FILE, sep=r"\s+").rename(columns={"#FID": "FID"})

    #GRM
    grm = pd.read_csv(GRM_FILE, sep=r"\s+", header=None)


    # Give rows/columns actual frog IDs

    #Sets matrix rows to IID
    grm.index = ids["IID"]

    #Sets columns to same IID
    grm.columns = ids["IID"]

    #Quick way to check diagonal (each indiv compared to iteslf)
    matrix = grm.to_numpy()
    diagonal = np.diag(matrix)
    
    #Should be very close to 1
    print(f"Mean Diagonal in GRM: {diagonal.mean():.2f}")

    return grm, ids



#After loading those, I'll need to analyse the PCA to see how it can be applied to the eventual GWAS
PCA_EIGENVAL_FILE = (GWAS_STRUCTURE_DIR / "population_structure_pca.eigenval")

def load_pca_variance():
    """
    Calculate how much of the genotype variation is explained by each of the 10 principal components.

    Basically, PCA identified axes along which the frogs vary genetically.
    Each PC has an associated eigenvalue which indicates how much variation exists along that axis.

    By comparing the eigenvalue of each PC relative to the sum of all eigenvalues, we can calculate
    how much genotype variation is explained by each PC.
    """

    #Theres no header so I assign the it a name of 'eigenvalue'
    eigenvalues = pd.read_csv(PCA_EIGENVAL_FILE, header=None, names=["eigenvalue"])

    #Creating a 'PC' column - basically 1-10 for each PC
    eigenvalues["PC"] = (range(1, len(eigenvalues) + 1))

    #Proportion of variation captured by each of the PCs (10).
    #Basically how much genetic variation it captures
    eigenvalues["relative_variance_pc10"] = (eigenvalues["eigenvalue"] / eigenvalues["eigenvalue"].sum())

    #Same effect-size but as a percentage. Helps figure out how much variation is explained by each PC.
    eigenvalues["relative_variance_pct_pc10"] = (eigenvalues["relative_variance_pc10"] * 100)

    return eigenvalues


def save_results(pca, variance):
    """
    Save analysis-ready PCA outputs for plotting.

    pca: One row per frog containing PC1-PC10 plus clutch, zoo and F_ROH metadata.
    variance: Eigenvalues and relative contribution of each PC among the first 10 PCs.
    """
    pca.to_csv(GWAS_STRUCTURE_DIR / "pca_with_metadata.csv", index=False)
    variance.to_csv(GWAS_STRUCTURE_DIR / "pca_eigenvalues.csv", index=False)


def main():

    pca = load_pca()
    grm, ids = load_grm()
    variance = load_pca_variance()

    print("PCA relative variance:")
    print(variance[["PC", "eigenvalue", "relative_variance_pct_pc10"]])

    save_results(pca,variance)



if __name__ == "__main__":
    main()