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

    #ID and GRM matrix are detached initially from PLINK output 
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

def load_pca_variance(grm):
    """
    Calculate how much of the genotype variation is explained by each of the 30 principal components.

    Basically, PCA identified axes along which the frogs vary genetically.
    Each PC has an associated eigenvalue which indicates how much variation exists along that axis.

    By comparing the eigenvalue of each PC relative to the sum of all eigenvalues, we can calculate
    how much genotype variation is explained by each PC.
    """

    #Theres no header so I assign the it a name of 'eigenvalue'
    eigenvalues = pd.read_csv(PCA_EIGENVAL_FILE, header=None, names=["eigenvalue"])

    #Creating a 'PC' column - basically 1-30 for each PC
    eigenvalues["PC"] = (range(1, len(eigenvalues) + 1))


    total_variance = np.trace(grm.to_numpy()) #Gets the trace of the GRM (sum of diagonals)

    #Tells us the proportion of variance that each PC contributes to
    eigenvalues['variance_explained_pct'] = eigenvalues["eigenvalue"] / total_variance * 100

    eigenvalues["cumulative_variance_pct"] = eigenvalues["variance_explained_pct"].cumsum()

    
    #Proportion of variation captured by each of the PCs (30).
    #Basically how much genetic variation it captures
    eigenvalues["relative_variance_pc30"] = (eigenvalues["eigenvalue"] / eigenvalues["eigenvalue"].sum())

    #Same effect-size but as a percentage. Helps figure out how much variation is explained by each PC.
    eigenvalues["relative_variance_pct_pc30"] = (eigenvalues["relative_variance_pc30"] * 100)

    return eigenvalues





def save_results(pca, pca_eigenvalues):
    """
    Save analysis-ready PCA outputs for plotting.

    pca: One row per frog containing PC1-PC30 plus clutch, zoo and F_ROH metadata.
    variance: Eigenvalues and relative contribution of each PC among the first30 PCs.
    """
    pca.to_csv(GWAS_STRUCTURE_DIR / "pca_with_metadata.csv", index=False)
    pca_eigenvalues.to_csv(GWAS_STRUCTURE_DIR / "pca_eigenvalues.csv", index=False)


def main():

    pca = load_pca()
    grm, ids = load_grm()
    pca_eigenvalues = load_pca_variance(grm)

    print("PCA relative variance:")
    print(pca_eigenvalues[["PC", "eigenvalue", "relative_variance_pct_pc30"]])

    print('-' * 50)


    """
    Since the GRM and PCA were made sing the same pruned dataset in PLINK2, we know that the trace of the GRM (Tr(K)) or sum of diagonals
    is equal to the sum of all eigenvalues in PCA.

    This allows us to calculate how much of the variance each of our PCs contribute.
    """

    total_variance = np.trace(grm.to_numpy())
    top30_variance = pca_eigenvalues["eigenvalue"].sum()

    print(f"GRM trace: {total_variance:.4f}")
    print(f"First 30 eigenvalues: {top30_variance:.4f}")
    print(f"Variance captured: {top30_variance / total_variance * 100:.2f}%")


    save_results(pca, pca_eigenvalues)



if __name__ == "__main__":
    main()