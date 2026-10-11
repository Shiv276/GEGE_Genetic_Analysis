from config import PLINK_PREFIX, PRIMARY_PRUNE_FILE, GWAS_STRUCTURE_DIR

#Borrowing functionality from earlier command-line adjacent work (kinship / LD processing)
from LD_pipeline import run_plink, PLINK2


#PCA
def run_structure_pca(n_pcs=30):
    """
    Calculate principal components from the LD-reduced SNP set.

    These PCs describe major genome-wide axes of genetic variation among the 910 frogs.
    """

    output_prefix = GWAS_STRUCTURE_DIR / "population_structure_pca"

    run_plink([
        PLINK2,
        "--bfile",
        PLINK_PREFIX, 
        "--extract",
        PRIMARY_PRUNE_FILE, #Extracting the .in file post-pruning
        "--pca", #Generating a PCA
        str(n_pcs), #Default var set to 10 - may adjust based on literature and genomic structure. UPDATE: Set to 30 to replicate original workflow from the authors
        "--out",
        output_prefix
    ])


#GRM
def run_grm():
    """
    Calculate the genomic relationship matrix (GRM).

    Creates a square matrix describing genome-wide covariance among all frogs.
    """

    output_prefix = GWAS_STRUCTURE_DIR /"genomic_relationship_matrix"

    run_plink([
        PLINK2,
        "--bfile",
        PLINK_PREFIX,
        "--extract",
        PRIMARY_PRUNE_FILE,
        "--make-rel", #Flag for a GRM
        "square", #Makes it fully symmetric and easier to work with (default is a triangular structure)
        "--out",
        output_prefix
    ])




def main():

    run_structure_pca(n_pcs=30)
    run_grm()


if __name__ == "__main__":
    main()