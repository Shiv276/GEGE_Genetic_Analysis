import subprocess
from pathlib import Path
from config import RAW_DATA, KINSHIP_DIR



PLINK2 = "plink2"

#These are ALL individuals in the study (~900)
FULL_PREFIX = RAW_DATA / "plink-infectionPsco"

#These contain one representive individual from each family / subpop. It's useful for LD
#pruning so that closely related family members don't confound linkage.
REPRESENTATIVE_PREFIX = RAW_DATA / "54_psco"



def run_plink(command):
    """
    Runs a PLINK command and shows all its useful info in terminal.
    """

    print("\nRunning:")

    #Joins all the arguments together (as strs) so Python calls the PLINK command correctly
    command = [str(x) for x in command]

    print(" ".join(command))
    print()

    subprocess.run(command, check=True)



#KING Kinship testing BEFORE Pruning
def run_all_snp_kinship():

    output_prefix = KINSHIP_DIR / "kinship_all"

    run_plink([
        PLINK2,
        "--bfile", FULL_PREFIX,
        "--make-king-table",
        "--out", output_prefix
    ])




#LD Pruning
def run_ld_pruning(window_kb=5000, r2_threshold=0.2):
    """
    Prnues for LD sing 5MB as described by the original author/s. r^2 threshold is kept 
    as an arg so that it KING output can be tested across different thresholds to assess
    its effect.
    """

    #The common name of the output files
    output_prefix = (KINSHIP_DIR /f"prune_{window_kb}kb_r2_{r2_threshold}")

    run_plink([
        PLINK2,
        "--bfile",
        REPRESENTATIVE_PREFIX, #Uses files containing one representative indiv from each family
        "--indep-pairwise", #PLINKs flag for LD pruning
        f"{window_kb}kb", #Likely to be 5000kb as per the authors analysis, however it could vary.
        '1', #Move the pruning window forward by 1 SNP each step
        str(r2_threshold), #r2 vals above this number will be excluded
        "--out", #Specifying output names
        output_prefix
        ])

    #Specifically had to change from 'output_prefix.with_suffix(".prune.in")' because Python interpreted the
    #'.' in 'prune_5000kb_r2_0.2.prune.in' as a filetype suffix
    return KINSHIP_DIR / (f"prune_{window_kb}kb_r2_{r2_threshold}.prune.in"), KINSHIP_DIR / (f"prune_{window_kb}kb_r2_{r2_threshold}.prune.out")


#KING Kinship testing post-LD pruning
def run_pruned_kinship(prune_file, label="5000kb_r2_0.2"):

    output_prefix = (KINSHIP_DIR /f"kinship_pruned_{label}")

    run_plink([
        PLINK2,
        "--bfile",
        FULL_PREFIX,
        "--extract",
        prune_file, #--extract on the prune_file separates individuals that did/didn't meet the threshold
        "--make-king-table",
        "--out",
        output_prefix
        ])



#Running the pipeline
def main():

    #Kinship using all available markers (QC-passed SNPs from the original authors)
    run_all_snp_kinship()

    #Determine an approximately LD-independent SNP set using one representative animal per pedigree family.


    #LD pruning analysis:

    #Testing a range of r2 thresholds to gauge stringency consequences. Might also adjust the sliding window too - need more reading
    thresholds = [0.1, 0.2, 0.5, 0.8]
    store_changes = []
    for r2 in thresholds:
        prune_in_file, prune_out_file = run_ld_pruning( window_kb=5000, r2_threshold=r2)

        with open(prune_in_file) as f:
            retained = sum(1 for _ in f)

        with open(prune_out_file) as f:
            removed = sum(1 for _ in f)

        store_changes.append(f"r2={r2}: {retained} retained; {removed} removed")
    
    
        # Calculate KING using this pruned SNP set
        run_pruned_kinship(prune_in_file, label=f"5000kb_r2_{r2}")

    [print(i) for i in store_changes]



if __name__ == "__main__":
    main()