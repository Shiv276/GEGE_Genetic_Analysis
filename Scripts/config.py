from pathlib import Path

"""
This primarily functions to assign directories for each major contribution,
and keep data organised.
"""

#Current dir (scripts) and then back one
PROJECT_DIR = Path(__file__).resolve().parents[1]

#Where / how I had the raw data saved
DRYAD_DIR = Path("../Psco-SNParray_Data/Psco-SNParray_Data copy").resolve()
RAW_DATA = DRYAD_DIR / "data"


#Setting up dirs for my further analysis
PROCESSED_DATA = PROJECT_DIR / "data" / "processed"
FIGURES = PROJECT_DIR / "outputs" / "figures"
TABLES = PROJECT_DIR / "outputs" / "tables"
KINSHIP_DIR = PROCESSED_DATA / "kinship"
GWAS_STRUCTURE_DIR = PROCESSED_DATA / "gwas_structure"


GWAS_STRUCTURE_DIR.mkdir(parents=True, exist_ok=True)
PROCESSED_DATA.mkdir(parents=True, exist_ok=True)
FIGURES.mkdir(parents=True, exist_ok=True)
TABLES.mkdir(parents=True, exist_ok=True)
KINSHIP_DIR.mkdir(parents=True, exist_ok=True)


#All the files I'll be using in my analysis
NE_FILE = RAW_DATA / "ld_ne_raw.txt"
ROH_FILE = RAW_DATA / "roh-allAnimals.hom"
METADATA_FILE = RAW_DATA / "list-infectionPsco910.txt"

#Currently moving forward with r^2 = 0.5 for LD pruning. May adjust later
PRIMARY_PRUNE_FILE = KINSHIP_DIR / "prune_5000kb_r2_0.2.prune.in"



PLINK_PREFIX = RAW_DATA / "plink-infectionPsco"
#*These should be accessed via; .with_suffix(".bed/.bim/.fam")


GENOME_SIZE_KB = 8190274

