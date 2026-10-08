# Data from: Design and application of a genome-wide SNP array to improve conservation outcomes in the critically endangered Southern Corroboree Frog

[https://doi.org/10.5061/dryad.ngf1vhj4b](https://doi.org/10.5061/dryad.ngf1vhj4b)

## Description of the data and file structure

Principle Investigator Contact Information

```
Name: Mikaeylah Davidson
Institution: James Cook University
Email: mikaeylah.davidson@jcu.edu.au
```

This dataset contains the data and code required to replicate analyses in Davidson, M. J. et al. Design and Application of a Genome-Wide SNP Array to Improve Conservation Outcomes in the Critically Endangered Southern Corroboree Frog. Molecular Ecology Resources 26, e70080 (2026). [https://doi.org/https://doi.org/10.1111/1755-0998.70080](https://doi.org/https://doi.org/10.1111/1755-0998.70080).

### Archive structure

All data, scripts and supplementary files are contained within: **Psco-SNParray_Data.zip**

The archive contains an R project file: **Psco_SNParray.Rproj**

Opening this project file in RStudio will preserve the relative file paths used by the analysis scripts.

The archive contains the following folders and files:

```
Psco-SNParray_Data.zip

├── 01_Array-design_scripts
├── 02_PopGen_scipts
├── data
├── figures
├── Psco_SNParray.R
└── Psco_SNParray.Rproj
```

### Scripts

#### 01_Array-design_scripts

Contains scripts used for the design of the custom SNP array. Open and follow the README_design.sh file within this folder. All files and labeled in numerical order, and their purpose laid out in the README_design.sh These scripts are provided as workflow documentation and do not require modification for replication of the analyses.

Files:

* `README_design.sh`\
  Description and workflow instructions for SNP array design analysis.
* `1_Genome.sh`\
  Genome preparation steps.
* `2_Quality-control.sh`\
  Sequence quality control workflow.
* `3_Align-merge.sh`\
  Sequence alignment and merging steps.
* `4_Variant-calling.sh`\
  Variant calling workflow.
* `5_Hard-filtering.sh`\
  Variant filtering steps.
* `6_SNP-filtering.sh`\
  SNP selection and filtering workflow.
* `7_Submission-formating.sh`\
  Formatting of SNP data for array submission.
* `8_Submitted-SNPs.txt`\
  Final SNP list submitted for array development.
* `A_samtools-index-file.slurm.sh`\
  Samtools indexing script.
* `B_bwa-genome-index.slurm.sh`\
  Genome indexing script.
* `C_fastqc-run.slurm.sh`\
  FastQC quality assessment script.
* `D_multiqc.slurm.sh`\
  MultiQC summary generation script.
* `E_fastp-array.slurm.sh`\
  Sequence trimming and quality filtering script.
* `F_runBWA.slurm.sh`\
  BWA alignment script.
* `G_merge-remDups.slurm.sh`\
  Alignment merging and duplicate removal script.
* `H_variantCalling-bcftools.slurm.sh`\
  Variant calling using bcftools.
* `I_71mer-format.py`\
  Python script used for SNP formatting.

#### 02_PopGen_scripts

Contains scripts used for the population genetic analyses part of the manuscript. Each script has explanations and works stepwise through the analysis. 

Files:

* `1_QC.sh`\
  Quality control workflow for population genomic analyses.
* `2_PopulationDiversity.sh`\
  Population diversity analyses.

#### Psco_SNParray.R

R script containing all statistical analyses and figure generation.

The script includes:

* data loading;
* quality control;
* population genetic analyses;
* figure generation.

### Data files

All data files required for analyses are contained within the **data** folder.

#### Genotype files

#### PLINK binary files

The following files comprise PLINK datasets:

1. Genotype data for all individuals included in the disease challenge experiment.

* `plink-infectionPsco.bed`
* `plink-infectionPsco.bim`
* `plink-infectionPsco.fam`

     2.Data for one representative individual from each pedigree family.

* `54_psco.bed`
* `54_psco.bim`
* `54_psco.fam`

The `.bed`file contains binary genotype data, while the accompanying `.bim` and `.fam` files contain variant and sample metadata required for interpreting the genotype data.

Users can access PLINK datasets using PLINK software ([https://www.cog-genomics.org/plink/). ](https://www.cog-genomics.org/plink/)The three associated PLINK files (`.bed`, `.bim`, and `.fam`) must remain together in the same directory. For example, the command below converts the dataset into VCF format for viewing or downstream analyses:

```
plink --bfile plink-infectionPsco --recode vcf --out infection_data
```

#### Variant Call Format files

These vcf files have already been generated, and can be found in;

* `psco_SNPs_rawdata.vcf`

Raw genotype calls for all individuals prior to downstream population genetic analyses.

* `infection-pscoONLY.recode.vcf`

Filtered genotype dataset containing individuals included in the infection experiment.

### Tabular data files and variable descriptions

#### Overall-data.xlsx

Excel workbook containing sample metadata, experimental information, phenotypic measurements, and SNP array quality control data for individuals included in the study.

The workbook contains three sheets:

1\. `DNA-SNP-info`

This sheet was used for sample tracking and SNP array quality control analyses. Contains sample metadata and SNP array quality control information.

Variables are:

* `ID`: Unique identifier for each individual sample.
* `species`: Species identification.
* `group`: Experimental group assigned to each individual.
* `experiment`: Experiment or dataset associated with the sample.
* `tissue_type`: Tissue type used for DNA extraction.
* `preservative`: Preservation method used for tissue storage.
* `plate_no`: SNP array plate identifier.
* `260_280`: DNA purity ratio based on absorbance at 260 nm and 280 nm.
* `260_230`: DNA purity ratio based on absorbance at 260 nm and 230 nm.
* `Conc_after_dilution`: DNA concentration following dilution.
* `Pass_Fail`: SNP array quality control pass/fail status.
* `DQC`: Dish quality control metric from SNP array processing.
* `QC_call_rate`: Genotype call rate following SNP array quality control.

The second sheet  

  2.`DNA-tissue-swabs`

Contains the sample metadata and quality control information for the five individuals used for comparing DNA quality and SNP array performance across different tissue types, as outlined in the manuscript. 

Variables are:

* `ID`: Unique identifier for each individual sample.
* `species`: Species identification.
* `group`: Experimental group assigned to each individual.
* `experiment`: Experiment or dataset associated with the sample.
* `tissue_type`: Tissue type used for DNA extraction.
* `preservative`: Preservation method used for tissue storage.

- `260_280`: DNA purity ratio based on absorbance at 260 nm and 280 nm.
- `260_230`: DNA purity ratio based on absorbance at 260 nm and 230 nm.
- `Conc_after_dilution`: DNA concentration following dilution.
- `Pass_Fail`: SNP array quality control pass/fail status.
- `DQC`: Dish quality control metric from SNP array processing.
- `QC_call_rate`: Genotype call rate following SNP array quality control.

The third sheet  

  3.`cross-species`

Contains the sample metadata and quality control information for  four species, with six individuals per species used for comparing the performance of the array across closely related species, as outlined in the manuscript. 

Variables are:

* `ID`: Unique identifier for each individual sample.
* `species`: Species identification.
* `DQC`: Dish quality control metric from SNP array processing.
* `QC_call_rate`: Genotype call rate following SNP array quality control.

#### SNeP.xlsx

Excel workbook containing tabulated output generated from SNeP analysis used to estimate historical effective population size (Ne) based on linkage disequilibrium.

The workbook contains two sheets, `F0` and `F1`. Both sheets contain the following variables:

* GenAgo: Estimated number of generations ago represented by the effective population size estimate.
* Ne: Estimated effective population size.
* dist: Distance between SNP markers used for linkage disequilibrium calculations.
* r2: Mean linkage disequilibrium (r²) value calculated between SNP pairs.
* r2SD: Standard deviation of r² estimates.
* items: Number of SNP pairs included in the calculation.

This file is loaded within the provided R script for downstream analyses and figure generation. The workbook does not require modification for replication of the analyses.

#### Population genetic output files

The following files contain output generated from population genetic analyses described in the provided scripts. These files are used for downstream statistical analyses and figure generation in R.

* `ld_adjacent.txt`: Output from linkage disequilibrium analyses between adjacent SNP markers. Used for downstream plotting.
* `LD-thinned5000kb.ld`: Output of LD estimates calculated using SNPs thinned at 5000 kb intervals.
* `ld_ne_raw.txt`: Output from linkage disequilibrium-based effective population size analyses.
* `het910.het`: Output from PLINK heterozygosity calculations.
* `roh-allAnimals.hom`: Output from runs of homozygosity analyses.
* `roh-allAnimals.hom.summary`: Summary statistics from runs of homozygosity analyses.

#### Supporting files

The following supporting files are required metadata used for downstream analysis and plotting in R. 

* `chr-size.txt`: Text file containing chromosome sizes used for genomic analyses.
* `Tiled-SNPs.txt`: List of SNP markers included on the custom SNP array, including genomic position information.
* `list-infection-animals.txt`: List of individuals included in the disease challenge experiment, excluding controls.
* `list-infectionPsco910.txt`: List of individuals included in infection analyses following quality control.
* `test-fams.txt`: Pedigree family information used for analysis validation.

#### Figures

The `figures` folder is empty. It is there as the output file for placing generating figures by running the analysis in `Psco_SNParray.R`

### Software requirements

All analyses were conducted on The University of Melbourne's high performance computing system Spartan. Required software packages and dependencies are specified within the relevant scripts. Users reproducing these analyses on alternative systems may need to modify resource specifications according to their computing environment.

All other statistical analyses and figure generation were conducted using R.

The analysis can be run using RStudio by opening:

* `Psco_SNParray.Rproj`

Which will load the following code file

* `Psco_SNParray.R `

Which contains explanations of each step and works through a stepwise manner.

Required packages are outlined and loaded within the script. Package versions may vary depending on the user's R installation.

##### Additional genomic resources

The *P. corroboree* genome and whole genome sequences are available to download from NCBI, under the BioProject PRJNA1088171. 
