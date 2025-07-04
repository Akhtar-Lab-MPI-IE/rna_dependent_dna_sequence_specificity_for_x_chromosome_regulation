[![DOI:10.5281/zenodo.15806913](https://img.shields.io/badge/DOI-10.5281/zenodo.5720285-yellow.svg)](https://doi.org/10.5281/zenodo.15806913)
# RNA-dependent DNA sequence specificity for X chromosome regulation. 

Code repository for Kulkarni et al., containing environments, code and information to reproduce the computational analysis.

## Setup

To re-run the analysis, you need to have a working (snakePipes)[https://github.com/maxplanck-ie/snakepipes] version 3.1.0 installation, available in a conda environment (i.e. named 'snakePipes_3.1.0'). Note that this includes having the environments ('createEnvs'), and have the dm6 genome available in there too.

Next, the raw data can be accessed from GEO:

 - Cut&Run: GSE301023  
 - FLASH: GSE301024  
 - RNA-seq: GSE301025  

Note that the processed files (peaks, bws for Cut&Run or count matrices for FLASH and RNA-seq are available in GEO too). The raw fastq-files should be organized in your working folder (WDIR) as such, and should retain their names as they are deposited in the GEO database.

```
WDIR/
  |---fqs/
       |---CNR/  (8 files, 4 samples)
       |---FLASH/ (24 files, 12 samples)
       |---RNA/ (36 files, 18 samples)
```

Additionally, public data is used as well. This is sourced from GEO accession GSE109901, and includes following files:

 - GSM2973501 	Dmel_female_MLE
 - GSM2973502 	Dmel_female_MSL2
 - GSM2973503 	Dmel_female_WOR
 - GSM2973504 	Dmel_male_MLE
 - GSM2973505 	Dmel_male_MSL2
 - GSM2973506 	Dmel_male_WOR

 Note that these contain two or three replicates per sample. They should be named similarly to the FLASH files included in this study, i.e.:

  - Dme_MSL2_female_rep*  
  - Dme_MSL2_male_rep*  
  - Dme_MLE_female_rep*  
  - DME_MLE_male_rep*  
  - Dme_WOR_female (no replicates for WOR)  
  - DME_WOR_male (no replicates for WOR)  

This brings the total number of files in the FLASH folder to 48, or 24 samples.

With this in place, the environment can be created. Relative from this repository, the conda environment can be created:

 > conda env create -f conf/env.yml -n xreg  
 > conda activate xreg  

Next, the config file (`conf/smk_config.yml`) needs to be updated. Note here that snakePipesEnv should correspond to the actual name of the conda environment that houses the snakePipes environment. Some files (fna, chromsizes, starindex, gtf) are provided through zenodo (10.5281/zenodo.15806913), the others are included in this repository.

  > snakePipesEnv: 'snakePipes_3.1.0'  
  > fna: 'path/to/genome.fa'  
  > chromsizes: 'path/to/genome.chrom.sizes'  
  > starindex: 'path/to/starix'  
  > gtf: 'path/to/gtf'  
  > gff: 'repodir/conf/Drosophila_melanogaster.BDGP6.91.canonical.gff3.gz'  
  > rar: 'repodir/conf/dm6_rar.bed'  
  > canchr: 'repodir/conf/can_chrs.txt'  
  > geneinfo: 'repodir/conf/genes.info'  
  > motifs: 'repodir/conf/motifs.meme'  
  > padj_cutoff: 0.01  

With this in place, the workflow is almost ready to run. Some system specific settings needs to be filled in the `run.py` file:

  > SMK_PROFILE: path to or name of your (snakemake profile)[https://snakemake.readthedocs.io/en/stable/executing/cli.html#profiles]  
  > WDIR: path to your working directory (note that fqs folder with content is located here)  

The analysis can subsequently be run with:

 > python run.py  
