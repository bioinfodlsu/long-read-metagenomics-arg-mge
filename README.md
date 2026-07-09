# Investigating the landscape of antibiotic resistance genes and their mobile genetic elements carriers across diverse taxa and environments using long-read metagenomics sequence data

This serves as a guide to run the analysis pipeline written in Snakemake.

## Installation for Snakemake
This Snakemake pipeline requires the package manager **Conda** and the workflow management system **Snakemake**.
Additional dependencies not handled by Snakemake are described in Section 1.3.

### 1.1. Install Miniconda 
```
$ curl -sL \
  "https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh" > \
  "Miniconda3.sh"
$ bash Miniconda3.sh
$ conda update conda
$ rm Miniconda3.sh
$ conda install wget
```

### 1.2. Install Mamba 
```
$ conda config --add channels conda-forge
$ conda update -n base --all
$ conda install -n base mamba
```

### 1.3. Install Snakemake
```
$ mamba create -c conda-forge -c bioconda -n snakemake snakemake
```
This creates an isolated enviroment containing the latest Snakemake. To activate it:
```
$ conda activate snakemake
```
To test snakemake:
```
$ snakemake --help
```

### 1.4. Install Additional Dependencies
Install git and gawk. We require gawk to process the filtering stage of our databases.
```
$ mamba install git
$ mamba install gawk
```

### 1.5. Download the pipeline
Download the online [repository](https://github.com/bioinfodlsu/long-read-metagenomics-arg-mge), or using the command line:
```
$ git clone https://github.com/bioinfodlsu/long-read-metagenomics-arg-mge
```

## 2. Quickstart Usage Guide

### 2.1. Running the pipeline
With the snakemake conda environment activated, you can call the pipeline from the top-level directory:
```
$ cd long-read-metagenomics-arg-mge
$ snakemake --use-conda --cores all
```
In case of errors encountered relating to the use of conda environments, please use the following command:
```
$ snakemake --use-conda --cores all --conda-frontend conda

```
