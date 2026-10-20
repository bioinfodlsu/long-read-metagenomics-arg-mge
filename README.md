# Expanding the Repertoire of Associations Between Antibiotic Resistance Genes and Putative Mobile Genetic Element Carriers Using Metagenomic Long Reads

This document explains how to run the analysis pipeline, which is written in Snakemake. The pipeline finds antibiotic resistance genes (ARGs) in long-read metagenomic reads. It then annotates the reads with nearby mobile genetic elements (MGEs), plasmid context and taxonomy.

The repository has two versions of the pipeline:

- **Main pipeline** ([workflow/Snakefile](workflow/Snakefile), [scripts/](scripts)): the end-to-end workflow. Outputs go to `results/`.
- **Modular pipeline** ([modular/](modular)): ARG, MGE, plasmid, taxonomy analysis step can be ran separately. Outputs go to `results_v1/`. See [Section 4](#4-modular-pipeline).

## Pipeline overview

1. **ARG identification**: align reads to CARD with LAST.
2. **Flanking-region extraction**: cut each ARG-carrying read into its non-ARG fragments (before, between and after the ARG hits). Fragments shorter than 50 bp are dropped.
3. **MGE alignment**: align those fragments to ISFinder, Tn3, Tn7, Tn554 and the integron database.
4. **Plasmid classification**: align the full ARG-carrying reads to the plasmid database.
5. **Taxonomy analysis**: classify the full ARG-carrying reads with Kraken2 (+ Bracken).

## 1. Installation

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
This creates an isolated environment containing the latest Snakemake. To activate it:
```
$ conda activate snakemake
```
To test Snakemake:
```
$ snakemake --help
```

### 1.4. Download the pipeline
Download the online [repository](https://github.com/bioinfodlsu/long-read-metagenomics-arg-mge), or use the command line:
```
$ git clone https://github.com/bioinfodlsu/long-read-metagenomics-arg-mge
```

## 2. Inputs

Both pipelines read the same `samples/`, `database/` and `metadata/` folders.

| Input | Location | Notes |
|---|---|---|
| Sample reads | `samples/{sample}.fastq` or `samples/{sample}.fasta` | One file per sample. Adding a file adds a sample. FASTQ files are converted to FASTA first. |
| Reference databases | `database/` | `card_protein_homolog_model.fasta`, `ISFinder.fasta`, `Tn3.fasta`, `Tn7.fasta`, `Tn554.fasta`, `integron.fasta`, `plasmid.fasta`. The LAST indexes are built into `database/last/db/` automatically. |
| Kraken2 database | `database/kraken2/gtdb_kraken2/` | GTDB Kraken2 database. **The pipeline does not build it**, so download it before running the taxonomy step. |
| CARD metadata | `metadata/aro_index.tsv` | |
| MGE metadata | `metadata/ISFinder_metadata.txt`, `metadata/integrase_headers.txt` | |
| Plasmid metadata | `metadata/typing.csv`, `metadata/plasmidfinder.csv` | |
| Read lengths | `metadata/avg_read_length.csv` | Average read length per sample, used by Bracken. **Add a row for each new sample.** |

Example `metadata/avg_read_length.csv`:
```
Samples,Avg_Read_Length
SRR28740373,3000
SRR28740374,2894
```

### Kraken2 database

Our Kraken2 database was built from the representative genome sequences of GTDB r89. Representative bacterial genomes with < 95% completeness or > 5% contamination were left out. The database was then downsampled with `--max-db-size`, capping its size at 64 GB.

**Taxonomy results depend on the Kraken2 database.** If you use a different GTDB release, a different set of genomes or a different size cap, the taxonomy columns and Bracken abundances may differ from ours.

## 3. Running the main pipeline

With the snakemake conda environment activated, call the pipeline from the top-level directory:
```
$ cd long-read-metagenomics-arg-mge
$ snakemake --use-conda --cores all
```
If you get errors about the conda environments, use this command instead:
```
$ snakemake --use-conda --cores all --conda-frontend conda
```
Its final output is `results/final_summary/{sample}_final_summary.csv`.

| Output | Contents |
|---|---|
| `results/arg_summary/{sample}_CARD_summary.csv` | ARG hits, one row per ARG-carrying read |
| `results/arg_is_merged/{sample}_ARG_IS_summary.csv` | + IS hits (composite transposons are labelled in the analysis notebook after filtering) |
| `results/all_summary/{sample}_all_summary.csv` | + Tn3, Tn7, Tn554 and integron columns |
| `results/final_summary/{sample}_final_summary.csv` | + plasmid and taxonomy columns |

### 3.1. Test run

Before running your own samples, run the small test sample in [test/](test) to check that the pipeline and databases are set up correctly. The folder contains:

- `test/test.fasta`: 10 ARG-carrying long reads from SRR32852254.
- `test/test_final_summary.csv`: tthe expected final summary for these reads. Your output should have the same reads and columns as this file, but the exact values may differ. The taxonomy columns depend on the Kraken2 database you use (see [Kraken2 database](#kraken2-database)).

The test needs the same databases as a full run, including the Kraken2 database (Section 2). `metadata/avg_read_length.csv` already has a row for the sample `test`.

1. Copy the test reads into `samples/`. The file name becomes the sample name, `test`:
   ```
   $ mkdir -p samples
   $ cp test/test.fasta samples/test.fasta
   ```
2. Run the pipeline for the test sample only:
   ```
   $ snakemake --use-conda --cores all results/final_summary/test_final_summary.csv
   ```
   Naming the output file runs only the rules that the test sample needs. It skips Bracken, so it doesn't build the Bracken k-mer distributions, which are slow to build. The first run takes longer than later ones, because it also creates the conda environments and builds the LAST indexes. Once the conda environments are created, the test run takes about 5 minutes on a machine with 20 cores and 64 GB of RAM.
3. Remove the test sample before running your own data, so that it isn't processed with them:
   ```
   $ rm samples/test.fasta
   ```

## 4. Modular pipeline

The modular pipeline in [modular/](modular) splits the analysis into independent modules. Each module is a script that takes its own inputs and writes one CSV keyed on `query_name`. A merge script joins them. The following can be done:

- run the whole pipeline, ARG → MGE → plasmid → taxonomy, in one command;
- run a single analysis (for example MGE or taxonomy) on all the reads of a sample, without the ARG step;

### 4.1. File Structure

```
modular/
├── workflow/
│   ├── Snakefile              # workflow for the modular scripts
│   └── envs/                  # conda environments (LAST, seqtk, seqkit, Kraken2, Bracken, process)
└── scripts/
    ├── arg_analysis.py        # ARG summary + flanking-fragment extraction
    ├── mge_analysis.py        # IS, Tn3, Tn7, Tn554 and integron annotation
    ├── plasmid_analysis.py    # plasmid hits per read
    ├── taxonomy_analysis.py   # Kraken2 lineage per read
    ├── merge_summaries.py     # joins module outputs on query_name
```

The modular pipeline has its own scripts and conda environments. It shares only `samples/`, `database/` and `metadata/` with the main pipeline.

### 4.2. Running the modular workflow

Run from the repository root, with `-s` pointing to the modular Snakefile:
```
$ snakemake -s modular/workflow/Snakefile --use-conda --conda-frontend conda --cores all
```

By default this runs the full pipeline. To run only part of the workflow, add a target name to the end of the command:

| Target | What it runs | Main output |
|---|---|---|
| `all` (default) | Full pipeline| `results_v1/final_summary/{sample}_final_summary.csv`, `results_v1/kraken/arg_reads/{sample}_genus.bracken` |
| `arg_only` | ARG identification and flanking-region extraction | `results_v1/arg_summary/{sample}_CARD_summary.csv`, `results_v1/arg_extracted/{sample}_CARD_last.fasta` |
| `mge_only` | MGE annotation of **all** reads (no ARG step) | `results_v1/mge_summary/reads/{sample}_mge_summary.csv` |
| `plasmid_only` | Plasmid classification of **all** reads | `results_v1/plasmid_summary/reads/{sample}_plasmid_summary.csv` |
| `taxonomy_only` | Kraken2 + Bracken on **all** reads | `results_v1/taxonomy_summary/reads/{sample}_taxonomy_summary.csv`, `results_v1/kraken/reads/{sample}_genus.bracken` |

For example:
```
$ snakemake -s modular/workflow/Snakefile --use-conda --conda-frontend conda --cores all mge_only
```


**Outputs of the full pipeline**, in order:

| Output | Contents |
|---|---|
| `results_v1/arg_summary/{sample}_CARD_summary.csv` | ARG hits, one row per ARG-carrying read |
| `results_v1/all_summary/{sample}_all_summary.csv` | + IS + Tn3, Tn7, Tn554 and integron columns |
| `results_v1/final_summary/{sample}_final_summary.csv` | + plasmid and taxonomy columns |

Logs for every rule are written to `results_v1/logs/`.

## 5. Downstream analysis and figures

All the analyses and figures in the study can be reproduced with [analysis_visualizations.ipynb](analysis_visualizations.ipynb). The notebook takes the final summaries from the pipeline. To fully replicate the results, the annotated data from our run can be found in Zenodo  [https://doi.org/10.5281/zenodo.23492757](https://doi.org/10.5281/zenodo.23492757).

For comparison, we downloaded the ProMGE and PLSDB data from their repositories and matched them with our results in the notebook:

- **ProMGE**: we mapped the ARO accessions of its ARG–MGE associations to CARD gene families in the latest CARD version. We then compared them with the genus, ARG and MGE combinations found in our samples.
- **PLSDB (v.2024 05 31 v2)**: we used only the AMR annotations made with RGI (CARD). We compared them with the ARGs on our putative plasmid reads, by gene family.

To run the notebook on your own samples, provide your own sample metadata in the same format as ours. For each sample, give its accession, source, average read length and total number of bases, as well as its number of reads. The normalization and the per-source comparisons use these values. Also update the file paths in the notebook to point to your own files.
