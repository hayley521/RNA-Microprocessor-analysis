# RNA Microprocessor sequencing analysis

Code used to process and analyze the sequencing datasets associated with the
RNA Microprocessor study. The repository is organized by assay and target so
that each deposited GEO dataset can be traced to its primary processing code.

## Repository layout

- `original_scripts/bulk_RNAseq/`: paired-end RNA-seq preprocessing,
  alignment, gene quantification, and rMATS alternative-splicing workflows.
- `original_scripts/small_RNAseq/`: small-RNA preprocessing and quantification,
  plus DESeq2 notebooks for DROSHA, DGCR8, and DICER1 experiments.
- `original_scripts/GROseq/`: GRO-seq preprocessing, rRNA depletion,
  alignment, HOMER counting, and browser-track generation.
- `scripts/`: deterministic export and validation utilities used to prepare
  the processed tables deposited with GEO.

## Datasets

The repository covers newly generated bulk RNA-seq, small RNA-seq, and GRO-seq
experiments for DROSHA, DGCR8, DICER1, and matched controls. GEO accession
numbers will be added after submission processing is complete.

## References and main software

- Organism: *Mus musculus*
- Genome assembly: GRCm38/mm10
- Gene annotation: GENCODE vM25
- Core command-line tools: Trim Galore, Cutadapt, STAR, Bowtie/Bowtie2,
  SAMtools, featureCounts, rMATS, and HOMER
- Main R packages: DESeq2, tidyverse, EnhancedVolcano, clusterProfiler,
  DOSE, org.Mm.eg.db, enrichplot, stringr, and repr

## Reuse

The shell and Slurm files preserve the commands and parameters used for the
original analyses. Paths beginning with `/home/` or `/wanglab1/` are
site-specific historical paths and must be replaced with local FASTQ,
reference, and output locations before reuse. Sample-to-condition mappings are
documented by the GEO metadata deposited with the sequencing data.

The Jupyter notebooks are distributed without executed outputs to keep the
repository compact and prevent result duplication. Processed count and
differential-analysis tables are deposited in GEO rather than committed here.

## Data and code availability

Sequencing data: NCBI Gene Expression Omnibus, accession numbers pending.

Code archive DOI: pending creation of the versioned Zenodo release.
