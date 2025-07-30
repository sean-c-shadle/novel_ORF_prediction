# Pipeline for novel (unannotated) translated ORF discovery with Ribotricer
## Dependencies
1. conda
2. snakemake (tested on version 9.3.3)

## Detailed instructions for use
1. get scripts `git clone https://github.com/sean-c-shadle/novel_ORF_prediction.git`
2. modify working directory config.yaml file to point to gtfs, fasta files, reference bed.
3. generate merged bam files (e.g. from duplicate ribolite embryo developmental timepoints). See script in scripts/merge_bams.sh for an example
4. create a samples.txt file (see example) with filenames minus suffix
5. Run with e.g. `snakemake --profile granite --use-conda  --jobs 100`

## Required file inputs
1. gtf files generated from stringtie output and a concatenated Ensembl gtf file with the stringtie output (adding more info on how this was generated later).
2. genome fasta file.
3. reference exon bed file (I used gencode downloaded from UCSC genome table browser).


## snakemake profile for slurm use
Here is an example of my config.yaml file in ~/.config/snakemake/granite
```
executor: slurm

latency-wait: 60 #seconds

default-resources:
    qos: "cairns-grn"
    slurm_partition: "cairns-grn"
    slurm_account: "cairns"
    runtime: 720  # minutes
```
## Snakemake workflow DAG
![DAG](https://github.com/sean-c-shadle/novel_ORF_prediction/blob/main/dag2.png)
