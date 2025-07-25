# Pipeline for novel (unannotated) translated ORF discovery with Ribotricer
## Dependencies
1. modules (version used)
    - bedtools v2.28.0
    - cufflinks v2.2.1
    - InterProScan version 5.63-95.0
2. [Transdecoder](https://github.com/TransDecoder/TransDecoder/releases) v5.7.1
3. conda
4. snakemake (tested on version 9.3.3)

## Detailed instructions for use
1. get scripts `git clone https://github.com/sean-c-shadle/novel_ORF_prediction.git`
2. unzip Trandecoder (see above) in the working directory
3. modify working directory config.yaml file to point to gtfs, fasta files, reference bed.
4. we merged bam files from duplicate ribolite embryo developmental timepoints. See script in scripts/merge_bams.sh
5. create a samples.txt file (see example) with filenames minus suffix
6. Run with `snakemake --profile granite --use-conda  --jobs 100`

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
