
configfile: "config.yaml"

#prefix for bam files, one per line.
with open("samples.txt") as file:
    samples = [line.strip() for line in file if line.strip()] 

rule all:
    input:
        expand("{sample}/{sample}_stringent_translating_ORFs.tsv", sample=samples),

rule prepare_ORFs:
    message: "prepping ORF file with Ribotricer"
    input:
        config["candidate_gtf"]
    output:
        "Embryo_candidate_orfs.tsv"
    conda:
        "envs/ribotricer_env.yaml"
    params:
        fasta = config["genome_fasta"]
    resources:
        mem_mb = 64000
    shell:
        "ribotricer prepare-orfs --gtf {input} --fasta {params.fasta} --prefix Embryo"

rule detect_translating_ORFs:
    message: "detecting translating ORFs with Ribotricer"
    input:
        bam = "{sample}.bam",
        gtf = rules.prepare_ORFs.output
    output:
        "{sample}/{sample}_stringent_translating_ORFs.tsv"
    conda:
        "envs/ribotricer_env.yaml"
    resources:
        mem_mb = 64000
    shell: 
        """
        mkdir -p {wildcards.sample}
        ribotricer detect-orfs --bam {input.bam} --ribotricer_index {input.gtf} --prefix {wildcards.sample}/{wildcards.sample}_stringent --min_valid_codons 30 --min_read_density 0.3 --min_valid_codons_ratio 0.1
        """
