
configfile: "config.yaml"

#prefix for bam files, one per line.
with open("samples.txt") as file:
    samples = [line.strip() for line in file if line.strip()] 

rule all:
    input:
        expand("{sample}/{sample}_domains.tsv", sample=samples),

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
rule make_bedfile:
    message: "making exon bed file for candidate ORFs"
    input:
        tsv = rules.detect_translating_ORFs.output,
        candidates = rules.prepare_ORFs.output
    output:
        bed = "{sample}/{sample}_candidates.bed"
    shell:
        """
        python scripts/get_longest_MSTRG.py {input.tsv} > {wildcards.sample}/{wildcards.sample}_longest_MSTRG.txt
        grep -Ff {wildcards.sample}/{wildcards.sample}_longest_MSTRG.txt {input.candidates} > {wildcards.sample}/{wildcards.sample}_candidates.txt
        python scripts/make_bed.py {wildcards.sample}/{wildcards.sample}_candidates.txt > {wildcards.sample}/{wildcards.sample}_candidates.bed
        """ 

rule extract_unique_genes:
    message: "extracting genes with all unique exons (not overlapping your annotation e.g. gencode mm10 exons)"
    input: 
        candidate_bed = rules.make_bedfile.output
    output:
        "{sample}/{sample}_unique_novel_unannotated_no_known_overlap.fasta"
    params:
        reference_bed = config["reference_bed"],
        stringtie_gtf = config["stringtie_gtf"],
        fasta = config["genome_fasta"] 
    shell:
        """
        module load bedtools
        module load cufflinks
        # -v    Only report those entries in A that have no overlap in B. Restricted by -f and -r.
        # -s    Force “strandedness”. That is, only report hits in B that overlap A on the same strand. By default, overlaps are reported without respect to strand.
        # -a    BAM/BED/GFF/VCF file “A”. Each feature in A is compared to B in search of overlaps. Use “stdin” if passing A with a UNIX pipe.
        bedtools intersect -v -s -a {input.candidate_bed} -b {params.reference_bed} | cut -f 4 | cut -f 1 -d _ | uniq | sort > {wildcards.sample}/{wildcards.sample}_unique_novel_unannotated.txt
        #report those that overlap A
        bedtools intersect -s -a {input.candidate_bed} -b {params.reference_bed} | cut -f 4 | cut -f 1 -d _ | uniq | sort > {wildcards.sample}/{wildcards.sample}_exon_overlap.txt 
        #print lines that are only found in file1 (*_unique_novel_unannotated_sorted.txt)
        comm -23  {wildcards.sample}/{wildcards.sample}_unique_novel_unannotated.txt {wildcards.sample}/{wildcards.sample}_exon_overlap.txt > {wildcards.sample}/{wildcards.sample}_unique_novel_unannotated_no_known_overlap.txt        
        #extract the novel genes from the stringtie MSTRGs in the gtf file:
        grep -Ff {wildcards.sample}/{wildcards.sample}_unique_novel_unannotated_no_known_overlap.txt {params.stringtie_gtf} > {wildcards.sample}/{wildcards.sample}_unique_novel_unannotated_no_known_overlap.gtf
        #extract fasta seqs
        gffread {wildcards.sample}/{wildcards.sample}_unique_novel_unannotated_no_known_overlap.gtf -g {params.fasta} -w {wildcards.sample}/{wildcards.sample}_unique_novel_unannotated_no_known_overlap.fasta
        """
rule transdecoder:
    message: "getting candidate novel/unannotated amino acid sequences"
    input:
        candidate_fasta = rules.extract_unique_genes.output
    output:
        "{sample}/{sample}_unique_novel_unannotated_no_known_overlap.fasta.transdecoder.pep"
    shell:
        """
        TransDecoder-TransDecoder-v5.7.1/TransDecoder.LongOrfs --complete_orfs_only -O {wildcards.sample} -t {input.candidate_fasta}
        TransDecoder-TransDecoder-v5.7.1/TransDecoder.Predict -O {wildcards.sample} -t {input.candidate_fasta}
        """
rule interproscan:
    message: "detecting candidate domains"
    input:
        pep = rules.transdecoder.output
    output:
        "{sample}/{sample}_domains.tsv"
    threads: 12
    resources:
        mem_mb = 64000
    shell:
        """
        #interproscan doesn't like * characters, need to remove them
        sed 's/\*//g' {input.pep} > {wildcards.sample}/{wildcards.sample}_fixed_unique_novel_unannotated_no_known_overlap.pep
        #running using only the pfam domain database-- can alter if desired with -appl 
        module load interproscan
        interproscan.sh --appl pfam -cpu 12 -i {wildcards.sample}/{wildcards.sample}_fixed_unique_novel_unannotated_no_known_overlap.pep -f tsv -o {wildcards.sample}/{wildcards.sample}_domains.tsv
        """
