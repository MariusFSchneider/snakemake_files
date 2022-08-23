configfile: "config.yaml"
import pandas as pd
samples_data = pd.read_csv(config['samples'], sep =";")

#SAMPLE = samples_data["sample"]
ACCESSIONS = samples_data["SRR"]
SAMPLES = samples_data["Sample Name"]
ANTIBODIES = samples_data["Experiment"]
LAYOUT = samples_data["modus"]
antibodies_used = config["antibodies"]
sl = samples_data["modus"]+ "/" + samples_data["SRR"]
biosample = samples_data["bio_sample"]
project = samples_data["bioproject"]
sp = biosample + '_' + project
sample_size = samples_data["library_size"]
#localrules: prefetch, get_index_files, fastqdump


rule all:
    input:
        #expand("02_mapped/{layout_srr}.sam", layout_srr= sl)
        #expand("02_mapped/SINGLE/{srr}.sam", srr = ACCESSIONS[LAYOUT =="SINGLE"]),
        #expand("02_mapped/PAIRED/{srr}.sam", srr = ACCESSIONS[LAYOUT =="PAIRED"])
        #expand("02_mapped/filtered/{srr}.sam", srr = ACCESSIONS)
        expand("02_mapped/sorted/{srr}.bam.bai", srr = ACCESSIONS),
        expand("02_mapped/input/{bio_id}.bam", bio_id = sp[ANTIBODIES =="ChIP-Seq input"])

def get_Sam(wildcards):
   return expand("02_mapped/{layout_srr}.sam", layout_srr = list(sl[samples_data['SRR']==wildcards.srr]))


def getINPUTS(wildcards):
    return expand("02_mapped/sorted/{srr}.bam", srr = ACCESSIONS.loc[(sp == wildcards.bio_id) & (ANTIBODIES == "ChIP-Seq input")])



rule get_index_files:
    params:
        genomeLocation=config['genomeFTP']
    output:
        "genome.fa"
    shell:
        """
        wget -O genome.fa.gz {params.genomeLocation}
        gunzip genome.fa.gz
        """

rule build_genome:
    input:
        "genome.fa"
    params:
        out = "genome"
    output:
        "genome.1.bt2"
    threads: 16
    resources:
        mem_mb=32000,
    shell: "bowtie2-build --threads {threads} {input} {params.out}"


rule prefetch_SINGLE:
    output:
        temp("01_raw/SINGLE/{srr}/{srr}.sra")
    params:
        "{srr} -O 01_raw/SINGLE"
    conda:
        "sra_chipseq.yaml"
    shell:
        "prefetch {params}"

rule prefetch_PAIRED:
    output:
        temp("01_raw/PAIRED/{srr}/{srr}.sra")
    params:
        "{srr} -O 01_raw/PAIRED"
    conda:
        "sra_chipseq.yaml"
    shell:
        "prefetch {params}"

rule fastqdump_SINGLE:
    input:
        "01_raw/SINGLE/{srr}/{srr}.sra"
    output:
        temp("01_raw/SINGLE/{srr}.fastq")
    params:
        args = "-S -O 01_raw/SINGLE/ -t 01_raw/SINGLE/",
        id_srr = "{srr}"
    conda:
        "sra_chipseq.yaml"
    threads: 16
    shell:
        "fasterq-dump {params.args} {params.id_srr} -e {threads}"


rule fastqdump_PAIRED:
    input:
        "01_raw/PAIRED/{srr}/{srr}.sra"
    output:
        temp("01_raw/PAIRED/{srr}_1.fastq"),
        temp("01_raw/PAIRED/{srr}_2.fastq")
    params:
        args = "-S -O 01_raw/PAIRED/ -t 01_raw/PAIRED/",
        id_srr = "{srr}"
    conda:
        "sra_chipseq.yaml"
    threads: 16
    shell:
        "fasterq-dump {params.args} {params.id_srr} -e {threads}"



rule bowtie2_map_SINGLE:
    input:
        read = "01_raw/SINGLE/{srr}.fastq",
        genome = "genome.1.bt2"
    output:
        temp("02_mapped/SINGLE/{srr}.sam")
    conda:
        "sra_chipseq.yaml"
    threads: 24
    resources:
        mem_mb=20000
    shell:
        "bowtie2 -p {threads} -x genome -U {input.read} -S {output}"

rule bowtie2_map_PAIRED:
    input:
        read_1 = "01_raw/PAIRED/{srr}_1.fastq",
        read_2 = "01_raw/PAIRED/{srr}_2.fastq",
        genome = "genome.1.bt2"
    output:
        temp("02_mapped/PAIRED/{srr}.sam")
    conda:
        "sra_chipseq.yaml"
    threads: 24
    resources:
        mem_mb=20000
    shell:
        "bowtie2 -p {threads} -x genome -1 {input.read_1} -2 {input.read_2} -S {output}"

rule filterSAM:
    input:
        "02_mapped/PAIRED/{srr}.sam"
    output:
        temp("02_mapped/filtered/{srr}.sam")
    conda:
        "sra_chipseq.yaml"
    threads: 8
    shell:
        "samtools view --threads {threads} -S -h {input} | grep -v 'XS:i:' > {output}"


rule get_BAM:
    input:
        "02_mapped/filtered/{srr}.sam"
    output:
        temp("02_mapped/BAM/{srr}.bam")
    conda:
        "sra_chipseq.yaml"
    threads: 8
    shell:
        "samtools view --threads {threads} -b {input} > {output}"

rule samtools_sort:
    input:
        get_Sam
    output:
        protected("02_mapped/sorted/{srr}.bam")
    conda:
        "sra_chipseq.yaml"
    shell:
    # why wildcards.sample? and not sample
        "samtools sort -T 02_mapped/sorted/{wildcards.srr} -O bam {input} > {output}"

rule samtools_index:
    input:
        "02_mapped/sorted/{srr}.bam"
    output:
        protected("02_mapped/sorted/{srr}.bam.bai")
    shell:
        "samtools index {input}"

rule merge_inputs:
    input:
        getINPUTS
    output:
        "02_mapped/input/{bio_id}.bam"
    conda:
        "sra_chipseq.yaml"
    threads: 8
    shell:
        "samtools merge --threads {threads} -o {output} {input}"
