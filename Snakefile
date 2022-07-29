configfile: "config.yaml"
import pandas as pd
samples_data = pd.read_csv(config['samples'], sep =";")

#SAMPLE = samples_data["sample"]
ACCESSIONS = samples_data["accession"]
SAMPLES = samples_data["samples"]

rule all:
    input:
        expand("02_mapped/{srr}.bam", srr = ACCESSIONS),
        expand("04_bigWigFiles/{sample}.bw", sample = SAMPLES),
        expand("03_calledPeaks/{sample}_summits.bed", sample = SAMPLES[SAMPLES != "input"])


def getTreatment(wildcards):
    return expand("02_mapped/{srr}.bam",
    srr = list(samples_data[samples_data["samples"]==wildcards.sample]["accession"]))
def getSample(wildcards):
    return expand("04_sorted/{srr}.bam",
    srr = list(samples_data[samples_data["samples"]==wildcards.sample]["accession"]))
def getIndex(wildcards):
    return expand("04_sorted/{srr}.bam.bai",
    srr = list(samples_data[samples_data["samples"]==wildcards.sample]["accession"]))


rule get_index_files:
    params:
        genomeLocation=config['genomeFTP']
    output:
        "genome.fa",

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
    shell: "bowtie2-build --threads {threads} {input} {params.out}"


rule prefetch:
    output:
        "01_raw/.prefetch/sra/{srr}.sra"
    params:
        "{srr} --max-size 50GB -O 01_raw"
    log:
        "01_raw/.prefetch/sra/{srr}.log"
    conda:
        "sra_chipseq.yaml"
    shell:
        """
        prefetch {params} > {log} 2>&1 && touch {output}
        """

rule fastqdump:
    input:
        "01_raw/.prefetch/sra/{srr}.sra"
    output:
        "01_raw/{srr}.fastq"
    params:
        args = "-S -O 01_raw/ -t 01_raw/",
        id_srr = "{srr}"
    log:
        "01_raw/{srr}.log"
    conda:
        "sra_chipseq.yaml"
    threads: 8
    shell:
        "(fasterq-dump {params.args} {params.id_srr} -e {threads})2> {log}"

rule bowtie2_map:
    input:
        read = "01_raw/{srr}.fastq",
        genome = "genome.1.bt2"
    params:
        genome = "genome"
    output:
        "02_mapped/{srr}.sam"
    conda:
        "sra_chipseq.yaml"
    log:
        "02_mapped/{srr}.log"
    threads: 8
    shell:
        "(bowtie2 -p {threads} -x {params.genome} -U {input.read} -S {output})2> {log}"

rule get_BAM:
    input:
        "02_mapped/{srr}.sam"
    output:
        "02_mapped/{srr}.bam"
    conda:
        "sra_chipseq.yaml"
    threads: 8
    shell:
        "samtools view --threads {threads} -bq 1 {input} > {output}"
rule samtools_sort:
    input:
        "02_mapped/{srr}.bam"
    output:
        "04_sorted/{srr}.bam"
    conda:
        "sra_chipseq.yaml"
    shell:
    # why wildcards.sample? and not sample
        "samtools sort -T 04_sorted/{wildcards.srr} "
        "-O bam {input} > {output}"

rule samtools_index:
    input:
        "04_sorted/{srr}.bam"
    output:
        "04_sorted/{srr}.bam.bai"
    shell:
        "samtools index {input}"

rule call_peaks:
    input:
        getTreatment
    output:
        "03_calledPeaks/{sample}_summits.bed"
    log:
        "03_calledPeaks/{sample}.log"
    params:
        control= list(samples_data[samples_data["samples"] =="input"]["accession"]),
        input = "-f BAM -g ",
        genome_size = config['genome_size'],
        qCutOff = config["q_cutOff_peakCall"]
    threads:2
    shell:
        "(macs2 callpeak -t {input} -c 02_mapped/{params.control}.bam {params.input} {params.genome_size} --outdir 03_calledPeaks/ -n {wildcards.sample} -q {params.qCutOff})2> log"

rule make_bigWig:
    input:
        sample= getSample,
        index = getIndex

    output:
        "04_bigWigFiles/{sample}.bw"
    params:
        normalization = config["bigWig_normalization"],
        binsize= config["binSize"]
    threads: 6
    log:
        "04_bigWigFiles/{sample}.log"
    shell:
        "(bamCoverage -b {input.sample} --normalizeUsing {params.normalization} --binSize {params.binsize} -o {output} --numberOfProcessors {threads})2> log"
