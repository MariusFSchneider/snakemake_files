configfile: "config.yaml"
import pandas as pd
from itertools import combinations
samples_data = pd.read_csv(config['samples'], sep =";")

antibodies_oi = ["H3K4me3", "H3K27me3", ]
antibodies_used = ["H3K4me3", "H3K27me3"]
#antibodies_used = config["antibodies"]
antibodies_combination = [*combinations(antibodies_used, 2)]


###add this to code to reduce samples and to tidy up table
samples_data = samples_data.loc[samples_data["Exeriment"].isin(antibodies_oi)]
samples_data['Sample Name'] = samples_data["Sample Name"].str.replace(',','_')
samples_data['Sample Name'] = samples_data["Sample Name"].str.replace('\s','_')


ACCESSIONS = samples_data["SRR"]
SAMPLES = samples_data["Sample"]
ANTIBODIES = samples_data["Experiment"]
LAYOUT = samples_data["modus"]
BIOSAMPLE = samples_data["biosample"]
BIOPROJECT = samples_data["project"]

sp = samples_data["bio_sample"]+ "_" + samples_data["project"]
sa = samples_data["Experiment"] + "/" + samples_data["SRR"]

#localrules: prefetch, get_index_files, fastqdump

rule all:
    input:
        expand("02_mapped/sorted/{srr}.bam.bai", srr = ACCESSIONS),
        expand("02_mapped/input/{bio_id}.bam", bio_id = sp[ANTIBODIES =="ChIP-Seq input"]),
        expand("03_calledPeaks/{sa}_summits.bed", sa = sa[ANTIBODIES != "ChIP-Seq input"]),
        expand("04_bigWigFiles/{srr}.bw",  srr = ACCESSIONS[ANTIBODIES =="ChIP-Seq input"]),
        
### add code to file        
        expand("{AB}_quantified.csv", AB = total_abs.drop.duplicates)
        

def get_Sam(wildcards):
   return expand("02_mapped/{layout}/{wildcards.srr}.sam", layout = list(LAYOUT[samples_data['SRR']==wildcards.srr]))


def getINPUTS(wildcards):
    return expand("02_mapped/sorted/{srr}.bam", srr = ACCESSIONS.loc[(sp == wildcards.bio_id) & (ANTIBODIES == "ChIP-Seq input")])

def getTreat(wildcards):
   return expand("02_mapped/sorted/{srr}.bam", srr = list(ACCESSIONS[sa==wildcards.sa]))

def getCtrl(wildcards):
    return expand("02_mapped/input/{bio_id}.bam", bio_id = list(sp[sa==wildcards.sa]))

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
    threads: 24
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
    threads: 24
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
    shell:
        "bowtie2 -p {threads} -x genome -U {input.read} -S 02_mapped/SINGLE/{wildcards.srr}.sam)2> 02_mapped/SINGLE/{wildcards.srr}_alignment.txt"
        "samtools view --threads {threads} -S -h {input} | grep -v 'XS:i:' > 02_mapped/SINGLE/{wildcards.srr}_filtered.sam"
        "samtools view --threads {threads} -b {input} > {output}"    
        "rm -f 02_mapped/SINGLE/{wildcards.srr}.sam"
        "rm -f 02_mapped/SINGLE/{wildcards.srr}_filtered.sam"
        
        
rule bowtie2_map_PAIRED:
    input:
        read_1 = "01_raw/PAIRED/{srr}_1.fastq",
        read_2 = "01_raw/PAIRED/{srr}_2.fastq",
        genome = "genome.1.bt2"
    output:
        temp("02_mapped/PAIRED/{srr}.bam")
    conda:
        "sra_chipseq.yaml"
    threads: 24
    shell:
        "(bowtie2 -p {threads} -x genome -1 {input.read_1} -2 {input.read_2} -S 02_mapped/PAIRED/{wildcards.srr}.sam)2> 02_mapped/PAIRED/{wildcards.srr}_alignment.txt"
        "samtools view --threads {threads} -S -h {input} | grep -v 'XS:i:' > 02_mapped/PAIRED/{wildcards.srr}_filtered.sam"
        "samtools view --threads {threads} -b {input} > {output}"    
        "rm -f 02_mapped/PAIRED/{wildcards.srr}.sam"
        "rm -f 02_mapped/PAIRED/{wildcards.srr}_filtered.sam"


rule samtools_sort_index:
    input:
        get_Bam
    output:
        BAM = protected("02_mapped/sorted/{srr}.bam")
    conda:
        "sra_chipseq.yaml"
    shell:
    # why wildcards.sample? and not sample
        "samtools sort -T 02_mapped/sorted/{wildcards.srr} -O bam {input} > {output}"
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

rule call_peaks:
    input:
        treat = getTreat,
        ctrl = getCtrl
    output:
        protected("03_calledPeaks/{sa}_summits.bed")
    params:
        input = "-f BAM -g ",
        genome_size = config['genome_size'],
        qCutOff = config["q_cutOff_peakCall"]
    threads:2
    shell:
        "macs2 callpeak -t {input.treat} -c {input.ctrl} {params.input} {params.genome_size} --outdir 03_calledPeaks/ -n {wildcards.sa} -q {params.qCutOff})"

rule make_bigWig:
    input:
        sample= "02_mapped/sorted/{srr}.bam",
        index = "02_mapped/sorted/{srr}.bam.bai"

    output:
        "04_bigWigFiles/{srr}.bw"
    params:
        normalization = config["bigWig_normalization"],
        binsize= config["binSize"]
    threads: 8
    shell:
        "bamCoverage -b {input.sample} --normalizeUsing {params.normalization} --binSize {params.binsize} -o {output} --numberOfProcessors {threads})"
    
