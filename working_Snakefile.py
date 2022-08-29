configfile: "config.yaml"
import pandas as pd
from itertools import combinations
samples_data = pd.read_csv(config['samples'], sep =";")

antibodies_oi = ["H3K4me3", "H3K27me3", "ChIP-Seq input"]
antibodies_used = ["H3K4me3", "H3K27me3"]
antibodies_used = config["antibodies"]
antibodies_combination = [*combinations(antibodies_used, 2)]
antibodies_combination = pd.DataFrame(antibodies_combination, columns =['AB1','AB2'])
ab1 = antibodies_combination.AB1
ab2 = antibodies_combination.AB2
combos =ab1 +"_" +ab2

###add this to code to reduce samples and to tidy up table
samples_data = samples_data.loc[samples_data["Experiment"].isin(antibodies_oi)]
samples_data['Sample Name'] = samples_data["Sample Name"].str.replace(',','_')
samples_data['Sample Name'] = samples_data["Sample Name"].str.replace('/s','_')


ACCESSIONS = samples_data["SRR"]
SAMPLES = samples_data["Sample Name"]
ANTIBODIES = samples_data["Experiment"]
LAYOUT = samples_data["modus"]
BIOSAMPLE = samples_data["bio_sample"]
BIOPROJECT = samples_data["project"]

sp = BIOSAMPLE + "_" + BIOPROJECT
sa = ANTIBODIES + "/" + ACCESSIONS
sl = LAYOUT +"/" + ACCESSIONS
#localrules: prefetch, get_index_files, fastqdump
samples = SAMPLES.drop_duplicates()
antibodies_combination['samples']= [samples for x in antibodies_combination['AB1']]
antibodies_combination = antibodies_combination.explode('samples')

allABS = ANTIBODIES[ANTIBODIES != "ChIP-Seq input"] + combos
allSAMP = ACCESSIONS[ANTIBODIES != "ChIP-Seq input"] + antibodies_combination['samples']
ALL_PEAKS = allABS +"/" + allSAMP

rule all:
    input:
        expand("02_mapped/sorted/{srr}.bam.bai", srr = ACCESSIONS),
        expand("02_mapped/input/{bio_id}.bam", bio_id = sp[ANTIBODIES =="ChIP-Seq input"]),
        expand("03_calledPeaks/{sa}_summits.bed", sa = sa[ANTIBODIES != "ChIP-Seq input"]),
        expand("04_bigWigFiles/{srr}.bw",  srr = ACCESSIONS[ANTIBODIES =="ChIP-Seq input"]),

 ## add code to file
        expand("05_quantified_signal/mono/{abs}_quantified.csv", abs = ANTIBODIES.drop_duplicates()),
        expand("06_merged_Samples/{ab}/{samples}.bw", ab = ANTIBODIES[ANTIBODIES != "ChIP-Seq input"].drop_duplicates(), samples = SAMPLES.drop_duplicates()),
        expand("03b_calledPeaks/{combination}/{samples}_summits.bed", combination = combos, samples = SAMPLES.drop_duplicates()),
        expand("05b_quantified_signal/{ab_combo}_quantified.csv", ab_combo = allABS)

### add code to file
def getPeakFile_combo(wildcards):
    return expand("03b_calledPeaks/{ab_combo}/{samples}_summits.bed", ab_combo = wildcards.ab_combo, samples =SAMPLES.drop_duplicates())


def getFile1(wildcards):
    return expand("06_merged_Samples/{ab}/{samples}.bw", ab = list(ab1[combos == wildcards.abcombo]), samples = wildcards.samples)

def getFile2(wildcards):
    return expand("06_merged_Samples/{ab}/{samples}.bw", ab = list(ab2[combos == wildcards.abcombo]), samples = wildcards.samples)


def Merge_SAMPLES(wildcards):
    return expand("02_mapped/sorted/{srr}.bam", srr = ACCESSIONS.loc[(ANTIBODIES == wildcards.ab) & (SAMPLES == wildcards.samples)])

def getPeakFile_single(wildcards):
    return expand("03_calledPeaks/{sa}_summits.bed", sa = sa[ANTIBODIES ==wildcards.abs])

def getBW(wildcards):
    return expand("04_bigWigFiles/{srr}.bw",  srr = ACCESSIONS[ANTIBODIES =="ChIP-Seq input"])



#### modify codes
def get_Sam(wildcards):
      return expand("02_mapped/{layout_srr}.sam", layout_srr = list(sl[ACCESSIONS==wildcards.srr]))


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
        "prefetch --ngc prj_phs710EA_test.ngc {params}"

rule prefetch_PAIRED:
    output:
        temp("01_raw/PAIRED/{srr}/{srr}.sra")
    params:
        "{srr} -O 01_raw/PAIRED"
    conda:
        "sra_chipseq.yaml"
    shell:
        "prefetch --ngc prj_phs710EA_test.ngc {params}"

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
    threads: 8
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
    threads: 8
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
    shell:
        "bowtie2 -p {threads} -x genome -1 {input.read_1} -2 {input.read_2} -S {output}"

rule filterSAM_paired:
    input:
        "02_mapped/PAIRED/{srr}.sam"
    output:
        temp("02_mapped/filtered/{srr}.sam")
    conda:
        "sra_chipseq.yaml"
    threads: 8
    shell:
        "samtools view --threads {threads} -S -h {input} | grep -v 'XS:i:' > {output}"

rule filterSAM_single:
    input:
        "02_mapped/SINGLE/{srr}.sam"
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

rule call_peaks:
    input:
        treat = getTreat,
        ctrl = getCtrl
    output:
        protected("03_calledPeaks/{sa}_summits.bed")
    conda:
        "sra_chipseq.yaml"
    params:
        input = "-f BAM -g ",
        genome_size = config['genome_size'],
        qCutOff = config["q_cutOff_peakCall"]
    threads:2
    shell:
        "macs2 callpeak -t {input.treat} -c {input.ctrl} {params.input} {params.genome_size} --outdir 03_calledPeaks/ -n {wildcards.sa} -q {params.qCutOff}"

rule make_bigWig:
    input:
        sample= "02_mapped/sorted/{srr}.bam",
        index = "02_mapped/sorted/{srr}.bam.bai"
    output:
        "04_bigWigFiles/{srr}.bw"
    conda:
        "sra_chipseq.yaml"
    params:
        normalization = config["bigWig_normalization"],
        binsize= config["binSize"]
    threads: 8
    shell:
        "bamCoverage -b {input.sample} --normalizeUsing {params.normalization} --binSize {params.binsize} -o {output} --numberOfProcessors {threads}"

rule quantify_peaks_single:
    input:
        peaks= getPeakFile_single,
        coverage = getBW
    output:
        "05_quantified_signal/mono/{abs}_quantified.csv"
    threads: 2
    conda:
        "sra_chipseq.yaml"
    shell:
        "Rscript scripts/Quantify_Regions.s 03_calledPeaks/{wildcards.abs}"

rule Merge_BWs:
    input:
        Merge_SAMPLES
    output:
        "06_merged_Samples/{ab}/{samples}.bw"
    conda:
        "sra_chipseq.yaml"
    params:
        normalization = config["bigWig_normalization"],
        binsize= config["binSize"]
    threads: 8
    shell:
        "samtools merge --threads {threads} -o 06_merged_Samples/{wildcards.ab}/{wildcards.samples}_unsorted.bam {input}"
        "samtools sort -T 06_merged_Samples/{wildcards.ab}/{wildcards.samples} -O 06_merged_Samples/{wildcards.ab}/{wildcards.samples}_unsorted.bam > 06_merged_Samples/{wildcards.ab}/{wildcards.samples}_sorted.bam"
        "samtools index 06_merged_Samples/{wildcards.ab}/{wildcards.samples}_sorted.bam"
        "bamCoverage -b 06_merged_Samples/{wildcards.ab}/{wildcards.samples}_sorted.bam --normalizeUsing {params.normalization} --binSize {params.binsize} -o {output} --numberOfProcessors {threads}"
        "rm -f 06_merged_Samples/{wildcards.ab}/{wildcards.samples}_unsorted.bam"
        "rm -f 06_merged_Samples/{wildcards.ab}/{wildcards.samples}_sorted.bam"
        "rm -f 06_merged_Samples/{wildcards.ab}/{wildcards.samples}_sorted.bam.bai"

rule find_overlaps:
    input:
        file1= getFile1,
        file2= getFile2
    output:
        "03b_calledPeaks/{abcombo}/{samples}_summits.bed"
    conda:
        "sra_chipseq.yaml"
    params:
        genome_size = config['genome_size'],
        qCutOff = config["q_cutOff_peakCall"]
    threads: 2
    shell:
        "bedtools unionbedg -i {input.file1} {input.file2} > 06_merged_Samples/{wildcards.abcombo}/{wildcards.samples}.bedGraph"
        "python scripts/takeLower.py 06_merged_Samples/{wildcards.abcombo}/{wildcards.samples}.bedGraph 06_merged_Samples/{wildcards.abcombo}/{wildcards.samples}.bw"
        "macs2 bdgpeakcall -i 06_merged_Samples/{wildcards.abcombo}/{wildcards.samples}.bw -c {params.qCutOff} -o 03_calledPeaks/{wildcards.abcombo}/{wildcards.samples}"

rule quantify_peaks_combo:
    input:
        peaks= getPeakFile_combo,
        coverage = getBW
    output:
        "05b_quantified_signal/{ab_combo}_quantified.csv"
    conda:
        "sra_chipseq.yaml"
    threads: 2
    shell:
        "Rscript scripts/Quantify_Regions.s 03b_calledPeaks/{wildcards.ab_combo}"
