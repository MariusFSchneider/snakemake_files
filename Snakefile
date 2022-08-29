configfile: "config.yaml"
import pandas as pd
from itertools import combinations
samples_data = pd.read_csv(config['samples'], sep =";")

antibodies_oi = ["H3K4me3", "H3K27me3", "ChIP-Seq input"]
antibodies_used = ["H3K4me3", "H3K27me3"]
#antibodies_used = config["antibodies"]
antibodies_combination = [*combinations(antibodies_used, 2)]
antibodies_combination = pd.DataFrame(antibodies_combination, columns =['AB1','AB2'])
ab1 = antibodies_combination.AB1
ab2 = antibodies_combination.AB2
combo =ab1 +"_" +ab2

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

sp = BIOSAMPLE + "_" + BIOPROJECT
sa = ANTIBODIES + "/" + ACCESSIONS

#localrules: prefetch, get_index_files, fastqdump

rule all:
    input:
        expand("02_mapped/sorted/{srr}.bam.bai", srr = ACCESSIONS),
        expand("02_mapped/input/{bio_id}.bam", bio_id = sp[ANTIBODIES =="ChIP-Seq input"]),
        expand("03_calledPeaks/{sa}_summits.bed", sa = sa[ANTIBODIES != "ChIP-Seq input"]),
        expand("04_bigWigFiles/{srr}.bw",  srr = ACCESSIONS[ANTIBODIES =="ChIP-Seq input"]),
        
 ## add code to file       
        expand("05_quantified_signal\{abs}_quantified.csv", abs = ANTIBODIES.drop.duplicates())       
#        expand("06_merged_Samples/{ab}/{samples}.bw", ab = ANTIBODIES.drop.duplicates(), samples = SAMPLES.drop.duplicates()"),
#        expand("03_calledPeaks/{abcombo}/{samples}_summits.bed", abcombo= combo, samples= SAMPLES.drop.duplicates()),
        expand("05_quantified_signal\{ab_combo}_quantified.csv", abcombo = combo)
        
### add code to file 
def getPeakFile_combo(wildcards):
    return expand("03_calledPeaks/{abcombo}/{samples}_summits.bed", abcombo = combo, samples =SAMPLES.drop.duplicates())


def getFile1(wildcards):
    return expand("06_merged_Samples/{ab}/{samples}.bw", ab = ab1[combo == wildcards.abcombo], samples = wildcards.samples2)

def getFile2(wildcards):
    return expand("06_merged_Samples/{ab}/{samples}.bw", ab = ab2[combo == wildcards.abcombo], samples = wildcards.samples2)


def Merge_SAMPLES(wildcards):
    return expand("02_mapped/sorted/{srr}.bam", srr = ACCESSIONS.loc[(ANTIBODIES == wildcards.ab) & (SAMPLES == wildcards.ab)])

def getPeakFile_single(wildcards):
    return expand(""03_calledPeaks/{sa}_summits.bed", sa = sa[ANTIBODIES ==wildcards.abs]")

def getBW(wildcards):
    return expand("04_bigWigFiles/{srr}.bw",  srr = ACCESSIONS[ANTIBODIES =="ChIP-Seq input"])    



#### modify codes
def get_Bam(wildcards):
   return expand("02_mapped/{layout}/{srr}.bam",  srr = ACCESSIONS[ANTIBODIES =="ChIP-Seq input"], layout = list(LAYOUT[ACCESSION==wildcards.srr]))


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
        BAMO = protected("02_mapped/sorted/{srr}.bam"),
        BAIO = protected("02_mapped/sorted/{srr}.bam.bai")
    conda:
        "sra_chipseq.yaml"
    shell:
        "samtools sort -T 02_mapped/sorted/{wildcards.srr} -O bam {input} > {output.BAMO}"
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
        "05_quantified_signal\{ab}_quantified.csv"
    threads: 2
    conda:
        "sra_chipseq.yaml"        
    shell:
        "Rscript scripts/Quantify_Regions.s {wildcards.ab}"
        
rule Merge_BWs:
    input: 
        Merge_SAMPLES
    output:
        temp("06_merged_Samples/{ab}/{samples}.bw")
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
        "03_calledPeaks/{abcombo}/{samples}_summits.bed"
    conda:
        "sra_chipseq.yaml"
    params:      
        genome_size = config['genome_size'],
        qCutOff = config["q_cutOff_peakCall"]
    threads: 2        
    shell:
        "bedtools unionbedg -i {input.file1} {input.file2} > 06_merged_Samples/{wildcards.abcombo}/{wildcards.samples}.bedGraph"
        "python scripts/takeLower.py 06_merged_Samples/{wildcards.abcombo}/{wildcards.samples}.bedGraph 06_merged_Samples/{wildcards.abcombo}/{wildcards.samples}.bw"
        "macs2 bdgpeakcall -i 06_merged_Samples/{wildcards.abcombo}/{wildcards.samples}.bw-c {params.qCutOff} -o 03_calledPeaks/{wildcards.abcombo}/{wildcards.samples}"
        
rule quantify_peaks_combo:
    input:
        peaks= getPeakFile_compo,
        coverage = getBW
    output:
        "05_quantified_signal\{ab_combo}_quantified.csv"
    conda:
        "sra_chipseq.yaml"        
    threads: 2
    shell:
        "Rscript scripts/Quantify_Regions.s {wildcards.ab_combo}"
    
