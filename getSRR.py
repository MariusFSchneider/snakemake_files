#import neccessary module
import pandas as pd
import subprocess
import os
import re
import glob


#define functions
def getSRR(accession):
    ''' parameters:
    input: accession string
    output: SRR accesions as list

    read fetched accession data and reads SRR accesion'''
    input = accession + ".txt"
    input_data = pd.read_csv(input, delimiter = ",")
    srr = input_data['Run']
    modus = input_data['LibraryLayout']
    project = input_data["BioProject"]
    sample_type = input_data["Sample"]
    sample_bio = input_data["BioSample"]
    lib_size = input_data["size"]
    return list(srr+ "_" + modus + '_' + project + '_' + sample_type + '_' + sample_bio + '_' + lib_size )



def getRUN_info(accession):
    ''' parameters:
    input: accession as string
    output: returns na or runs getSRR function
    error: Value error
    interface to run esearch program, save fetched data as temporary file, runs getSRR function and remove temporary file'''
    path = accession + ".txt"
    subprocess.call(("esearch -db sra -query "+ accession+ " | efetch -format runinfo >" + path), shell = True)
    if os.stat(path).st_size == 0:
        return "na"
    else:
        return getSRR(accession)

    subprocess.call(("rm -f " +  path), shell = True)


def getAccession(link):
    '''parameters
    input: link as string
    output: returns accession as string
    find by Reges for SRX string '''
    return re.search('SRX[0-9]+', link).group(0)


def getRUN_info_bySRA(link):
    ''' parameters:
    input: link including accesion as string
    output: returns na or runs getSRR function
    error: Value error

    interface to run esearch program, save fetched data as temporary file, runs getSRR function'''
    try:
        accession = getAccession(link)
        path = accession + ".txt"
        if not os.path.exists(path):
            subprocess.call(("esearch -db sra -query "+ accession+ " | efetch -format runinfo >" + path), shell = True)
        return getSRR(accession)
    except:
        return "na"

def getRUN_info_byGEO(geo):
    ''' parameters:
    input: accession as string
    output: returns na or runs getSRR function
    error: Value error

    interface to run esearch program, save fetched data as temporary file, runs getSRR function'''
    try:
        accession2 = geo
        path2 = accession2 + ".txt"
        if not os.path.exists(path):
            subprocess.call(("esearch -db sra -query "+ accession2+ " | efetch -format runinfo >" + path), shell = True)
        return getSRR(accession2)
    except:
        return "na"


def getRUN_info(link_oi, geo_oi):
    '''parameters:
    input: lin as string or GEO accession as string
    output returns na or runs getSRR function

     try to get SRR accesion either by SRA or GEO
    '''
    srr_result = getRUN_info_bySRA(link_oi)
    if srr_result == "na":
        srr_result2 = getRUN_info_byGEO(geo_oi)
        if srr_result == "na":
            return "na"
        else:
            return srr_result2
    else:
        return srr_result


## create output folder

outdir = './output'
if not os.path.exists(outdir):
    os.mkdir(outdir)

## loop through samples files

SAMPLES = pd.DataFrame()
for file in glob.glob("*.csv"):
    print(file)
    samples  = pd.read_csv(file ,delimiter= ",")
    SAMPLES = SAMPLES.append(samples)

SAMPLES = SAMPLES.drop_duplicates()


SAMPLES['SRR'] = [getRUN_info(x,y) for x, y in zip(SAMPLES["SRA Run Selector"], SAMPLES["# GEO Accession"])]


SAMPLES.to_csv("output/all_samples.csv", sep = ";")
##conver list in col to multiple rows
SAMPLES =SAMPLES.explode('SRR')

## remove na
SAMPLES = SAMPLES[SAMPLES['SRR'] != "na"]
## split SRR to SRR and mode column:

SAMPLES[['SRR', 'modus','project','sample_id','bio_sample','library_size']] = SAMPLES['SRR'].str.split('_', expand=True)
SAMPLES = SAMPLES.drop_duplicates()

SAMPLES.to_csv("output/unique_samples.txt", sep = ";")



data_selected = SAMPLES[SAMPLES['Experiment'].isin(["ChIP-Seq input","H3K4me3","H3K27me3"])]


test_data = data_selected[data_selected['Sample Name'].isin(["spleen","neurosphere cultured cells, ganglionic eminence derived"])]

test_data.to_csv("output/test_samples.csv", sep = ";")
