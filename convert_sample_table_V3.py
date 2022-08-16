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
    return list(srr+ "_" + modus)



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


def getRUN_info(link, geo):
    ''' parameters:
    input: accession as string
    output: returns na or runs getSRR function
    error: Value error

    interface to run esearch program, save fetched data as temporary file, runs getSRR function and remove temporary file   '''
    try:
        accession = getAccession(link)
        try:
            path = accession + ".txt"
            subprocess.call(("esearch -db sra -query "+ accession+ " | efetch -format runinfo >" + path), shell = True)
            return getSRR(accession)
            subprocess.call(("rm -f " +  path), shell = True)
        except:
            accesion= geo
            path = accession + ".txt"
            subprocess.call(("esearch -db sra -query "+ accession+ " | efetch -format runinfo >" + path), shell = True)
            return getSRR(accession)
        #    return(accesion)
            subprocess.call(("rm -f " +  path), shell = True)
    except:
        return "na"


## loop through samples files

SAMPLES = pd.DataFrame()
for file in glob.glob("*.csv"):
    print(file)
    samples  = pd.read_csv(file ,delimiter= ",")
    SAMPLES = SAMPLES.append(samples)

SAMPLES = SAMPLES.drop_duplicates()
#SAMPLES['SRR'] = [getRUN_info(x) for x in SAMPLES[["SRA Run Selector","GEO Accession"]]]

SAMPLES['SRR'] = [getRUN_info(x,y) for x, y in zip(SAMPLES["SRA Run Selector"], SAMPLES["# GEO Accession"])]


SAMPLES.to_csv("all_samples.txt", sep = "\t")
##conver list in col to multiple rows
SAMPLES =SAMPLES.explode('SRR')

## remove na
SAMPLES = SAMPLES[SAMPLES['SRR'] != "na"]
## split SRR to SRR and mode column:

SAMPLES[['SRR', 'modus']] = SAMPLES['SRR'].str.split('_', expand=True)
SAMPLES = SAMPLES.drop_duplicates()

SAMPLES.to_csv("unique_samples.txt", sep = "\t")

data_selected = SAMPLES[SAMPLES['Experiment'].isin(["ChIP-Seq input","H3K4me3","H3K27me3"])]


test_data = data_selected[data_selected['Sample Name'].isin(["IMR90 cell line","CD34 mobilized primary cells"])]

test_data.to_csv("test_samples.csv", sep = ";")
