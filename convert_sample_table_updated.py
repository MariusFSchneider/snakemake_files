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


def getRUN_info_byLink(link):
    ''' parameters:
    input: accession as string
    output: returns na or runs getSRR function
    error: Value error

    interface to run esearch program, save fetched data as temporary file, runs getSRR function and remove temporary file   '''
    try:
        accession = getAccession(link)
        path = accession + ".txt"
        subprocess.call(("esearch -db sra -query "+ accession+ " | efetch -format runinfo >" + path), shell = True)
        return getSRR(accession)
        subprocess.call(("rm -f " +  path), shell = True)
    except:
        return "na"

## loop through samples files

SAMPLES = pd.DataFrame()
for file in glob.glob("*.csv"):
    print(file)
    samples  = pd.read_csv(file ,delimiter= ",")
    samples['SRR'] = [getRUN_info_byLink(x) for x in samples["SRA Run Selector"]]
    SAMPLES = SAMPLES.append(samples)

SAMPLES.to_csv("all_samples.txt", sep = "\t")
##conver list in col to multiple rows
SAMPLES =SAMPLES.explode('SRR')

## remove na
SAMPLES = SAMPLES[SAMPLES['SRR'] != "na"]
## split SRR to SRR and mode column:

SAMPLES[['SRR', 'modus']] = SAMPLES['SRR'].str.split('_', expand=True)
SAMPLES = SAMPLES.drop_duplicates()

SAMPLES.to_csv("unique_samples.txt", sep = "\t")
