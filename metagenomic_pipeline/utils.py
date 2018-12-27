# utility functions
import os
from Bio import Entrez
import xml.etree.ElementTree as ET

from config import *


# download_sra_data.py
OUTPUT_PATH = "SRA_data/"
## Optional. You can either put the project id here or through the command line
INPUT_KEY = "PRJNA298936"
# INPUT_KEY = ["PRJNA298936"]
# INPUT_KEY = 'SRS1435466'
# INPUT_KEY = ["SRS1107428", "SRS1107429"]


# blast_preprocess.py
FASTA_PATH = "SRA_data/fasta/"
OUTPUT_DIR = "blast/"
## output format
OUTFMT = '6 qseqid sseqid sscinames pident length mismatch gapopen qstart qend sstart send evalue bitscore'


# postprocessing.py
# TOP_N = 1
PROTEIN_BATCH_SIZE = "all"
## these blast format have to be the same as above in order to process to results
BLAST_FORMAT = "qseqid sseqid sscinames pident length mismatch gapopen qstart qend sstart send evalue bitscore"


# combine_result.py
ANALYSIS_PATH = "analysis/"
COMBINE_OUTPUT_NAME = ""


Entrez.email = EMAIL

def process_input(INPUT_KEY):
    blast_flag = False

    if type(INPUT_KEY) != list:
        INPUT_KEY = [x.strip() for x in INPUT_KEY.split(',')]

    temp = []
    for key in INPUT_KEY:
        temp += key.split(',')

    INPUT_KEY = [key for key in temp if key != ""]
    if "blast" in INPUT_KEY:
        blast_flag = True
        INPUT_KEY.remove('blast')

    INPUT_KEY = [key.replace(",", "") for key in INPUT_KEY]

    return INPUT_KEY, blast_flag


def _classify_single_input_type(input_keyword):
    if "PRJNA" in input_keyword:
        return "bioproject"
    else:
        return "sra"

def classify_input_type(input_keyword):
    if (type(input_keyword) == list):
        input_type = [_classify_single_input_type(key) for key in input_keyword]
        if len(set(input_type)) == 1:
            if input_type[0] == "sra":
                return "sra_list"
            if input_type[0] == "bioproject":
                return "bioproject_list"
        else:
            return "mix_list"
    else:
        return _classify_single_input_type(input_keyword)

def _get_sra_from_bioproject(bioproject_id):
    out = os.popen(f"wget -qO- 'http://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?save=efetch&db=sra&rettype=runinfo&term={bioproject_id}' | grep {bioproject_id} | cut -f25 -d','").read()
    return list(set(out.split('\n')[:-1]))

def get_sra_list(INPUT_KEY):
    input_type = classify_input_type(INPUT_KEY)

    if input_type == "bioproject":
        sra_list = _get_sra_from_bioproject(INPUT_KEY)
        print(f"Input is bioproject. Got {len(sra_list)} related sras.")
    elif input_type == "sra_list":
        sra_list = INPUT_KEY
        print(f"Input is list of sra. Got {len(sra_list)} sras.")
    elif input_type == "bioproject_list":
        if len(INPUT_KEY) == 1:
            sra_list = _get_sra_from_bioproject(INPUT_KEY[0])
            print(f"Input is bioproject. Got {len(sra_list)} related sras.")
        else:
            print("Not support multiple bioproject! Please input one by one")
    else:
        sra_list = [INPUT_KEY]

    return sra_list

def get_latest_run_version(runs):
    print("run ", runs)
    if(len(runs) > 1):
        latest_date = ""
        for run_id in runs:
            handle = Entrez.efetch(db="sra", id=run_id, rettype="gb", retmode="xml")
            root = ET.fromstring(handle.read())
            published_date = root.find('EXPERIMENT_PACKAGE/RUN_SET/RUN').attrib['published']
            print(f'run id: {run_id} is published on {published_date}')
            if(published_date>latest_date):
                lastest_run = run_id
                latest_date = published_date
        print(f'\t\tselected {lastest_run}')
    else:
        lastest_run = runs[0]
    return lastest_run

