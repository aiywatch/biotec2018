import os
import xml.etree.ElementTree as ET
from Bio import Entrez, SeqIO
from sys import argv
from datetime import datetime
import sys

import blast_preprocess
from config import *
from utils import *

""" 
This file is used to download SRA data from NCBI. 
Input: Bioproject ID, SRA ID, or List of SRA IDs.
Output: SRR data in both FastA and FastQ formats.
"""


def create_folder_and_path(INPUT_KEY, OUTPUT_PATH="SRA_data/", separate_bioproject=False):

    FASTQ_OUTPUT_PATH = OUTPUT_PATH + "fastq/"
    FASTA_OUTPUT_PATH = OUTPUT_PATH + "fasta/"

    if not os.path.exists(OUTPUT_PATH):
        os.mkdir(OUTPUT_PATH)
    if not os.path.exists(FASTQ_OUTPUT_PATH):
        os.mkdir(FASTQ_OUTPUT_PATH)
        os.system(f"touch {FASTQ_OUTPUT_PATH}fastq.log")
    if not os.path.exists(FASTA_OUTPUT_PATH):
        os.mkdir(FASTA_OUTPUT_PATH)

    input_type = classify_input_type(INPUT_KEY)


    if separate_bioproject & (input_type == "bioproject"):
        if (type(INPUT_KEY) == list):
            bioproject_id = INPUT_KEY[0]
        else:
            bioproject_id = INPUT_KEY

        FASTQ_OUTPUT_PATH += bioproject_id + "/"
        FASTA_OUTPUT_PATH += bioproject_id + "/"

        if not os.path.exists(FASTQ_OUTPUT_PATH):
            os.mkdir(FASTQ_OUTPUT_PATH)
        if not os.path.exists(FASTA_OUTPUT_PATH):
            os.mkdir(FASTA_OUTPUT_PATH)
    
    return FASTQ_OUTPUT_PATH, FASTA_OUTPUT_PATH


def fastq_to_fasta(input_file, FASTQ_OUTPUT_PATH, FASTA_OUTPUT_PATH):
    output_file = input_file.split('.')[0] + '.fasta'
    in_handle = open(FASTQ_OUTPUT_PATH+input_file, "r")
    out_handle = open(FASTA_OUTPUT_PATH+output_file, "w")

    SeqIO.convert(in_handle, "fastq", out_handle, "fasta")

    in_handle.close()
    out_handle.close()
    print(f'{output_file} has been created.')



def download_sra_data(sra_list, FASTQ_OUTPUT_PATH, FASTA_OUTPUT_PATH, blast_flag):

    # with open(FASTQ_OUTPUT_PATH+"fastq.log", "r") as f:
    #     loaded_run_ids = f.read().split('\n')

    loaded_run_ids = [filename.split('.')[0] for filename in os.listdir(FASTA_OUTPUT_PATH)]

    sra_to_run = {}

    for i, sra in enumerate(sra_list):
        
        # Get run ids from sra
        print(f'{i} start sra id: {sra} ')

        out = os.popen(f"wget -qO- 'http://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?save=efetch&db=sra&rettype=runinfo&term={sra}' | grep {sra} | cut -f1 -d','").read()
        runs = out.split('\n')[:-1]
        # if i % 2 == 0:
        #     out = os.popen(f"wget -qO- 'http://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?save=efetch&db=sra&rettype=runinfo&term={sra}' | grep {sra} | cut -f1 -d','").read()
        #     runs = out.split('\n')[:-1]
        # else:
        #     out = os.popen(f"wget -qO- 'http://www.ebi.ac.uk/ena/data/warehouse/filereport?accession={sra}&result=read_run' | grep -v accesssion | cut -f6").read()
        #     runs = out.split('\n')[1:-1]

        print(f'\tFound {len(runs)} runs: {runs}.')
        
        # Select most recent run id
        lastest_run = get_latest_run_version(runs)
        sra_to_run[sra] = lastest_run
        
        if lastest_run not in loaded_run_ids:
            # Download Fastq/Fasta from run id
            print(f'\tDownloading Fastq for {lastest_run}...')
            os.system(f"fastq-dump --outdir {FASTQ_OUTPUT_PATH} {lastest_run}")

            # log for finished downloading
            with open(FASTQ_OUTPUT_PATH+"fastq.log", "a") as f:
                f.write(lastest_run+"\n")

            fastq_to_fasta(f'{lastest_run}.fastq', FASTQ_OUTPUT_PATH, FASTA_OUTPUT_PATH)
            print('......Finished.....\n')
        else:
            print(f"\tThis SRA, {lastest_run}, is already stored locally")

        if blast_flag:
            blast_preprocess.split_and_blast(lastest_run, blast_flag)

    return sra_to_run

def create_log_file(INPUT_KEY, sra_list, sra_to_run):
    filename = datetime.now().strftime('download_log_%Y_%m_%d_%H_%M.log')
    with open(OUTPUT_PATH+filename, "w") as f:
        f.write("Input " + ", ".join(INPUT_KEY) + "\n")
        f.write("SRA list: " + ", ".join(sra_list) + "\n")
        f.write("run list: " + ", ".join(sra_to_run.values()) + "\n")

        f.write("Detail:\n")

        for sra, run in sra_to_run.items():
            f.write(f"{sra} : {run}\n")


def run(INPUT_KEY):

    INPUT_KEY, blast_flag = process_input(INPUT_KEY)

    FASTQ_OUTPUT_PATH, FASTA_OUTPUT_PATH = create_folder_and_path(INPUT_KEY)

    # get sra ids to download (input may be Bioproject, it's needed to be converted)
    sra_list = get_sra_list(INPUT_KEY)
    sra_to_run = download_sra_data(sra_list, FASTQ_OUTPUT_PATH, FASTA_OUTPUT_PATH, blast_flag)
    

    create_log_file(INPUT_KEY, sra_list, sra_to_run)

    print(", ".join(sra_to_run.values()), "are downloaded")


Entrez.email = EMAIL



if __name__ == '__main__':
    if len(argv) <= 1:
        print("Please enter input key (Bioproject ID or SRA ID(s)) after the filename")
        sys.exit()
    
    INPUT_KEY = argv[1:]

    run(INPUT_KEY)






