from Bio import SeqIO
import math
import os
from sys import argv
import sys

from config import *
from utils import *

""" 
This file takes SRR data in FastA from download_sra_data.py to blast.
It firstly splits the FastA files into chunks and blasts them paralelly.
Input: List of SRR to blast.
Output: 1. Splitted FastA files for blast in blast/ folder.
		2. run_blast script which used to run blast

Optional: If keyword "blast" is in the input command line, 
the program will call blast command automatically. (HPC required)
"""


def create_dir_and_path(one_key):
	srr = one_key.replace(".fasta", "").replace(",", "")
	output_path = f"{OUTPUT_DIR}splitted_{srr}/"


	if not os.path.exists(OUTPUT_DIR):
	    os.mkdir(OUTPUT_DIR)
	if not os.path.exists(output_path):
	    os.mkdir(output_path)

	return srr, output_path

def _build_script_file(output_path, srr, SPLIT_NUMBER):

	path = os.getcwd() + '/' + output_path
	blast_file = f"""#!/bin/bash

#BSUB -J exampleblast[1-{SPLIT_NUMBER}]
#BSUB -o {path+srr}.out.%J.%I

module load biobuilds
export BLASTDB={BLASTDB_PATH}
blastn -out {path+srr}_part${{LSB_JOBINDEX}}.out \
-outfmt '{OUTFMT}' \
-query {path+srr}_part${{LSB_JOBINDEX}}.fasta \
-db {DB_PATH} -evalue {E_VALUE} -max_target_seqs {MAX_TARGET_SEQ}
"""

	with open(output_path + "run_blast.sh", "w") as f:
		f.write(blast_file)


def split_fasta(srr, output_path):
	def _write_part_file(write_frag, i, output_path):
	    SeqIO.write(write_frag, 
	                f"{output_path + srr}_part{i}.fasta", 
	                "fasta")


	input_file = FASTA_PATH+srr
	total_records = sum(1 for record in SeqIO.parse(input_file+".fasta", "fasta"))
	records_per_file = math.ceil( total_records / SPLIT_NUMBER )
	write_frag = []

	for i, record in enumerate(SeqIO.parse(input_file+".fasta", "fasta")):
	    write_frag.append(record)
	    
	    if ( (i+1) % records_per_file == 0) or (i==(total_records-1)):
	        _write_part_file(write_frag, (i//records_per_file)+1, output_path )
	        write_frag = []

	_build_script_file(output_path, srr, SPLIT_NUMBER)

	print(f'Finished splitting {srr}')

def split_and_blast(INPUT_KEY, blast_flag):

	if type(INPUT_KEY) != list:
		INPUT_KEY = [x.strip() for x in INPUT_KEY.split(',')]

	for one_key in INPUT_KEY:
		srr, output_path = create_dir_and_path(one_key)
		split_fasta(srr, output_path)
		if blast_flag:
			os.system(f"echo 'bsub < {output_path}run_blast.sh'")
			os.system(f"bsub < {output_path}run_blast.sh")



def run(INPUT_KEY):

	INPUT_KEY, blast_flag = process_input(INPUT_KEY)
	split_and_blast(INPUT_KEY, blast_flag)



if __name__ == '__main__':
	if len(argv) <= 1:
		print("Please enter input key (SRA ID(s)) after the filename")
		sys.exit()

	INPUT_KEY = argv[1:]

	run(INPUT_KEY)