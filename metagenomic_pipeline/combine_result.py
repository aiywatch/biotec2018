from sys import argv
from datetime import datetime
import os.path
import pandas as pd
import shutil
import os
import sys

from config import *
from utils import *

"""
This file is run automatically after postprocessing.py
This file is used to combine results from several results of bioproject or list of sras.
Input: List of SRR
Output: combined result in analysis/ folder : 
	1. <PROJECT_ID>_protein : all protein sequences in the project
	2. summary_gene_<PROJECT_ID> : gene count in the project
	3. summary_species_<PROJECT_ID> : spcies count in the project
"""


def combine_result(RUN_NUMBERS, type):
	result_dict = {}

	for RUN_NUMBER in RUN_NUMBERS:
	    run_folder = f"{ANALYSIS_PATH}{RUN_NUMBER}/"
	    with open(f"{run_folder+RUN_NUMBER}_{type}.tsv") as f:
	        for line in f:
	            gene, count = line.replace("\n", "").split('\t')
	            if gene in result_dict:
	                result_dict[gene] += int(count)
	            else:
	                result_dict[gene] = int(count)

	return result_dict

def save_result_to_file(RUN_NUMBERS, result_dict, result_type, INPUT_KEY, input_type):

	if input_type == "bioproject_list":
		filename = datetime.now().strftime(f'summary_{result_type}_{INPUT_KEY[0]}')
	else:
		filename = datetime.now().strftime(f'summary_{result_type}_%Y_%m_%d_%H_%M')

	with open(ANALYSIS_PATH + filename, 'w') as f:
		if input_type == "bioproject_list":
			f.write(f"Input bioproject : {INPUT_KEY[0]}\n\n")

		srrs_str = ", ".join(RUN_NUMBERS)
		f.write(f"Input SRRs : {srrs_str}\n\n")
		for key, value in sorted(result_dict.items(), key=lambda x: x[1], reverse=True):
			f.write(f"{key}\t{value}\n")
	print(f"{filename} has been created!")

def combine_proteins(RUN_NUMBERS, proj_id):
	print("Start combining protein sequences...")

	total_prot = pd.DataFrame()
	for run_id in RUN_NUMBERS:
	    prot_f = f"{ANALYSIS_PATH}{run_id}/protein/protein_seq.csv"
	    if os.path.isfile(prot_f):
	        prot_df = pd.read_csv(prot_f, index_col=0)
	        total_prot = total_prot.append(prot_df)

	count_unique = total_prot.groupby(["protein_id", "protein"]).count().reset_index()

	if count_unique.shape[0]>0:
		with open(ANALYSIS_PATH+proj_id+"_protein.fasta", "w") as f:
			for index, row in count_unique.iterrows():
				f.write(f">{row['protein_id']}-c={row['gene']}\n")
				f.write(f"{row['protein']}\n")
	            
	print("Combining protein sequences finished!")

def move_folder(RUN_NUMBERS, proj_id):
	""" move all folders and files from previous steps to BIOPROJECT folder """
	proj_dir = ANALYSIS_PATH+proj_id
	if not os.path.exists(proj_dir):
		os.makedirs(proj_dir)

	for run_id in RUN_NUMBERS:
		source_dir = ANALYSIS_PATH+run_id
		if os.path.exists(source_dir):
			shutil.move(source_dir, proj_dir)

	if os.path.exists(f"{ANALYSIS_PATH+proj_id}_protein.fasta"):
		shutil.move(f"{ANALYSIS_PATH+proj_id}_protein.fasta", proj_dir)
	if os.path.exists(f"{ANALYSIS_PATH}summary_gene_{proj_id}"):
		shutil.move(f"{ANALYSIS_PATH}summary_gene_{proj_id}", proj_dir)
	if os.path.exists(f"{ANALYSIS_PATH}summary_species_{proj_id}"):
		shutil.move(f"{ANALYSIS_PATH}summary_species_{proj_id}", proj_dir)


def combine_all_results(RUN_NUMBERS, INPUT_KEY, input_type):
	result_types = ["species", "gene"]
	for result_type in result_types:
		result_dict = combine_result(RUN_NUMBERS, result_type)
		save_result_to_file(RUN_NUMBERS, result_dict, result_type, INPUT_KEY, input_type)

	if input_type == "bioproject_list":
		combine_proteins(RUN_NUMBERS, INPUT_KEY[0])
		move_folder(RUN_NUMBERS, INPUT_KEY[0])



# RUN_NUMBERS = ["SRR2579826", "SRR2592314"]

if __name__ == '__main__':
	if len(argv) <= 1:
		print("Please enter input key (SRA ID(s)) after the filename")
		sys.exit()

	INPUT_KEY = argv[1:]

	RUN_NUMBERS, blast_flag = process_input(INPUT_KEY)

	combine_all_results(RUN_NUMBERS)



