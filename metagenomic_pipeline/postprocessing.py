import pandas as pd
import numpy as np
from Bio import Entrez
import xml.etree.ElementTree as ET
import csv
import re
import glob
from urllib.error import HTTPError
from concurrent import futures
import os
from multiprocessing import Pool, Manager, cpu_count
from sys import argv
import sys

from combine_result import combine_all_results
from config import *
from utils import *


"""
This file is used to do postprocessing and analyzing on blast output.
Input: Bioproject ID, SRA ID, or List of SRA IDs. These have to be already blasted
        from the previous steps
Output: 
    1. blast/clean/<SRR_ID>.out is a cleaned blast output
    2. analysis/<SRR_ID>/ this folder includes gene count, 
    species count, protein sequnces, unprocess results
    3. cache/ keeps cache when downloading data from NCBI
"""


### PART 1: PREPROCESSING STEP - CREATE FOLDER, BUILD CLEANED BLAST OUTPUT FILE

def bioproject_to_run(bioproject_list):
    sra_list = get_sra_list(INPUT_KEY)
    run_list = []

    for i, sra in enumerate(sra_list):
        
        # Get run ids from sra
        print(f'{i} start sra id: {sra} ')

        out = os.popen(f"wget -qO- 'http://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?save=efetch&db=sra&rettype=runinfo&term={sra}' | grep {sra} | cut -f1 -d','").read()
        runs = out.split('\n')[:-1]
        print(f'\tFound {len(runs)} runs: {runs}.')
        
        # Select most recent run id
        lastest_run = get_latest_run_version(runs)
        run_list += [lastest_run]
        # print(runs)

    return run_list

def reduce_blast_output_sorted(RUN_NUMBER, BLAST_OUTPUT_PATH, BLAST_CLEAN_PATH, TOP_N=1):
    """ This function is used to build a cleaned blast output file (select Top N),
        remove unneccessary features, and sort them by gb """
    blast_map = BLAST_FORMAT.split()


    if not os.path.isdir(BLAST_CLEAN_PATH):
        os.mkdir(BLAST_CLEAN_PATH)

    blast_output_floder = "/".join(BLAST_OUTPUT_PATH.split('/')[:-1])

    NUMBER_OF_FILES = len(glob.glob1(blast_output_floder, f"{RUN_NUMBER}_part*.out"))

    last_run_id = ""
    seqs = []
    for i in range(1, NUMBER_OF_FILES+1):
        with open(f"{BLAST_OUTPUT_PATH}{i}.out") as file:
            for line in file:
                line_list = line.split('\t')
                run_id = line_list[0]
                if (last_run_id != run_id):
                    count = 0
                if (count < TOP_N):
                    # new_line = "\t".join(line_list[0:4] + line_list[10:12])
                    qseqid = line_list[blast_map.index("qseqid")]
                    sseqid = line_list[blast_map.index("sseqid")]
                    sscinames = line_list[blast_map.index("sscinames")]
                    pident = line_list[blast_map.index("pident")]
                    sstart = line_list[blast_map.index("sstart")]
                    send = line_list[blast_map.index("send")]

                    clean_row = "\t".join([qseqid, sseqid, sscinames, pident, sstart, send])
                    seqs.append(clean_row)
                    count += 1
                last_run_id = run_id
    sorted_seqs = sorted(seqs, key=lambda row: row.split('\t')[1].split('|')[3])

    with open(f"{BLAST_CLEAN_PATH+RUN_NUMBER}.out", "w") as output_file:
        output_file.write("\n".join(sorted_seqs))
    print("cleaned blast output file is created at", f"{BLAST_CLEAN_PATH+RUN_NUMBER}.out")

def create_folder_and_clean_blast_output(OUTPUT_FOLDER, PROTEIN_OUTPUT_PATH, BLAST_OUTPUT_PATH, BLAST_CLEAN_PATH, RUN_NUMBER):
    ## Creating neccessary folders
    if not os.path.isdir("cache"):
        os.mkdir("cache")
        os.system("touch cache/run_id.txt")

    for folder_path in [OUTPUT_FOLDER, PROTEIN_OUTPUT_PATH]:
        accumulated_path = ""
        for folder in folder_path.split('/'):
            accumulated_path += folder + "/"
            if not os.path.isdir(accumulated_path):
                os.mkdir(accumulated_path)

    ## Create cleaned blast output

    if not os.path.isfile(f"{BLAST_CLEAN_PATH+RUN_NUMBER}.out"):
        reduce_blast_output_sorted(RUN_NUMBER, BLAST_OUTPUT_PATH, BLAST_CLEAN_PATH)


### PART 2: DOWNLOAD SPECIES INFO FROM NCBI STEP
def get_unique_gbs_from_blast_output(input_file):
    gb_to_gi = {}
    gbs = set()

    with open(input_file) as csvfile:
        reader = csv.reader(csvfile, delimiter='\t')
        for i, row in enumerate(reader):
            gi = row[1].split('|')[1]
            gb = row[1].split('|')[3]
            gbs.add(gb)
            gb_to_gi[gb] = gi

    print(f"Got {len(gbs)} unique gbs from {i+1} sequences")
    return gbs, gb_to_gi, i+1

def get_cache_gbs():
    with open("cache/run_id.txt", "r") as f:
        cache_gbs = f.read().split('\n')
    return cache_gbs

def get_gbs_storage_info(gbs):
    cache_gbs = get_cache_gbs()
    gbs_from_local = gbs.intersection(cache_gbs)
    gbs_from_external = gbs.difference(cache_gbs)
    print(f"local:{len(gbs_from_local)}  external:{len(gbs_from_external)}")
    
    return gbs_from_local, gbs_from_external


def __extract_partial_feature(gb, specie, gbfeature):
    row = []
    feature_dict = {}

    for gbqualifier in gbfeature.iter('GBQualifier'):
        
        if gbqualifier.find("GBQualifier_value") is not None:
            name = gbqualifier.find("GBQualifier_name").text
            value = gbqualifier.find("GBQualifier_value").text
            feature_dict[name] = value
    
    feature_dict['GBFeature_key'] = gbfeature.find('GBFeature_key').text
    if ('product' in feature_dict) or (feature_dict['GBFeature_key'] == 'misc_feature'):
        
        location_text =  gbfeature.find('GBFeature_location').text.replace("<", "").replace(">", "")
        extracted_locations = re.findall(r'(\d+)\.\.(\d+)',location_text)
        locations = np.reshape(extracted_locations, (-1, 2) ).astype(np.int).tolist()

        row = [gb, specie, locations, feature_dict['GBFeature_key']]
        if 'translation' in feature_dict:
            row.extend([feature_dict['product'],feature_dict['protein_id'],feature_dict['translation']])
        elif 'product' in feature_dict:
            row.append(feature_dict['product'])
        elif (feature_dict['GBFeature_key'] == 'misc_feature') & ('note' in feature_dict):
            row.append(feature_dict['note'])

    return row

def __extract_single_feature(gbseq):
    seq_features = []
    
    gb = gbseq.find('GBSeq_accession-version').text
    specie = gbseq.find("GBSeq_organism").text

    for gbfeature in gbseq.iter('GBFeature'):
        row = __extract_partial_feature(gb, specie, gbfeature)
        if len(row) > 4:
            seq_features.append(row)
    if len(seq_features) == 0:
        seq_features.append([gb, specie])
    
    return seq_features

def _extract_features(root):
    # Extract information from XML
    seq_features = []

    for gbseq in root.iter("GBSeq"):
        partial_seq_features = __extract_single_feature(gbseq)
        seq_features = seq_features + partial_seq_features

    
    return seq_features

def _cache(partial_seq_features, batch_gbs):
    with open("cache/seq_features.txt", "a") as f:
        for seq in partial_seq_features:
            f.write(str(seq) + "\n")
    with open("cache/run_id.txt", "a") as f:
        for gb in batch_gbs:
            f.write(str(gb) + "\n")

def load_external_feature_from_gbs(gbs, gb_to_gi):
    seq_features = []
    error = []

    gbs = list(gbs)
    EPOST_BATCH_SIZE = 100
    for start_index in range(0, len(gbs), EPOST_BATCH_SIZE):
#     for start_index in range(0, EPOST_BATCH_SIZE, EPOST_BATCH_SIZE):
        print(f'index: {start_index}-{start_index+EPOST_BATCH_SIZE-1}')
        
        # Fetch data from NCBI
        batch_gbs = gbs[start_index:start_index+EPOST_BATCH_SIZE]
        batch_gis = [gb_to_gi[gb] for gb in batch_gbs]
#         print(batch_gis)
        try:
            search_results = Entrez.read(
                Entrez.epost("nucleotide", id=",".join(batch_gis)))

            webenv = search_results["WebEnv"]
            query_key = search_results["QueryKey"]

            handle = Entrez.efetch(db="nucleotide", retmode="xml", 
                           retmax=EPOST_BATCH_SIZE, webenv=webenv, query_key=query_key,)
            root = ET.fromstring(handle.read())
        except (AttributeError, HTTPError, Exception) as e:
            print('Error on ', e)
            error += batch_gbs
            continue

        # Extract information from XML
        print('extracting ...')
        partial_seq_features = _extract_features(root)
        seq_features = seq_features + partial_seq_features

        _cache(partial_seq_features, batch_gbs)
    return error
    # return seq_features, error

def download_species_feature_from_ncbi(gbs, gb_to_gi):
    gbs_from_local, gbs_from_external = get_gbs_storage_info(gbs)

    while (len(gbs_from_external) != 0):
        # seq_features_external, error = load_external_feature_from_gbs(gbs_from_external, gb_to_gi)
        error = load_external_feature_from_gbs(gbs_from_external, gb_to_gi)
        gbs_from_local, gbs_from_external = get_gbs_storage_info(gbs)
        print(f"round finished with: {error}")

def load_local_features(gbs_from_local):
    seq_features_local = []
    gbs_with_protein_local = set()
    with open("cache/seq_features.txt", "r") as f:
        for line in f:
            row = eval(line)
            if row[0] in gbs_from_local:
                seq_features_local.append(row)

    global COLUMN_NAMES
    COLUMN_NAMES = ["gb", "species", "locations", "type", "gene", "protein_id", "protein"]
    sorted_matched_gbs = pd.DataFrame(seq_features_local, columns=COLUMN_NAMES).sort_values("gb")
#     sorted_matched_gbs_ids = sorted_matched_gbs["gb"].values

    return sorted_matched_gbs



### PART 3: FILTER SEQUENCES
def _is_in_location(query_location, subject_location):
    start_query, end_query = min(query_location), max(query_location)
    start_subject, end_subject = min(subject_location), max(subject_location)
    
    return (start_query in range(start_subject, end_subject+1)) and (end_query in range(start_subject, end_subject+1))

def _is_in_locations(query_location, subject_locations):
    for subject_location in subject_locations:
        if _is_in_location(query_location, subject_location):
            return True
    return False

def get_next_df(sorted_matched_gbs, last_index, sorted_matched_gbs_ids):
    first_index = last_index
    first_gb = sorted_matched_gbs_ids[first_index]
    df_len = sorted_matched_gbs.shape[0]

    while (last_index<df_len) and (first_gb == sorted_matched_gbs_ids[last_index]):
        last_index += 1
    return sorted_matched_gbs[first_index: last_index], last_index, first_gb

def get_df_from_gb(query_gb, sorted_matched_gbs, last_index, last_query_gb, sorted_matched_gbs_ids):
    
    df = pd.DataFrame([])
    gb = last_query_gb
    
    if last_query_gb == "":
        df, last_index, gb = get_next_df(sorted_matched_gbs, last_index, sorted_matched_gbs_ids)

    while (query_gb > gb):
        df, last_index, gb = get_next_df(sorted_matched_gbs, last_index, sorted_matched_gbs_ids)
    
    return df, last_index, gb

def _get_location_length(locations):
    length = 0
    for location in locations:
        min_loc, max_loc = min(location), max(location)
        length += max_loc - min_loc
    return length

def get_seq_info(line_number):

    line_range = range(line_number[0], line_number[1]+1)
    
    with open(input_file) as csvfile:
        reader = csv.reader(csvfile, delimiter='\t')
        
        last_query_gb = ""
        last_index = 0
        for i, row in enumerate(reader):
            if i not in line_range:
                continue
            
            if i % 5000 == 0:
                print(i, len(seqs))
            gb = row[1].split('|')[3]
            
            if last_query_gb != gb:
                df, last_index, last_query_gb = get_df_from_gb(gb, sorted_matched_gbs, last_index, last_query_gb, sorted_matched_gbs_ids)
                if df.shape[0] > 0:
                    found_seqs = df
            
            if (last_query_gb == gb) and (found_seqs.shape[0] > 1):

                query_location = (int(row[4]), int(row[5]))
                row_indexes = found_seqs['locations'].apply(lambda subject_locations: _is_in_locations(query_location, subject_locations))
                
                result_count = row_indexes.sum()
                
                if result_count > 1:
                    seq_opts = found_seqs[row_indexes]
                    seq_no_misc = seq_opts[seq_opts["type"].apply(lambda seq_type: "misc" not in seq_type)]
                    
                    
                    if seq_no_misc.shape[0] == 1:
                        seqs.extend(seq_no_misc.values)
                    else:
                        seq_with_protein = seq_no_misc[seq_no_misc["protein"].notnull()]
                        if seq_with_protein.shape[0] == 1:
                            seqs.extend(seq_with_protein.values)
                        elif seq_with_protein.shape[0] > 1:
                            longest_protein = seq_with_protein.loc[seq_with_protein["protein"].apply(lambda protein: len(protein)).idxmax()]
                            seqs.append(longest_protein.values)
                        else:
                            longest_seq = seq_opts.loc[seq_opts["locations"].apply(lambda locations: _get_location_length(locations)).idxmax()]
                            seqs.append(longest_seq.values)
                elif result_count == 1:
                    seqs.extend(found_seqs[row_indexes].values)
                elif result_count == 0:
                    seqs.append(found_seqs[["gb", "species"]].iloc[0].values)
            elif (last_query_gb == gb) and (found_seqs.shape[0] == 1):
                seqs.extend(found_seqs.values)
            else:
                print("error: not found this gb: ", i, gb, last_query_gb)
    return seqs

def split_line_by_cores(total_lines, number_of_cores):
    line_per_core = total_lines//number_of_cores
    file_indices = []
    for i in range(number_of_cores):
        start_index = i*line_per_core
        if i == (number_of_cores-1):
            end_index = total_lines
        else:
            end_index = (i+1)*line_per_core - 1
        file_indices.append((start_index, end_index))
    return file_indices

def filter_species_data(sorted_matched_gbs, total_lines, number_of_cores=40):
    ## filter species data by location
    global seqs, sorted_matched_gbs_ids
    
    seqs = Manager().list()
    sorted_matched_gbs_ids = sorted_matched_gbs["gb"].values
    file_indices = split_line_by_cores(total_lines, number_of_cores)

    with Pool(number_of_cores) as p:
        seq = p.map(get_seq_info, file_indices)

    seqs_df = pd.DataFrame(list(seqs), columns=COLUMN_NAMES)

    return seqs_df

### PART 4: EXPORT DATA
def export_data(seqs_df, OUTPUT_FOLDER, RUN_NUMBER, PROTEIN_BATCH_SIZE=10000):
    ## count species and genes
    species_count = seqs_df['species'].value_counts()
    species_count.to_csv(f"{OUTPUT_FOLDER + RUN_NUMBER}_species.tsv", sep='\t')
    gene_count = seqs_df['gene'].value_counts()
    gene_count.to_csv(f"{OUTPUT_FOLDER + RUN_NUMBER}_gene.tsv", sep='\t')
    seqs_df.to_csv(f'{OUTPUT_FOLDER + RUN_NUMBER}.csv')

    
    ## export protein sequences
    protein_seq_df = seqs_df[seqs_df['protein'].notnull()][["protein_id", "protein", "gene"]]
    count_unique = protein_seq_df.groupby(["protein_id", "protein"]).count().reset_index()
    if (count_unique.shape[0]>0):
        count_unique.to_csv(f"{PROTEIN_OUTPUT_PATH}/protein_seq.csv")
    print(f"found protein {protein_seq_df.shape[0]} seqs of which {count_unique.shape[0]} seqs are unique")

    if (count_unique.shape[0]>0) and ((PROTEIN_BATCH_SIZE == "all") or (PROTEIN_BATCH_SIZE == None)):
        with open(f"{PROTEIN_OUTPUT_PATH}/protein_seq.fasta", "w") as f:
            for index, row in count_unique.iterrows():
                f.write(f">{row['protein_id']}-c={row['gene']}\n")
                f.write(f"{row['protein']}\n")
    else:
        f = None
        for index, row in count_unique.iterrows():
            if (index % PROTEIN_BATCH_SIZE) == 0:
                batch_number = (index // PROTEIN_BATCH_SIZE)+1
                if f:
                    f.close()
                f = open(f"{PROTEIN_OUTPUT_PATH}/protein_seq_part{batch_number}.fasta", "w")
            # print(index, row)
            # print("aaa", row['index'])
            f.write(f">{row['protein_id']} count={row['gene']}\n")
            f.write(f"{row['protein']}\n")
        if f:
            f.close()


## PARAMETER
# TOP_N = 1
# PROTEIN_BATCH_SIZE = 10000
# NUMBER_OF_CORES = 40

# RUN_NUMBERS = "SRR3501889"
# RUN_NUMBERS = "SRR2579826"



if __name__ == '__main__':
    if len(argv) <= 1:
        print("Please enter input key (Bioproject ID or SRA ID(s)) after the filename")
        sys.exit()

    INPUT_KEY = argv[1:]
        
    INPUT_KEY, _ = process_input(INPUT_KEY)
    input_type = classify_input_type(INPUT_KEY)

    if input_type == "bioproject_list":
        RUN_NUMBERS = bioproject_to_run(INPUT_KEY)
    else:
        RUN_NUMBERS = INPUT_KEY
        
    Entrez.email = EMAIL


    for RUN_NUMBER in RUN_NUMBERS:
        print(f"start processing on {RUN_NUMBER}")

        BLAST_OUTPUT_PATH = f'blast/splitted_{RUN_NUMBER}/{RUN_NUMBER}_part'
        BLAST_CLEAN_PATH = f"blast/clean/"
        OUTPUT_FOLDER = f"analysis/{RUN_NUMBER}/"
        PROTEIN_OUTPUT_PATH = f"{OUTPUT_FOLDER}protein/"


        input_file = f'{BLAST_CLEAN_PATH+RUN_NUMBER}.out'
        create_folder_and_clean_blast_output(OUTPUT_FOLDER, PROTEIN_OUTPUT_PATH, BLAST_OUTPUT_PATH, BLAST_CLEAN_PATH, RUN_NUMBER)

        ## fetching species data
        gbs, gb_to_gi, total_lines = get_unique_gbs_from_blast_output(input_file)
        download_species_feature_from_ncbi(gbs, gb_to_gi)


        NUMBER_OF_CORES = cpu_count()-2
        sorted_matched_gbs = load_local_features(gbs)
        seqs_df = filter_species_data(sorted_matched_gbs, total_lines, NUMBER_OF_CORES)
        export_data(seqs_df, OUTPUT_FOLDER, RUN_NUMBER, PROTEIN_BATCH_SIZE)

    combine_all_results(RUN_NUMBERS, INPUT_KEY, input_type)

