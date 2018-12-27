import re
from sys import argv

# from utils import *

def count_function(BLAST_OUTPUT_PATH):
    qseqid_index = PROTEIN_BLAST_FORMAT.split().index("qseqid")
    stitle_index = PROTEIN_BLAST_FORMAT.split().index("stitle")

    function_count = {}
    with open(BLAST_OUTPUT_PATH) as f:
        for line in f:
            line_list = line.replace("\n", "").split("\t")
            qseqid = line_list[qseqid_index]
            stitle = line_list[stitle_index]
            count = int(re.findall(r"c=(\d+)", qseqid)[0])
            function = stitle.split("|")[0]
            if function in function_count:
                function_count[function] += count
            else:
                function_count[function] = count
    return function_count

def export_function_count(function_count):
    name = BLAST_OUTPUT_PATH.split('/')[0]

    with open(f"{ANALYSIS_PATH}{name}_function_count.tsv", "w") as f:
        for function, count in function_count.items():
            f.write(f"{function}\t{count}\n")
    print(f"{ANALYSIS_PATH}{name}_function_count.tsv has been saved")

def run(BLAST_OUTPUT_PATH):
	function_count = count_function(BLAST_OUTPUT_PATH)
	export_function_count(function_count)




ANALYSIS_PATH = "analysis/"

# BLAST_OUTPUT_PATH = 'out_test'
PROTEIN_BLAST_FORMAT = "qseqid stitle pident length mismatch gapopen evalue bitscore"




if __name__ == '__main__':
    if len(argv) <= 1:
        print("Please enter input key (blast output) after the filename")
        sys.exit()

    BLAST_OUTPUT_PATH = argv[1]

    run(BLAST_OUTPUT_PATH)
