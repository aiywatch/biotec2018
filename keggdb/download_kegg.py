import requests
import os
from multiprocessing import Pool


def download_org(organism):
    if "Prokaryotes" in organism:
        organism = organism.split('\t')
        org_id, org_code, org_name, org_tax = organism
        print(org_name)
        
        response = requests.get(f'http://rest.kegg.jp/list/{org_id}')
        f_protein_fn.write(response.text)
        proteins = response.text.split('\n')
        protein_ids = [(line.split('\t')[0]) for line in proteins if line != ""]
        protein_len = len(protein_ids)
        print(protein_len)
        with open(f"download/aaseq/aaseq_{org_id}.fasta", "w") as f_seq:
            for i in range(0, protein_len, 10):
                id_list = "+".join(protein_ids[i:i+10])
                response = requests.get(f'http://rest.kegg.jp/get/{id_list}/aaseq')
                f_seq.write(response.text)

if not os.path.exists("download"):
    os.makedirs("download")
if not os.path.exists("download/aaseq"):
    os.makedirs("download/aaseq")


response = requests.get('http://rest.kegg.jp/list/organism')

with open("download/organism.tsv", "w") as f:
    f.write(response.text)


organisms = response.text.split('\n')
f_protein_fn = open("download/kegg_protein_fn.tsv", "w")

with Pool(2) as p:
    p.map(download_org, organisms)

f_protein_fn.close()


