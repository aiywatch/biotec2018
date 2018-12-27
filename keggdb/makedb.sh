mkdir keggdb
cat download/aaseq/* > keggdb/aaseq.fasta
makeblastdb -in keggdb/aaseq.fasta -input_type fasta -dbtype prot -out keggdb/aaseq