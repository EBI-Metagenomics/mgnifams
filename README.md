# mgnifams

Execute these workflows in the following order:
1. preprocess_input
Starting from /nfs/production/rdf/metagenomics/users/vangelis/plp_flatfiles_pgsql_2/sequence_explorer_protein.csv.bz2
this workflow unzips the file (bzip -dk) and then removes the header, so it can be split into text chunks in the next workflow.

2. initiate proteins
There are two different options for the same file to start this workflow.
One is situated in /nfs/production/rdf/metagenomics/users/vangelis/plp_flatfiles_pgsql_2/sequence_explorer_protein_no_header.csv
and is preprocessed manually (bzip -dk and then remove header)
and the other is by running the preprocess_input workflow, which produces the same file in 
/nfs/production/rdf/metagenomics/users/vangelis/mgnifams//data/output/input/sequence_explorer_protein_no_header.csv"
The workflow then splits the input csv into 10 chunks and slices the known pfams, keeping sequences >= 50 amino acids.
