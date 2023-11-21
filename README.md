# mgnifams

Execute these workflows in the following order:

1. **preprocess_input**

Starting from */nfs/production/rdf/metagenomics/users/vangelis/plp_flatfiles_pgsql_2/sequence_explorer_protein.csv.bz2*
this workflow unzips the file (bzip -dk) and then removes the header,
so it can be split into text chunks to process in parallel in the next workflow.

Run with: **nextflow run workflows/preprocess_input/main.nf -c subworkflows/preprocess_input/nextflow.config -dsl2 -profile slurm -with-tower -resume**

2. **initiate proteins**

There are two different path options for the same input file for this workflow.
One is */nfs/production/rdf/metagenomics/users/vangelis/plp_flatfiles_pgsql_2/sequence_explorer_protein_no_header.csv*
which has been preprocessed manually (pbzip2 -dk and then removed header)
and the other is the output of the preprocess_input workflow 
*/nfs/production/rdf/metagenomics/users/vangelis/mgnifams/data/output/input/sequence_explorer_protein_no_header.csv*. 
The workflow then splits the input CSV into file chunks of 50M rows/sequences and slices out the already known pfams,
keeping sequences >= 50 amino acids in a fasta file (mgnifams_input.fa).

Run with: **nextflow run workflows/initiate_proteins/main.nf -c subworkflows/initiate_proteins/nextflow.config -dsl2 -profile slurm -with-tower -resume**
