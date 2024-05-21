# MGnifams

Starting from MGnify proteins, MGnifams aims to generate microbial sequence families to expand the currently known protein family space and to investigate potential novel functionalities.

## Nextflow pipeline

slurm:  
nextflow run main.nf -c conf/end-to-end.config -profile slurm -with-tower  
local:  
nextflow run main.nf -c conf/end-to-end.config -profile local

![alt text](images/end-to-end.jpg)

The end-to-end MGnifams pipeline chains the subworkflows of three thematically different workflows; setup_clusters, generate_families and annotate_families. Running the whole end-to-end pipeline in one step is not realistic with the current number of input sequences, due to the lengthy execution of the generate_families subworkflow. This would find use for estimated input sets of 1-10M sequences instead of 700M. After the pipeline finishes, the export_data workflow must be executed to produce all necessary CSV files to be then imported in an sqlite database (mgnifams.sqlite3). Following, the bin/helper/append_blobs_sqlite.py must be run, to append all blob file items to the db. Then, the db must be copied to either the mgnifams-site repo for local testing, or directly to ifs (/nfs/public/rw/metagenomics/mgnifams/dbs) to be finally deployed online with k8s.

A more realistic scenario is breaking the MGnifams pipeline into four main workflows (each consisting of its respective same-coloured subworkflows below) and executing them one after the other as shown below:

![alt text](images/workflows.png)

After the export_data produces all necessary output tables, do the following:  
* import CSVs into sqlite db  
* append blobs to db  
* copy db to site/ifs  
* host online with k8s  

### 1. setup_clusters
slurm:  
nextflow run workflows/setup_clusters/main.nf -profile slurm -with-tower  
or local:  
nextflow run workflows/setup_clusters/main.nf -profile local  

This is the first workflow to be executed before the main family generation. It consists of three subworkflows; preprocess_input, initiate_proteins and execute_clustering. In a nutshell, this workflow converts the initial input (see below) into family-generation-ready input.

The initial input for this pipeline is the output file of the protein-landing-page data generation pipeline, sequence_explorer_protein.csv (e.g., /nfs/production/rdf/metagenomics/users/vangelis/plp_flatfiles_pgsql_4/sequence_explorer_protein.csv). In case this file is compressed, there are two different decompression modes available; gz and bz2. Set the --compress_mode parameter accordingly. Then, the known pfam domains are sliced off from proteins and we filter the remaining proteins to be above a given length threshold with the min_sequence_length parameter (e.g., >=100 AA). 


### 2. generate_families

### 3. annotate_families

### 4. export_dta

### Post workflows

**Exporting database CSVs**

After all subworkflows have finished executing run:

**python bin/export/export_mgnifams_csv.py <mgnifams_out_dir>**

to generate ready to import data tables

**Querying biome and pfam data for site**

The next two scripts query the protein PGSQL db (plp) to gather and parse information regarding the family underlying protein biomes and pfams.
The generated data are used to build biome sunburst, and domain architecture plots respectively on the mgnifams-site.

A db_config.ini file with secrets must be passed in the scripts.

```
[database]
dbname = ***
user = ***
password = ***
host = ***
port = ***
```

**python3 bin/post-processing/query_biome_csvs.py bin/db_config.ini /home/vangelis/Desktop/Projects/mgnifams-site-data_backup/families/updated_refined_families.tsv 0**

**python3 bin/post-processing/query_pfam_jsons.py bin/db_config.ini /home/vangelis/Desktop/Projects/mgnifams-site-data_backup/families/updated_refined_families.tsv 0**

**Translating seed_msa and msa proteins to MGYPS**

**python bin/post-processing/translate_msa_mgyps.py /home/vangelis/Desktop/Projects/mgnifams-site-data_backup/families/seed_msa**

**python bin/post-processing/translate_msa_mgyps.py /home/vangelis/Desktop/Projects/mgnifams-site-data_backup/families/msa**

# End-to-end pipeline

This is mainly used locally for testing. Chains all aforementioned modules.

**nextflow run main.nf -c conf/end-to-end.config -profile local -resume**

**nextflow run main.nf -c conf/end-to-end.config -profile slurm -resume -N vangelis@ebi.ac.uk**

# Anti bus-factor 1 measures
Currently, extra documentation can be found in my google doc: https://docs.google.com/document/d/1eeglnQb9M-D0iK9AFbTypLYvvKHeUg6XtzmlKN874k4/edit
