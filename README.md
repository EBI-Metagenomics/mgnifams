# MGnifams

Starting from MGnify proteins, MGnifams aim to: i) generate microbial sequence families to expand the currently known protein family space, ii) aid in their characterisation and iii) investigate potential novel functionalities.

## Nextflow pipeline

slurm:  
nextflow run main.nf -c conf/end-to-end.config -profile slurm -with-tower  
local:  
nextflow run main.nf -c conf/end-to-end.config -profile local

![alt text](images/end-to-end.jpg)

The end-to-end MGnifams pipeline chains the subworkflows of four thematically different workflows; setup_clusters, generate_nonredundant_families, annotate_families and export_data. After the pipeline finishes, all necessary CSV files to be imported in an sqlite database (mgnifams.sqlite3) will have been produced in the output tables/ folder. Following, the bin/helper/append_blobs_sqlite.py must be run, to append all blob file items to the sqlite db. Then, the db must be copied to either the mgnifams-site repo for local testing, or directly to ifs (/nfs/public/rw/metagenomics/mgnifams/dbs) to be finally deployed online with k8s.

Another option (not recommended) is running the MGnifams pipeline’s four main workflows sequentially, each consisting of its respective same-coloured subworkflows below:

![alt text](images/workflows.png)

After all necessary CSV output tables are produced by the pipeline, do the following:  
* import CSVs into sqlite db  
* append blobs to db  
* copy db to site/ifs  
* host online with k8s  

To run workflows other than the main.nf in the root directory of the project, an environment variable must be set for the latest versions of Nextflow:

export NXF_SINGULARITY_HOME_MOUNT=true

“Changed in version 23.07.0-edge: Nextflow no longer mounts the home directory when launching an Apptainer container. To re-enable the old behavior, set the environment variable NXF_APPTAINER_HOME_MOUNT to true.“ https://nextflow.io/docs/latest/container.html

### 1. setup_clusters

slurm:  
nextflow run workflows/setup_clusters/main.nf -profile slurm -with-tower  
or local:  
nextflow run workflows/setup_clusters/main.nf -profile local  

This is the first workflow to be executed before the main family generation. It consists of three subworkflows; preprocess_input, initiate_proteins and execute_clustering. In a nutshell, this workflow converts the initial input (see below) into family-generation-ready input.

The initial input for this pipeline is the output file of the protein-landing-page data generation pipeline, sequence_explorer_protein.csv (e.g., /nfs/production/rdf/metagenomics/users/vangelis/plp_flatfiles_pgsql_4/sequence_explorer_protein.csv). In case this file is compressed, there are two different decompression modes available; gz and bz2. Set the --compress_mode parameter accordingly. Then, the known pfam domains are sliced off from proteins and we filter the remaining proteins to be above a given length threshold with the min_sequence_length parameter (e.g., >=100 AA).

### 2. generate_nonredundant_families

slurm:  
nextflow run workflows/generate_nonredundant_families/main.nf -profile slurm -with-tower -resume  
local:  
nextflow run workflows/generate_nonredundant_families/main.nf -profile local  

This workflow is the essence of MGnifams and is responsible for converting initial clusters into nonredundant protein families. The clusters from the previous workflow are chunked (minimum_members threshold=50, for clusters to keep) and then, along with the mgnifams_input.fa file they are fed into the generate_families_parallel subworkflow, which iteratively recruits sequences in the families, for each clusters’ chunk. The results are then pooled and checked for transmembrane families (flagged and removed) and redundancy among families (keeping uniques) via the flag_transmembrane and remove_redundancy subworkflows. DeepTMHMM is used to predict TM protein regions. If TM regions (alpha or beta sheet) are more or equal to 0.4 of the total fam rep length, then the family is mapped as transmembrane and discarded. We chose a more relaxed threshold to flag families rather than the strict measure DeepTMHMM uses to flag them (lowest encountered so far 0.0840). To achieve removing redundancy, an hhsuite compatible database is created, incorporating all MGnifams, and then through hhblits, an initial filtering is carried out based on e-value scores. Then for each pair of remaining families a Jaccard index score is calculated for the base MGYPs in the families. If this score is equal or more than 0.5, the per Amino Acid Jaccard Similarities is also calculated. In more detail, if the rep from the first family is found in the second family, a per AA set is constructed for each family for all same-to-rep MGYPs in the family (subslices of the same protein/slice). If the second family doesn’t contain the rep of the first, we try the reverse (second family’s rep in the first family). If again it doesn’t exist, we do nothing. Else, we calculate the AA Jaccard Index; if it is equal or above 0.95 we remove the second family as redundant, else if it is above equal or 0.5 we write out the similarity edge between the two families, else we do nothing. The remaining families are then assigned a unique integer ID.

### 3. annotate_families

slurm:  
nextflow run workflows/annotate_families/main.nf -profile slurm -with-tower -resume  
local:  
nextflow run workflows/annotate_families/main.nf -profile local  

This workflow is responsible for pulling both model and structural annotations for MGnifams. The first subworkflow, reformat_msa, is used to reformat the MSA files to be usable for the downstream subworkflows. Then, distant Pfam annotations are searched through hhsuite/hhblits for the model through the annotate_models subworkflow. In parallel, the predict_structures subworkflow predicts the family representative structures (first sequence of full msa). In some cases, some very long sequences can’t receive sufficient GPU virtual memory on the cluster to predict their structures. These will show in the pdb*_scores.txt file as: 24/05/25 22:16:10 | INFO | root | Failed (CUDA out of memory) on sequence 80 of length 1180. The EXTRACT_LONG_FA gathers these sequences and runs the prediction on the CPU via the ESMFOLD_CPU module. The results of both GPU and CPU predictions are then merged and fed in the annotate_structures subworkflow. This subworkflow is responsible for identifying structural homologs by using foldseek against the PDB, AlphaFolDB and ESM databases.

### 4. export_data

slurm:  
nextflow run workflows/export_data/main.nf -profile slurm -with-tower -resume  
local:  
nextflow run workflows/export_data/main.nf -profile local  

The final workflow, export_data, creates all the CSV tables and BLOB files with all required data and metadata for the MGnifams database. This consists of two different execution units; the first one is parsing files from the outputs of the pipeline into MGnifam CSV tables and the second one is querying the MGnify Proteins database (PGSQL) for additional post-processing information regarding underlying biomes and domain architectures of the families. The result CSV tables include; mgnifam.csv, mgnifam_proteins.csv, mgnifam_folds.csv and mgnifam_pfams.csv. The result post-processing files include two id-to-name mapping files (biomes and pfams from MGnify Proteins database), the query results for each family’s proteins for metadata against the MGnify Proteins database and finally the respective biome and domain results that are going to be appended as BLOBs in the mgnifams database, along with other families generated from previous workflows (MSAs, HMM, CIF, etc.).

A db_config.ini filepath with secrets must be set in the export_data nextflow.config.

```
[database]
dbname = ***
user = ***
password = ***
host = ***
port = ***
```

## Final steps

Manually execute the next steps to finalise setting up the MGnifams database and online website.

### Loading in sqlite

Step 1: Create the SQLite database from the schema  
sqlite3 DB/mgnifams.sqlite3 < DB/schema.sqlite

Step 2: Import data from CSV files  
Import the CSV table files into the database.  
For example, through datagrip, right click on each table and import respective file.

Step 3: Append BLOBs to db  
python bin/helper/append_blobs_sqlite.py

Step 4: Test the mgnifams-site locally  
python manage.py collectstatic --noinput
python manage.py migrate --fake
python manage.py runserver 0.0.0.0:8000

### Hosting with k8s

From within the main mgnifams-site repo:  
Update Docker image  
sudo systemctl start docker  
sudo docker build -t quay.io/microbiome-informatics/mgnifams_site:ebi-wp-k8s-hl .

Push to quay.io  
sudo docker login quay.io  
sudo docker push quay.io/microbiome-informatics/mgnifams_site:ebi-wp-k8s-hl

Move sqlite3 DB from local machine to /nfs/public/rw/metagenomics/mgnifams/dbs  
slurm:  
salloc -t 3:30:00 --mem=8G -p datamover

wormhole send mgnifams_site/dbs/mgnifams.sqlite3

This needs to be added to ~/.zshrc:  
MIT_BASERC="/hps/software/users/rdf/metagenomics/service-team/repos/mi-automation/team_environments/codon/baserc.sh"

if [ -f $MIT_BASERC ]; then  
  . $MIT_BASERC  
fi  
mitload miniconda; conda activate wormhole

wormhole receive code-id (e.g., wormhole receive 8-saturday-endorse)

chmod 775 mgnifams.sqlite3 after moving the db there

k8s:  
kubectl apply -f ebi-wp-k8s-hl.yaml

restarts:  
kubectl rollout restart deployment mgnifams-site

## Anti bus-factor 1 measures

Currently, extra documentation can be found in my google doc: https://docs.google.com/document/d/1eeglnQb9M-D0iK9AFbTypLYvvKHeUg6XtzmlKN874k4/edit
