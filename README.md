# MGnifams

The end-to-end pipeline is hard to execute at once, due to the family generation step which is heavy and clunky. Execute these workflows in the following order instead:

1. **preprocess_input**

The starting point of this workflow is the latest version of the sequence_explorer_protein file. This is the result CSV file of the flatfile parsing pipeline that produces the mgnify proteins db content (plp -MGnify90v2024_03). It can be in .csv.gz, .csv.bz2 or uncompressed file format (.csv) and contains data headers. The latest file version (v4) can be found at /nfs/production/rdf/metagenomics/users/vangelis/plp_flatfiles_pgsql_4/sequence_explorer_protein.csv and contains 717,738,164 (~718M) proteins.

Depending on the extension of the file, pass arguments --compress_mode 'gz' and --sequence_explorer_protein_path /nfs/production/rdf/metagenomics/users/vangelis/plp_flatfiles_pgsql_4/sequence_explorer_protein.csv.gz

The workflow decompresses the file (if compressed) and then removes the header. Run this workflow from the main mgnifams directory with the command:

**nextflow run workflows/preprocess_input/main.nf -c nextflow_workflows.config -c workflows/preprocess_input/nextflow.config -dsl2 -profile slurm**

To save storage space, we don’t store this intermediate file, but instead need to pass the resulting workDir result path in the next workflow (initiate_proteins). Example: /hps/nobackup/rdf/metagenomics/service-team/users/vangelis/work_mgnifams/7c/a25c89b1269d58b2cdfbab35cd58dc/sequence_explorer_protein_no_header.csv

2. **initiate_proteins**

After the initial input file has been optionally decompressed and had it’s header removed, we pass it to the initiate_proteins workflow. The output file of preprocess_input workflow must be passed in the initiate_proteins workflow nextflow.config.
The initiate_proteins subworkflow then in parallel chunks of 50M (default value) proteins:

(i) Slices off MGnify sequence regions with known pfam annotations

(ii) Keeps sliced sequences with AA >= min_sequence_length (default: 100)

(iii) Concatenates and exports results in output/mgnifams_input.fa which will then be clustered on the next step

This is done through a python script named bin/filter_unannotated_slices_fasta.py. For the latest version and a chunking value of 50M, there are 15 parallel chunks, and each finishes in ~20'.

Run with:

**nextflow run workflows/initiate_proteins/main.nf -c nextflow_workflows.config -dsl2 -profile slurm**

3. **execute_clustering**

For execute_clustering, we are using the latest version (15) of mmseqs tools. createdb and linclust updated in nf-core, but in MGnifams we will use ‘compressed 0’ arg instead. For createtsv we are using are own script.

linclust parameters:
sequence identity = 0.5
coverage = 0.9, both

Run with: **nextflow run workflows/execute_clustering/main.nf -c workflows/execute_clustering/nextflow.config -dsl2 -profile slurm -with-tower -resume**

4. **generate_families**

This is the main algorithm for the generation of MGnifams. This step is better run with LSF because we don't need to specify a time limit.
The create_families subworkflow calls the refine_families.py to iteratively recruit sequences starting from the largest family towards the smallest (minimum_members threshold 50).

**Restart strategy**

To continue from where the last family finished on the server execution, the files named **updated_mgnifams_dict.fa**, **updated_refined_families.tsv** and **updated_discarded_clusters.txt** must be manually copied/moved from the work folder to **/nfs/production/rdf/metagenomics/users/vangelis/mgnifams/data/output/input/families**. After removing or renaming (for backup) the three old files there, the three new ones must be renamed to **mgnifams_dict.fa**, **refined_families.tsv** and **discarded_clusters.txt** respectively. Also, make sure that the **“iteration”** parameter in nextflow.config in the generate_families subworkflow is set to the numeric value of the last properly checked family (e.g. if last generated family is mgnifam111, set the iteration param to 111). Add -resume to nextflow run.

Run with: **nextflow run workflows/generate_families/main.nf -c workflows/generate_families/nextflow.config -dsl2 -profile slurm -with-tower -resume**

5. **reformat_msa**

Results from the previous step must be manually moved to:

/nfs/production/rdf/metagenomics/users/vangelis/mgnifams/data/output/families/
in batches (all families hard to finish successfully).

Example from within the work folder:

cp msa/* /nfs/production/rdf/metagenomics/users/vangelis/mgnifams/data/output/families/msa/

cp domtblout/* /nfs/production/rdf/metagenomics/users/vangelis/mgnifams/data/output/families/domtblout/

cp seed_msa/* /nfs/production/rdf/metagenomics/users/vangelis/mgnifams/data/output/families/seed_msa/

cp hmm/* /nfs/production/rdf/metagenomics/users/vangelis/mgnifams/data/output/families/hmm/

cp updated_* /nfs/production/rdf/metagenomics/users/vangelis/mgnifams/data/output/input/families/

Reformating output **seed_msa_sto** and hmmalign **msa_sto** Stockholm MSAs in simple fas format to be able to annotate models and visualize/download with MSAViewer on the MGnifams site.

Run with: **nextflow run workflows/reformat_msa/main.nf -c workflows/reformat_msa/nextflow.config -dsl2 -profile slurm -with-tower -resume**

6. **annotate_models**

Only using **seed_msa** folder in this step.

Annotating hmm models with HHblits against HH Pfams.

-E evalue cutoff for output does not work properly. This is why we are filtering with our own code to report results.

Run with: **nextflow run workflows/annotate_models/main.nf -c workflows/annotate_models/nextflow.config -dsl2 -profile slurm -with-tower -resume**

7. **predict_annotate_structures**

Producing, PDB, CIF and CFZ files in the predict step.

foldseek results filtered by e-value 0.001: -e 0.001 (default not working, need to state).

No additional filtering carried out.

Command to extract high-quality pLDDT candidates from ESMFold log results:
grep -E 'pLDDT [7-9][0-9]\.|pLDDT 100' pdb2_scores.txt | awk -F' |, ' '{print $11, $15, $16}'

# Post workflow

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

**nextflow run main.nf -dsl2 -c nextflow_end-to-end.config -profile local -resume**

**nextflow run main.nf -dsl2 -c nextflow_end-to-end.config -profile slurm -resume -N vangelis@ebi.ac.uk**

# Anti-bus factor measures
Currently, extra documentation can be found in my google doc: https://docs.google.com/document/d/1eeglnQb9M-D0iK9AFbTypLYvvKHeUg6XtzmlKN874k4/edit
