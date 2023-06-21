#!/bin/bash

# Clear subfolders to REDO
rm -rf mmseqs_dbs/*
rm -rf mmseqs_clus/*
rm pgsql_clusters.tsv

# Create a database for your initial sequences
mmseqs createdb pgsql.fasta mmseqs_dbs/pgsql_db

# mmseqs cluster mmseqs_dbs/pgsql_db mmseqs_clus/pgsql_clu mmseqs_clus/tmp --min-seq-id 0.7 --cov-mode 1 -c 0.8
# Perform fast clustering with linclust
mmseqs linclust mmseqs_dbs/pgsql_db mmseqs_clus/pgsql_clu mmseqs_clus/tmp --min-seq-id 0.7 --cov-mode 1 -c 0.8

# TODO remove for speed
# Human readable, TODO test for greater number of sequences
mmseqs createtsv mmseqs_dbs/pgsql_db mmseqs_clus/pgsql_clu pgsql_clusters.tsv

# TODO
# Create a database for your new sequences
# mmseqs createdb new_proteins.fasta mmseqs_dbs/new_db

# Merge the old and new databases
# mmseqs mergedbs mmseqs_dbs/pgsql_db combined_db mmseqs_dbs/pgsql_db mmseqs_dbs/new_db

# mmseqs clusterupdate mmseqs_dbs/pgsql_db combined_db pgsql_clu new_clu tmp --min-seq-id 0.7 --cov-mode 1 -c 0.8
