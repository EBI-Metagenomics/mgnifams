#!/usr/bin/env python3

import argparse
import configparser
import os
import psycopg2
from concurrent.futures import ThreadPoolExecutor, as_completed

BATCH_SIZE = 50000

def read_config(mgnprotein_db_config_file):
    config = configparser.ConfigParser()
    config.read(mgnprotein_db_config_file)
    return dict(config.items('database'))

def extract_sequence_id(sequence_id_with_region):
    parts = sequence_id_with_region.split('/')
    if len(parts) > 1:
        return parts[0].split('_')[0]
    elif "_" in sequence_id_with_region:
        return sequence_id_with_region.split('_')[0]

    return sequence_id_with_region

def load_family_proteins(family_proteins_file):
    families = {}
    with open(family_proteins_file, 'r') as file:
        for line in file:
            family_id, sequence_id_with_region = line.strip().split('\t')
            sequence_id = extract_sequence_id(sequence_id_with_region)
            if family_id not in families:
                families[family_id] = []
            families[family_id].append(sequence_id)
    return families

def fetch_protein_metadata(db_params, all_sequences):
    """Bulk-fetch metadata for all sequences in batches. Returns {mgyp: metadata}."""
    protein_metadata = {}
    conn = psycopg2.connect(**db_params)
    try:
        cursor = conn.cursor()
        sequence_list = list(all_sequences)
        for i in range(0, len(sequence_list), BATCH_SIZE):
            batch = sequence_list[i:i + BATCH_SIZE]
            cursor.execute(
                "SELECT mgyp, metadata FROM sequence_explorer_protein WHERE mgyp = ANY(%s::bigint[])",
                (batch,)
            )
            for mgyp, metadata in cursor.fetchall():
                protein_metadata[str(mgyp)] = metadata
        cursor.close()
    finally:
        conn.close()
    return protein_metadata

def write_family_file(family_id, sequences, protein_metadata, output_dir):
    output_tsv = f"{output_dir}/{family_id}.tsv"
    with open(output_tsv, 'w') as file:
        for seq in set(sequences):
            if seq in protein_metadata:
                metadata = protein_metadata[seq]
                meta_b = metadata.get('b', '')
                meta_p = metadata.get('p', '')
                file.write(f"{seq}\t{meta_b}\t{meta_p}\n")

def query_sequence_explorer_protein(db_params, family_proteins_file, output_dir, threads):
    families = load_family_proteins(family_proteins_file)

    # Deduplicate across all families so each protein is fetched exactly once
    all_sequences = {seq for seqs in families.values() for seq in seqs}
    protein_metadata = fetch_protein_metadata(db_params, all_sequences)

    # Fan out to per-family files in parallel
    with ThreadPoolExecutor(max_workers=threads) as executor:
        futures = {
            executor.submit(write_family_file, family_id, sequences, protein_metadata, output_dir): family_id
            for family_id, sequences in families.items()
        }
        for future in as_completed(futures):
            future.result()  # Re-raise any exceptions

def query_sequence_explorer_biome(cursor):
    sql_query = "SELECT id, name FROM sequence_explorer_biome"
    cursor.execute(sql_query)
    rows = cursor.fetchall()

    output_tsv = "biome_mapping.tsv"
    with open(output_tsv, 'w') as file:
        for row in rows:
            file.write(f"{row[0]}\t{row[1]}\n")

def query_sequence_explorer_pfam(cursor):
    sql_query = "SELECT id, name FROM sequence_explorer_pfam"
    cursor.execute(sql_query)
    rows = cursor.fetchall()

    output_tsv = "pfam_mapping.tsv"
    with open(output_tsv, 'w') as file:
        for row in rows:
            file.write(f"{row[0]}\t{row[1]}\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Query the PostgreSQL database MGnifams proteins data.")
    parser.add_argument("--mgnprotein_db_config_file", help="Path to the configuration file for the database secrets")
    parser.add_argument("--family_proteins_file", help="Path to the tsv file with families and respective proteins")
    parser.add_argument("--output_dir", help="Output directory with family TSV files containing per MGnify sequence biome and pfam annotations")
    parser.add_argument("--threads", type=int, default=1, help="Number of parallel threads for writing family TSV files")

    args = parser.parse_args()

    db_params = read_config(args.mgnprotein_db_config_file)

    if not os.path.exists(args.output_dir):
        os.mkdir(args.output_dir)

    query_sequence_explorer_protein(db_params, args.family_proteins_file, args.output_dir, args.threads)

    conn   = psycopg2.connect(**db_params)
    cursor = conn.cursor()
    query_sequence_explorer_biome(cursor)
    query_sequence_explorer_pfam(cursor)
    cursor.close()
    conn.close()
