#!/usr/bin/env python3

import argparse
import sqlite3
import os
import gzip

def test_connection(conn): # Debugging only
    try:
        cursor = conn.cursor()
        cursor.execute("SELECT COUNT(*) FROM mgnifam")
        row_count = cursor.fetchone()[0]
        print("Connection successful! Number of rows in the table: ", row_count)
    except sqlite3.Error as e:
        print("Connection failed:", e)

def construct_file_path(pipeline_results, biome_results, domain_results, col_iter, file_name):
    file_directory = {
        0: f"{pipeline_results}/generate_families/families/seed_msa",
        1: f"{pipeline_results}/generate_families/families/hmm",
        2: f"{pipeline_results}/generate_families/families/rf",
        3: f"{pipeline_results}/structures/esmfold/cif",
        4: biome_results,
        5: domain_results,
        6: f"{pipeline_results}/annotation/reps/s4pred/json",
        7: f"{pipeline_results}/annotation/reps/deeptmhmm/json"
    }
    
    directory = file_directory.get(col_iter)
    return os.path.join(directory, file_name)

def get_blob_column(col_iter):
    blob_col_name = {
        0: "seed_msa_blob",
        1: "hmm_blob",
        2: "rf_blob",
        3: "cif_blob",
        4: "biome_blob",
        5: "domain_blob",
        6: "s4pred_blob",
        7: "tm_blob"
    }
    
    return blob_col_name.get(col_iter)

def read_file(file_path):
    try:
        with open(file_path, 'rb') as file:
            return file.read()
    except (OSError, IOError) as e:
        print(f"Error reading file {file_path}: {e}")
        return None

def update_blob_column(db, column_name, blob_data, row_id):
    try:
        conn = sqlite3.connect(db)
        cursor = conn.cursor()
        query = f"UPDATE mgnifam SET {column_name} = ? WHERE id = ?"
        cursor.execute(query, (sqlite3.Binary(blob_data), row_id))
        conn.commit()
        conn.close()
        print(f"Updated {column_name} for row {row_id}")
    except sqlite3.Error as e:
        print(f"Failed to update {column_name} for row {row_id}: {e}")

def process_row(db, pipeline_results, biome_results, domain_results, row):
    family_id = str(row[0])

    # Expected file suffixes for each column
    file_suffixes = [
        ".fas.gz", # seed_msa (to decompress)
        ".hmm",    # hmm
        ".txt",    # rf
        ".cif",    # cif
        ".csv",    # biomes
        ".json",   # domains
        ".json",   # s4pred
        ".json"    # tm
    ]

    for i, suffix in enumerate(file_suffixes):
        file_name = family_id + suffix
        file_path = construct_file_path(pipeline_results, biome_results, domain_results, i, file_name)
        blob_column = get_blob_column(i)

        # Special case: decompress seed_msa
        if i == 0:
            try:
                with gzip.open(file_path, 'rb') as f:
                    blob_data = f.read()
            except (OSError, IOError) as e:
                print(f"Error reading/decompressing seed MSA {file_path}: {e}")
                blob_data = None
        else:
            blob_data = read_file(file_path)

        if blob_data:
            update_blob_column(db, blob_column, blob_data, family_id)
        else:
            print(f"Skipping update for {blob_column} of row {family_id} due to read error")

def import_files(db, pipeline_results, biome_results, domain_results):
    conn = sqlite3.connect(db)
    cursor = conn.cursor()
    cursor.execute("SELECT id FROM mgnifam")
    rows = cursor.fetchall()
    conn.close()

    for row in rows:
        process_row(db, pipeline_results, biome_results, domain_results, row)

def main():
    parser = argparse.ArgumentParser(description="Append BLOBs to mgnifam table from pipeline result files")
    parser.add_argument("--db", help="Path to the SQLite database file")
    parser.add_argument("--pipeline_results", help="Path to the pipeline results directory")
    parser.add_argument("--biome_results", help="Path to the biome distribution CSVs")
    parser.add_argument("--domain_results", help="Path to the domain architecture JSONs")

    args = parser.parse_args()

    import_files(args.db, args.pipeline_results, args.biome_results, args.domain_results)

if __name__ == "__main__":
    main()
