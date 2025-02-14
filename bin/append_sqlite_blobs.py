#!/usr/bin/env python3

import sqlite3
import sys
import os

def test_connection(conn): 
    try:
        cursor = conn.cursor()
        cursor.execute("SELECT COUNT(*) FROM mgnifam")
        row_count = cursor.fetchone()[0]
        print("Connection successful! Number of rows in the table: ", row_count)
    except sqlite3.Error as e:
        print("Connection failed:", e)

def construct_file_path(pipeline_results_dir, family_dir, col_iter, file_column):
    file_directory = {
        0: "structures/cif",
        1: f"{family_dir}/seed_msa",
        2: f"{family_dir}/msa",
        3: f"{family_dir}/hmm",
        4: f"{family_dir}/rf",
        5: "post-processing/biome_results",
        6: "post-processing/domain_results"
    }
    
    directory = file_directory.get(col_iter)
    return os.path.join(pipeline_results_dir, directory, file_column)

def get_blob_column(col_iter):
    blob_col_name = {
        0: "cif_blob",
        1: "seed_msa_blob",
        2: "msa_blob",
        3: "hmm_blob",
        4: "rf_blob",
        5: "biomes_blob",
        6: "domain_architecture_blob"
    }
    
    return blob_col_name.get(col_iter)

def read_file(file_path):
    try:
        with open(file_path, 'rb') as file:
            return file.read()
    except (OSError, IOError) as e:
        print(f"Error reading file {file_path}: {e}")
        return None

def update_blob_column(db_path, column_name, blob_data, row_id):
    try:
        conn = sqlite3.connect(db_path)
        cursor = conn.cursor()
        query = f"UPDATE mgnifam SET {column_name} = ? WHERE id = ?"
        cursor.execute(query, (sqlite3.Binary(blob_data), row_id))
        conn.commit()
        conn.close()
        print(f"Updated {column_name} for row {row_id}")
    except sqlite3.Error as e:
        print(f"Failed to update {column_name} for row {row_id}: {e}")

def process_row(db_path, pipeline_results_dir, family_dir, row):
    row_id = row[0]
    file_columns = row[1:]

    for i, file_column in enumerate(file_columns):
        if file_column is not None: 
            file_path = construct_file_path(pipeline_results_dir, family_dir, i, file_column)
            blob_column = get_blob_column(i)
            blob_data = read_file(file_path)
            
            if blob_data:
                update_blob_column(db_path, blob_column, blob_data, row_id)
            else:
                print(f"Skipping update for row {row_id} due to read error")

def import_files(db_path, pipeline_results_dir, family_dir):
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()
    cursor.execute("SELECT id, cif_file, seed_msa_file, msa_file, hmm_file, rf_file, biomes_file, domain_architecture_file FROM mgnifam")
    rows = cursor.fetchall()
    conn.close()

    for row in rows:
        process_row(db_path, pipeline_results_dir, family_dir, row)

def main():
    if len(sys.argv) != 3:
        print("Usage: python3 append_sqlite_blobs.py <db.sqlite3> <pipeline_results_dir>")
        sys.exit(1)

    db_path              = sys.argv[1]
    pipeline_results_dir = sys.argv[2]
    family_dir_name      = "families"

    # Import files in parallel
    import_files(db_path, pipeline_results_dir, family_dir_name)
    
if __name__ == "__main__":
    main()
