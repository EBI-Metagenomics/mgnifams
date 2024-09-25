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

def construct_file_path(base_dir, family_dir, col_iter, file_column):
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
    return os.path.join(base_dir, directory, file_column)

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

def update_blob_columns(conn, task):
    if len(task[0][1]) > 2**31 - 1:
            raise ValueError(f"BLOB size for row {task[0][2]} exceeds the maximum allowed limit of 2GB")
    try:
        cursor = conn.cursor()
        query = f"UPDATE mgnifam SET {task[0][0]} = ? WHERE id = ?"
        cursor.execute(query, (sqlite3.Binary(task[0][1]), task[0][2]))
        conn.commit()
        print(f"Updated {task[0][0]} for row {task[0][2]}")
    except sqlite3.Error as e:
        print(f"Failed to update row {task[0][2]}: {e}")

def process_row(base_dir, family_dir, row):
    row_id = row[0]
    file_columns = row[1:]
    tasks = []

    for i, file_column in enumerate(file_columns):
        if file_column is not None:
            file_path = construct_file_path(base_dir, family_dir, i, file_column)
            blob_column = get_blob_column(i)
            blob_data = read_file(file_path)

            if blob_data:
                tasks.append((blob_column, blob_data, row_id))
            else:
                print(f"Skipping update for row {row_id} due to read error")

    return tasks

def import_files(conn, base_dir, family_dir, start_row_id=0):
    cursor = conn.cursor()
    cursor.execute("SELECT id, cif_file, seed_msa_file, msa_file, hmm_file, rf_file, biomes_file, domain_architecture_file FROM mgnifam WHERE id >= ?", (start_row_id,))
    rows = cursor.fetchall()

    for row in rows:
        tasks = process_row(base_dir, family_dir, row)  # Process each row
        for task in tasks:
            update_blob_columns(conn, [task])  # Write each row immediately


def main():
    if len(sys.argv) != 5:
        print("Usage: python3 append_blobs_sqlite_single_thread.py <db.sqlite3> <output_dir> <families_dir> <start_row_id>")
        sys.exit(1)

    db_path      = sys.argv[1]
    base_dir     = sys.argv[2]
    family_dir   = sys.argv[3]
    start_row_id = int(sys.argv[4])

    # Open the database connection once
    conn = sqlite3.connect(db_path)

    # Import files, starting from the specified row ID
    import_files(conn, base_dir, family_dir, start_row_id=start_row_id)

    # Close the database connection
    conn.close()

if __name__ == "__main__":
    main()
