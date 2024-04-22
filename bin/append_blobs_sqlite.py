# Manually add the _blob columns first in sql:
# ALTER TABLE mgnifam
# ADD COLUMN cif_blob BLOB;

# ALTER TABLE mgnifam
# ADD COLUMN seed_msa_blob BLOB;

# ALTER TABLE mgnifam
# ADD COLUMN msa_blob BLOB;

# ALTER TABLE mgnifam
# ADD COLUMN hmm_blob BLOB;

# ALTER TABLE mgnifam
# ADD COLUMN rf_blob BLOB;

# ALTER TABLE mgnifam
# ADD COLUMN biomes_blob BLOB;

# ALTER TABLE mgnifam
# ADD COLUMN domain_architecture_blob BLOB;

import sqlite3
import os

def test_connection(conn): # test
    try:
        cursor = conn.cursor()
        cursor.execute("SELECT COUNT(*) FROM mgnifam")
        row_count = cursor.fetchone()[0]
        print("Connection successful! Number of rows in the table: ", row_count)
    except sqlite3.Error as e:
        print("Connection failed:", e)

def construct_file_path(base_dir, col_iter, file_column):
    # Mapping file columns to their respective directories
    file_directory = {
        0: "cif",
        1: "families/seed_msa",
        2: "families/msa",
        3: "families/hmm",
        4: "families/rf",
        5: "biome_sunburst",
        6: "pfams"
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
    with open(file_path, 'rb') as file:
        return file.read()

def update_blob_column(conn, column_name, blob_data, row_id):
    try:
        cursor = conn.cursor()
        query = f"UPDATE mgnifam SET {column_name} = ? WHERE id = ?"
        cursor.execute(query, (sqlite3.Binary(blob_data), row_id))
        conn.commit()
    except sqlite3.Error as e:
        print(f"Failed to update {column_name} for row {row_id}: {e}")

def import_files(conn, base_dir):
    cursor = conn.cursor()
    cursor.execute("SELECT id, cif_file, seed_msa_file, msa_file, hmm_file, rf_file, biomes_file, domain_architecture_file FROM mgnifam")
    
    for row in cursor.fetchall():
        row_id = row[0]
        file_columns = row[1:]  # Exclude the id column

        for i, file_column in enumerate(file_columns):
            if file_column is not None:  # Check if the file column is not NULL
                file_path = construct_file_path(base_dir, i, file_column)
                blob_column = get_blob_column(i)
                blob_data = read_file(file_path)
                update_blob_column(conn, blob_column, blob_data, row_id)
        

# Connect to the SQLite database
conn = sqlite3.connect('/home/vangelis/Desktop/Projects/mgnifams/DB/mgnifams.sqlite3')
# test_connection(conn)

base_dir = '/home/vangelis/Desktop/Projects/mgnifams-site-data_backup'

import_files(conn, base_dir)

# Close the connection
conn.close()
