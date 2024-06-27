import sqlite3
import os
from concurrent.futures import ThreadPoolExecutor, as_completed

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

def process_row(db_path, base_dir, family_dir, row):
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

def import_files(db_path, base_dir, family_dir, max_workers=8):
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()
    cursor.execute("SELECT id, cif_file, seed_msa_file, msa_file, hmm_file, rf_file, biomes_file, domain_architecture_file FROM mgnifam")
    rows = cursor.fetchall()
    conn.close()

    tasks = []
    for row in rows:
        tasks.extend(process_row(db_path, base_dir, family_dir, row))

    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = [executor.submit(update_blob_column, db_path, *task) for task in tasks]
        for future in as_completed(futures):
            future.result()  # This will raise exceptions if any occurred during processing

# Define database path and directories
db_path    = '/home/vangelis/Desktop/Projects/mgnifams/DB/mgnifams.sqlite3' # /nfs/production/rdf/metagenomics/users/vangelis/mgnifams/output/tables/mgnifams.sqlite3
base_dir   = '/home/vangelis/Desktop/Projects/mgnifams/output'              # /nfs/production/rdf/metagenomics/users/vangelis/mgnifams/output
family_dir = 'families'

# Import files in parallel
import_files(db_path, base_dir, family_dir)
