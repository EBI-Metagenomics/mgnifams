import sqlite3
import sys
import os
from concurrent.futures import ThreadPoolExecutor, as_completed

def get_basenames(directory):
    # Ensure the directory exists
    if not os.path.exists(directory):
        raise FileNotFoundError(f"The directory {directory} does not exist.")

    # Get the list of basenames without the extensions
    basenames = [os.path.splitext(os.path.basename(file))[0] 
                for file in os.listdir(directory) 
                if os.path.isfile(os.path.join(directory, file))]

    return basenames

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

def import_files(db_path, base_dir, family_dir, basenames, max_workers=8):
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()
    cursor.execute("SELECT id, cif_file, seed_msa_file, msa_file, hmm_file, rf_file, biomes_file, domain_architecture_file FROM mgnifam")
    rows = cursor.fetchall()
    conn.close()

    tasks = []
    print(f"Using {max_workers} parallel import jobs\n")
    for row in rows:
        row_id = str(row[0])
        if row_id in basenames:
            tasks.extend(process_row(db_path, base_dir, family_dir, row))

    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = [executor.submit(update_blob_column, db_path, *task) for task in tasks]
        for future in as_completed(futures):
            future.result()  # This will raise exceptions if any occurred during processing

def main():
    if len(sys.argv) != 5:
        print("Usage: python3 ${params.scriptDir}/post-processing/append_blobs_sqlite_parallel.py <db.sqlite3> <output_dir> <families_dir> <threads>")
        sys.exit(1)

    db_path     = sys.argv[1]
    base_dir    = sys.argv[2]
    family_dir  = sys.argv[3]
    max_workers = int(sys.argv[4])

    # Get chunks to check from domains chunked folder
    domain_chunk_dir = os.path.join(base_dir, "post-processing", "domain_results")
    basenames = get_basenames(domain_chunk_dir)
    print(basenames)
    # Import files in parallel
    import_files(db_path, base_dir, family_dir, basenames, max_workers)
    
if __name__ == "__main__":
    main()
