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
        0: "families/seed_msa",
        1: "families/msa"
    }
    
    directory = file_directory.get(col_iter)
    return os.path.join(base_dir, directory, file_column)

def get_blob_column(col_iter):
    blob_col_name = {
        0: "seed_msa_blob",
        1: "msa_blob"
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
    cursor.execute("SELECT id, seed_msa_file, msa_file FROM mgnifam")
    
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
