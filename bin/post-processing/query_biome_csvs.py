import argparse
import configparser
import psycopg2
import pandas as pd
import os

def read_config(filename='bin/db_config.ini'):
    config = configparser.ConfigParser()
    config.read(filename)
    return dict(config.items('database'))

def execute_protein_query(cursor, query_items, mgyp_b_dir, b_counts_dir):
    mgyf_id = query_items[0][0]  # Get mgyf_id from the first item
    output_csv = f"{mgyp_b_dir}/{mgyf_id}.csv"

    unique_col2_values = list(set(item[1] for item in query_items))
    sql_query = "SELECT * FROM sequence_explorer_protein WHERE mgyp IN ({})".format(
        ",".join([f"'{col2}'" for col2 in unique_col2_values])
    )
    print(mgyf_id)
    # print(sql_query)
    cursor.execute(sql_query)
    rows = cursor.fetchall()

    with open(output_csv, 'w') as file:
        for row in rows:
            metadata = row[4]
            if 'b' in metadata:
                mgyp = row[0]
                b_ids = metadata['b']
                for b_id_pair in b_ids:
                    first_number = b_id_pair[0]
                    file.write(f"{mgyp},{first_number}\n")

    # Calculate and save b_counts_df for each family
    df = pd.DataFrame({'b_id': [b_id_pair[0] for row in rows for b_id_pair in row[4].get('b', [])]})
    b_counts = df['b_id'].value_counts()
    b_counts_df = b_counts.reset_index()
    b_counts_df.columns = ['b_id', 'count']
    b_counts_df.to_csv(os.path.join(b_counts_dir, f"{mgyf_id}_b_counts.csv"), index=False)

def parse_col2(col2):
    parts = col2.split('/')
    if len(parts) > 1:
        first_part = parts[0].split('_')[0]
        return first_part
    elif "_" in col2:
        return col2.split('_')[0]
    return col2

def is_above_family_id(fam_name, above_family_id):
    fam_id = fam_name.replace("mgnifam", "")
    fam_id = float(fam_id)

    if fam_id > int(above_family_id):
        return True
    else:
        return False

def query_sequence_explorer_protein(cursor, edge_list_file, above_family_id, mgyp_b_dir, b_counts_dir):
    with open(edge_list_file, 'r') as file:
        previous_col1 = None
        query_items = []
        for line in file:
            col1, col2 = line.strip().split('\t')
            if (is_above_family_id(col1, above_family_id)):
                if col1 != previous_col1:
                    if query_items:  # Execute previous query if items exist
                        execute_protein_query(cursor, query_items, mgyp_b_dir, b_counts_dir)
                    query_items = []
                    previous_col1 = col1
                parsed_col2 = parse_col2(col2)
                query_items.append((col1, parsed_col2))

        # Execute the last query
        if query_items:
            execute_protein_query(cursor, query_items, mgyp_b_dir, b_counts_dir)

def execute_biome_translation_query(cursor, b_id):
    sql_query = f"SELECT name FROM sequence_explorer_biome WHERE id = {b_id}"
    cursor.execute(sql_query)
    result = cursor.fetchall()
    return(result)

def get_parent(biome_path):
    parts = biome_path.split(':')
    if len(parts) <= 1:
        return ''
    parent = ':'.join(parts[:-1])
    return parent

def get_label(biome_path):
    parts = biome_path.split(':')
    if len(parts) == 0:
        return ''  # No label found
    return parts[-1]

def append_parents(parent_names, biome_names, out_file):
    if ('root' not in biome_names):
        out_file.write("root,root,,0\n")
        biome_names.append("root")
    for parent_name in parent_names:
        while parent_name not in biome_names and parent_name != '':
            biome_names.append(parent_name)
            grandparent_name = get_parent(parent_name)
            label = get_label(parent_name)
            out_file.write(f"{parent_name},{label},{grandparent_name},0\n")
            parent_name = grandparent_name

def query_sequence_explorer_biome(cursor, counts_dir, out_dir):
    files = os.listdir(counts_dir)
    for file_name in files:
        file_path = os.path.join(counts_dir, file_name)
        df = pd.read_csv(file_path)
        biome_names = []
        parent_names = []
        file_name = file_name.replace("mgnifam", "mgnfam")
        print(file_name)
        with open(os.path.join(out_dir, file_name), 'w') as out_file:
            out_file.write("ids,labels,parents,counts\n")
            for index, row in df.iterrows():
                b_id = row['b_id']
                result = execute_biome_translation_query(cursor, b_id)
                biome_name = result[0][0]
                parent_name = get_parent(biome_name)
                label = get_label(biome_name)
                biome_names.append(biome_name)
                parent_names.append(parent_name)
                out_file.write(f"{biome_name},{label},{parent_name},{row['count']}\n")
            
            append_parents(parent_names, biome_names, out_file)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Query the PostgreSQL database and create biome sunburst count CSVs.")
    parser.add_argument("config_file", help="Path to the configuration file for the database secrets")
    parser.add_argument("edge_list_file", help="Path to the edge list file with two columns")
    parser.add_argument("above_family_id", help="Threshold for family id, to not recalculate same ones")
    # python3 bin/post-processing/query_biome_csvs.py bin/db_config.ini /home/vangelis/Desktop/Projects/mgnifams-site-data_backup/families/updated_refined_families.tsv 0
    args = parser.parse_args()

    biomes_dir = "tmp/biomes"
    os.makedirs(biomes_dir, exist_ok=True)
    tmp_dir = os.path.join(biomes_dir, "tmp")
    os.makedirs(tmp_dir, exist_ok=True)
    mgyp_b_dir = os.path.join(tmp_dir, "mgyp_b")
    b_counts_dir = os.path.join(tmp_dir, "b_counts")
    os.makedirs(mgyp_b_dir, exist_ok=True)
    os.makedirs(b_counts_dir, exist_ok=True)
    out_dir = os.path.join(biomes_dir, "out")
    os.makedirs(out_dir, exist_ok=True)
    
    db_params = read_config(args.config_file)
    conn = psycopg2.connect(**db_params)
    cursor = conn.cursor()

    print("query_sequence_explorer_protein")
    query_sequence_explorer_protein(cursor, args.edge_list_file, args.above_family_id, mgyp_b_dir, b_counts_dir)
    print("query_sequence_explorer_biome")
    query_sequence_explorer_biome(cursor, b_counts_dir, out_dir)

    cursor.close()
    conn.close()
