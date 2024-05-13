import argparse
import configparser
import os
import psycopg2

def read_config(db_config_file):
    config = configparser.ConfigParser()
    config.read(db_config_file)
    return dict(config.items('database'))

def extract_sequence_id(sequence_id_with_region):
    parts = sequence_id_with_region.split('/')
    if len(parts) > 1:
        return parts[0].split('_')[0]
    elif "_" in sequence_id_with_region:
        return sequence_id_with_region.split('_')[0]

    return sequence_id_with_region

def write_out(family_id, rows):
    output_csv = f"{output_dir}/{family_id}.csv"
    with open(output_csv, 'w') as file:
        for row in rows:
            mgyp = row[0]
            metadata = row[1]
            meta_b = metadata.get('b', 'null')
            meta_p = metadata.get('p', 'null')
            file.write(f"{mgyp},{meta_b},{meta_p}\n")

def execute_query(cursor, family_id, query_sequences):
    unique_query_sequences = list(set(query_sequences))
    sql_query = "SELECT mgyp, metadata FROM sequence_explorer_protein WHERE mgyp IN ({})".format(
        ",".join([f"'{query_sequence}'" for query_sequence in unique_query_sequences])
    )
    # Debugging:
    # print(mgyf_id)
    # print(sql_query)
    cursor.execute(sql_query)
    rows = cursor.fetchall()
    write_out(family_id, rows)

def query_sequence_explorer_protein(cursor, family_proteins_file):
    with open(family_proteins_file, 'r') as file:
        previous_family_id = None
        query_sequences = []
        for line in file:
            family_id, sequence_id_with_region = line.strip().split('\t')
            if family_id != previous_family_id:
                if query_sequences:  # Execute previous query if items exist
                    execute_query(cursor, previous_family_id, query_sequences)
                query_sequences = []
                previous_family_id = family_id
            sequence_id = extract_sequence_id(sequence_id_with_region)
            query_sequences.append(sequence_id)

        # Execute the last query
        if query_sequences:
            execute_query(cursor, previous_family_id, query_sequences)

def query_sequence_explorer_biome(cursor):
    sql_query = "SELECT id, name FROM sequence_explorer_biome"
    cursor.execute(sql_query)
    rows = cursor.fetchall()

    output_tsv = "biome_mapping.tsv"
    with open(output_tsv, 'w') as file:
        for row in rows:
            file.write(f"{row[0]}\t{row[1]}\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Query the PostgreSQL database MGnifams proteins data.")
    parser.add_argument("db_config_file",       help="Path to the configuration file for the database secrets")
    parser.add_argument("family_proteins_file", help="Path to the tsv file with families and respective proteins")
    
    args = parser.parse_args()

    db_params = read_config(args.db_config_file)
    conn      = psycopg2.connect(**db_params)
    cursor    = conn.cursor()

    global output_dir
    output_dir = "query_results"
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
        
    query_sequence_explorer_protein(cursor, args.family_proteins_file)
    query_sequence_explorer_biome(cursor)

    cursor.close()
    conn.close()
