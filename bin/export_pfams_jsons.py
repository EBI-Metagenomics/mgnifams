import argparse
import configparser
import psycopg2
import os
import pandas as pd
import ast
import json

def read_config(filename='bin/db_config.ini'):
    config = configparser.ConfigParser()
    config.read(filename)
    return dict(config.items('database'))

def extract_mgyp(protein_name):
    parts = protein_name.split('/')
    if len(parts) > 1:
        first_part = parts[0].split('_')[0]
        return first_part
    elif "_" in protein_name:
        return protein_name.split('_')[0]
    return protein_name

def execute_query(cursor, query_items, mgyp_p_dir):
    mgyf_id = query_items[0][0]  # Get mgyf_id from the first item
    output_csv = f"{mgyp_p_dir}/{mgyf_id}.tsv"

    unique_col2_values = list(set(item[1] for item in query_items))
    sql_query = "SELECT mgyp, metadata FROM sequence_explorer_protein WHERE mgyp IN ({})".format(
        ",".join([f"'{col2}'" for col2 in unique_col2_values])
    )
    print(mgyf_id)
    # print(sql_query)
    cursor.execute(sql_query)
    rows = cursor.fetchall()

    with open(output_csv, 'w') as file:
        for row in rows:
            metadata = row[1]
            if 'p' in metadata:
                mgyp = row[0]
                file.write(f"{mgyp}\t{metadata['p']}\n")

def is_above_family_id(fam_name, above_family_id):
    fam_id = fam_name.replace("mgnifam", "")
    fam_id = float(fam_id)

    if fam_id > int(above_family_id):
        return True
    else:
        return False

def query_sequence_explorer_protein(cursor, edge_list_file, above_family_id, mgyp_p_dir):
    with open(edge_list_file, 'r') as file:
        previous_col1 = None
        query_items = []
        for line in file:
            col1, col2 = line.strip().split('\t')
            if (is_above_family_id(col1, above_family_id)):
                if col1 != previous_col1:
                    if query_items:  # Execute previous query if items exist
                        execute_query(cursor, query_items, mgyp_p_dir)
                    query_items = []
                    previous_col1 = col1
                parsed_col2 = extract_mgyp(col2)
                query_items.append((col1, parsed_col2))

        # Execute the last query
        if query_items:
            execute_query(cursor, query_items, mgyp_p_dir)

def calculate_mgnifams_start(protein_id):
    number_of_underscores = protein_id.count('_')
    if (number_of_underscores == 0): # 250671917
        start = 1
    elif (number_of_underscores == 1): # 250671917/1_50
        parts = protein_id.split('/')
        start = parts[1].split('_')[0]
    elif (number_of_underscores == 2): # 250671917_50_150
        start = protein_id.split('_')[1]
    elif (number_of_underscores == 3): # 250671917_50_200/2_34
        start = int(protein_id.split('_')[1])
        region = protein_id.split('/')[1].split('_')
        start = start + int(region[0]) - 1

    return start

def get_edgelist_family_subset(clusters_df, family_name):
    subset_clusters_df = clusters_df[clusters_df['family_name'] == family_name].copy()
    subset_clusters_df['mgyp'] = subset_clusters_df['protein_name'].apply(extract_mgyp)
    subset_clusters_df['mgnifams_start'] = subset_clusters_df['protein_name'].apply(calculate_mgnifams_start)
    return subset_clusters_df

def construct_domain_architecture(mgyp, pfams, matched_rows):
    pfams = ast.literal_eval(pfams)
    fam_names = []
    start_points = []
    for pfam in pfams:
        fam_names.append(pfam[0])
        start_points.append(pfam[3])

    for index, row in matched_rows.iterrows():
        fam_names.append(row['family_name'])
        start_points.append(int(row['mgnifams_start']))

    # Sort both arrays based on start_points
    sorted_data = sorted(zip(start_points, fam_names))
    start_points, fam_names = zip(*sorted_data)

    domain_architecture = '\t'.join(map(str, fam_names))

    return domain_architecture

def string_to_hex_color(s):
    hash_val = 0
    for char in s:
        hash_val = ord(char) + ((hash_val << 5) - hash_val)
    
    color = '#'
    for i in range(3):
        value = (hash_val >> (i * 8)) & 0xFF
        color += ('00' + format(value, 'x'))[-2:]

    return color

def write_out_json(element_counts, output_filename):
    architecture_containers = []
    for architecture_text, count in element_counts.items():
        domains = []
        for domain in architecture_text.split('\t'):
            domains.append({"name": domain, "color": string_to_hex_color(domain)})
        architecture_containers.append({"architecture_text": str(count), "domains": domains})

    output_json = {"architecture_containers": architecture_containers}

    # Convert to JSON string
    with open(output_filename, 'w') as f:
        json.dump(output_json, f, indent=4)

def construct_pfams_json(edge_list_file, mgyp_p_dir, json_id_dir):
    clusters_df = pd.read_csv(edge_list_file, delimiter='\t', header=None, names=['family_name', 'protein_name'])

    files = os.listdir(mgyp_p_dir)
    for file_name in files:
        print(file_name)
        family_name = file_name.split(".")[0]
        subset_clusters_df = get_edgelist_family_subset(clusters_df, family_name)
        file_path = os.path.join(mgyp_p_dir, file_name)
        family_domain_architectures = []
        counter = 0
        with open(file_path, 'r') as file:
            for line in file:
                if counter % 1000 == 0:
                    print(counter)
                # if counter == 10:
                #     break
                parts = line.strip().split('\t')
                if len(parts) == 2:
                    mgyp = parts[0]
                    matched_rows = subset_clusters_df[subset_clusters_df['mgyp'] == mgyp]
                    if not matched_rows.empty:
                        pfams = parts[1]
                        domain_architecture = construct_domain_architecture(mgyp, pfams, matched_rows)
                        family_domain_architectures.append(domain_architecture)
                counter += 1

        element_counts = pd.Series(family_domain_architectures).value_counts()
        family_name = family_name.replace("mgnifam", "mgnfam")
        output_filename = os.path.join(json_id_dir, f"{family_name}_domains.json")
        write_out_json(element_counts, output_filename)

def subset_json(json_items, threshold = 10):
    return json_items['architecture_containers'][:threshold]

def execute_pfam_translation_query(cursor, pfam_id):
    sql_query = f"SELECT name FROM sequence_explorer_pfam WHERE id = '{pfam_id}'"
    cursor.execute(sql_query)
    result = cursor.fetchall()
    return(result)

def append_link(domain):
    pass

def write_out_final_json(json_object, output_filename):
    final_json_object = {"architecture_containers": json_object}
    with open(output_filename, 'w') as f:
        json.dump(final_json_object, f, indent=4)

def hex_to_rgb(hex_color):
    # Remove '#' from the beginning of the hex color
    hex_color = hex_color.lstrip('#')
    
    # Convert the hex color to RGB
    r = int(hex_color[0:2], 16)
    g = int(hex_color[2:4], 16)
    b = int(hex_color[4:6], 16)
    
    return (r, g, b)

def calculate_luminosity(rgb):
    # Convert RGB to linear RGB
    def linearize(color):
        c = color / 255.0
        if c <= 0.03928:
            return c / 12.92
        return ((c + 0.055) / 1.055) ** 2.4
    
    r, g, b = rgb
    # Calculate relative luminance
    luminance = 0.2126 * linearize(r) + 0.7152 * linearize(g) + 0.0722 * linearize(b)
    return luminance

def decide_font_color(hex_color):
    # Convert hex color to RGB
    rgb = hex_to_rgb(hex_color)
    
    # Calculate luminosity
    luminosity = calculate_luminosity(rgb)
    
    # Decide font color based on luminosity
    if luminosity > 0.2:
        return 'black'
    else:
        return 'white'

def construct_name(mgnifam_id):
    name = mgnifam_id.replace('mgnifam', 'MGYF')
    number_str = name[4:]
    number = int(number_str)
    formatted_number = '{:010d}'.format(number)
    final_name = 'MGYF' + formatted_number
    
    return final_name

def translate_pfams(cursor, json_id_dir, out_dir):
    files = os.listdir(json_id_dir)
    for file_name in files:
        print(file_name)
        output_filename = os.path.join(out_dir, file_name)
        if os.path.exists(output_filename):
            print(f"Skipping {file_name} as it already exists in {out_dir}")
            continue

        full_path = os.path.join(json_id_dir, file_name)
        with open(full_path, 'r') as file:
            json_data = json.load(file)
            top_10_architecture_containers = subset_json(json_data, 10)
            for architecture_container in top_10_architecture_containers:
                for domain in architecture_container['domains']:
                    if 'mgnifam' not in domain['name']:
                        domain['link'] = f'https://www.ebi.ac.uk/interpro/entry/pfam/{domain["name"]}/domain_architecture/'
                        domain['name'] = execute_pfam_translation_query(cursor, domain['name'])[0][0]
                        domain['color'] = string_to_hex_color(domain['name'])
                    else:
                        domain['link'] = f'http://mgnifams-demo.mgnify.org/details/?id={construct_name(domain["name"])}'
                    domain['font_color'] = decide_font_color(domain['color'])

            write_out_final_json(top_10_architecture_containers, output_filename)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Query the PostgreSQL database and create pfam + mgnifam domain architecture JSON files.")
    parser.add_argument("config_file", help="Path to the configuration file for the database secrets")
    parser.add_argument("edge_list_file", help="Path to the edge list file with two columns")
    parser.add_argument("above_family_id", help="Threshold for family id, to not recalculate same ones")
    # python3 bin/export_pfams_jsons.py bin/db_config.ini /home/vangelis/Desktop/Projects/mgnifams-site-data_backup/families/updated_refined_families.tsv 0
    args = parser.parse_args()

    pfams_dir = "tmp/pfams"
    os.makedirs(pfams_dir, exist_ok=True)
    tmp_dir = os.path.join(pfams_dir, "tmp")
    os.makedirs(tmp_dir, exist_ok=True)
    mgyp_p_dir = os.path.join(tmp_dir, "mgyp_p")
    json_id_dir = os.path.join(tmp_dir, "json_id")
    os.makedirs(mgyp_p_dir, exist_ok=True)
    os.makedirs(json_id_dir, exist_ok=True)
    out_dir = os.path.join(pfams_dir, "out")
    os.makedirs(out_dir, exist_ok=True)

    db_params = read_config(args.config_file)
    conn = psycopg2.connect(**db_params)
    cursor = conn.cursor()

    print("query_sequence_explorer_protein")
    query_sequence_explorer_protein(cursor, args.edge_list_file, args.above_family_id, mgyp_p_dir)
    print("construct_pfams_json")
    construct_pfams_json(args.edge_list_file, mgyp_p_dir, json_id_dir)
    print("translate_pfams")
    translate_pfams(cursor, json_id_dir, out_dir)

    cursor.close()
    conn.close()
