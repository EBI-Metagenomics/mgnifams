import argparse
import os
import pandas as pd
import ast
import json
from concurrent.futures import ThreadPoolExecutor

def extract_mgyp(protein_name):
    parts = protein_name.split('/')

    if len(parts) > 1:
        first_part = parts[0].split('_')[0]
        return first_part
    elif "_" in protein_name:
        return protein_name.split('_')[0]

    return protein_name

def calculate_mgnifam_start(protein_name):
    number_of_underscores = protein_name.count('_')

    if (number_of_underscores == 0): # 250671917
        start = 1
    elif (number_of_underscores == 1): # 250671917/1_50
        parts = protein_name.split('/')
        start = parts[1].split('_')[0]
    elif (number_of_underscores == 2): # 250671917_50_150
        start = protein_name.split('_')[1]
    elif (number_of_underscores == 3): # 250671917_50_200/2_34
        start  = int(protein_name.split('_')[1])
        region = protein_name.split('/')[1].split('_')
        start  = start + int(region[0]) - 1

    return start

def get_refined_families_subset(clusters_df, family_id):
    refined_families_subset                  = clusters_df[clusters_df['family_id'] == family_id].copy()
    refined_families_subset['mgyp']          = refined_families_subset['protein_name'].apply(extract_mgyp)
    refined_families_subset['mgnifam_start'] = refined_families_subset['protein_name'].apply(calculate_mgnifam_start)

    return refined_families_subset

def construct_domain_architecture(pfams, matched_rows):
    pfams        = ast.literal_eval(pfams)
    fam_names    = []
    start_points = []

    for pfam in pfams:
        fam_names.append(pfam[0])
        start_points.append(pfam[3])

    for index, row in matched_rows.iterrows():
        fam_names.append(str(row['family_id']))
        start_points.append(int(row['mgnifam_start']))

    # Sort both arrays based on start_points
    sorted_data             = sorted(zip(start_points, fam_names))
    start_points, fam_names = zip(*sorted_data)
    domain_architecture     = '\t'.join(map(str, fam_names))

    return domain_architecture

def construct_solo_domain_architecture(matched_rows):
    fam_names = []

    for index, row in matched_rows.iterrows():
        fam_names.append(str(row['family_id']))

    domain_architecture = '\t'.join(map(str, fam_names))

    return domain_architecture

def count_domain_architectures(file_path, family_id, refined_families_df):
    refined_families_subset     = get_refined_families_subset(refined_families_df, family_id)
    family_domain_architectures = []

    with open(file_path, 'r') as file:
        for line in file:
            parts        = line.strip().split('\t')
            mgyp         = parts[0]
            matched_rows = refined_families_subset[refined_families_subset['mgyp'] == mgyp]
            
            if len(parts) == 3: # aka pfams not empty
                pfams               = parts[2]
                domain_architecture = construct_domain_architecture(pfams, matched_rows)
            else: # no pfams
                domain_architecture = construct_solo_domain_architecture(matched_rows)
            
            family_domain_architectures.append(domain_architecture)

    domain_architecture_counts = pd.Series(family_domain_architectures).value_counts()

    return domain_architecture_counts

def string_to_hex_color(s):
    hash_val = 0

    for char in s:
        hash_val = ord(char) + ((hash_val << 5) - hash_val)
    
    color = '#'

    for i in range(3):
        value = (hash_val >> (i * 8)) & 0xFF
        color += ('00' + format(value, 'x'))[-2:]

    return color

def construct_architecture_json(domain_architecture_counts):
    architecture_containers = []

    for architecture_text, count in domain_architecture_counts.items():
        domains = []

        for domain in architecture_text.split('\t'):
            domains.append({"id": domain, "color": string_to_hex_color(domain)})
        architecture_containers.append({"architecture_text": str(count), "domains": domains})

    output_json = {"architecture_containers": architecture_containers}

    return output_json

def subset_json(json_items, threshold = 15):
    return json_items['architecture_containers'][:threshold]

def construct_name(mgnifam_id):
    formatted_number = '{:010d}'.format(mgnifam_id)
    mgnifam_id       = 'MGYF' + formatted_number
    
    return mgnifam_id

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

def translate_architecture(architecture_json, pfam_mapping_df):
    translated_top_json = subset_json(architecture_json)
    
    for architecture_container in translated_top_json:
        for domain in architecture_container['domains']:
            if 'PF' in domain['id']: # pfam
                domain['link']  = f'https://www.ebi.ac.uk/interpro/entry/pfam/{domain["id"]}'
                domain['name']  = pfam_mapping_df.loc[pfam_mapping_df['id'] == domain["id"], 'name'].values[0]
                domain['color'] = string_to_hex_color(domain['name'])
            else: # mgnifam
                domain['link'] = f'http://mgnifams-demo.mgnify.org/details/?id={construct_name(int(domain["id"]))}'
                domain['name'] = f'MGnifam{domain["id"]}'
            domain['font_color'] = decide_font_color(domain['color'])

    return translated_top_json

def write_out(translated_json, out_file):
    translated_json = {"architecture_containers": translated_json}
    
    with open(out_file, 'w') as out_file:
        json.dump(translated_json, out_file, indent=4)

def process_tsv(file_path, family_id, refined_families_df, pfam_mapping_df, outdir):
    domain_architecture_counts = count_domain_architectures(file_path, family_id, refined_families_df)
    architecture_json = construct_architecture_json(domain_architecture_counts)
    translated_top_json = translate_architecture(architecture_json, pfam_mapping_df)

    out_file = os.path.join(outdir, f"{family_id}.json")
    write_out(translated_top_json, out_file)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Parse the protein query TSV results into domain architecture JSONs.")
    parser.add_argument("postprocessing_dir", help="Path to the output post-processing files generated by the query_mgnprotein_db module")
    parser.add_argument("refined_families_file", help="Path to the family-proteins TSV file with two columns")
    parser.add_argument("max_workers", help="num threads")
    # python3 bin/post-processing/parse_domains.py /home/vangelis/Desktop/Projects/mgnifams/output/post-processing /home/vangelis/Desktop/Projects/mgnifams/output/families/refined_families.tsv
    args = parser.parse_args()

    query_results_dir = os.path.join(args.postprocessing_dir, "query_results")
    if not os.path.exists(query_results_dir):
        raise FileNotFoundError(f"The folder {query_results_dir} does not exist.")

    pfam_mapping_path = os.path.join(args.postprocessing_dir, "pfam_mapping.tsv")
    if not os.path.exists(pfam_mapping_path):
        raise FileNotFoundError(f"The file {pfam_mapping_path} does not exist.")

    pfam_mapping_df = pd.read_csv(pfam_mapping_path, sep='\t', header=None, names=['id', 'name'])

    if not os.path.exists(args.refined_families_file):
        raise FileNotFoundError(f"The file {args.refined_families_file} does not exist.")

    query_results_files = os.listdir(query_results_dir)
    refined_families_df = pd.read_csv(args.refined_families_file, sep='\t', header=None, names=['family_id', 'protein_name'])

    max_workers = int(args.max_workers)
    print(f"Using {max_workers} parallel parse jobs\n")

    outdir = "domain_results"
    os.makedirs(outdir, exist_ok=True)

    # Use ThreadPoolExecutor to run tasks in parallel
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = []
        for tsv in query_results_files:
            file_path = os.path.join(query_results_dir, tsv)
            family_id = int(os.path.splitext(os.path.basename(file_path))[0])
            futures.append(executor.submit(process_tsv, file_path, family_id, refined_families_df, pfam_mapping_df, outdir))

        # Optionally, wait for all futures to complete
        for future in futures:
            future.result()
