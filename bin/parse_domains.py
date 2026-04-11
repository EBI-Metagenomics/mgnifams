#!/usr/bin/env python3

import argparse
import os
import ast
import json
from collections import Counter, defaultdict

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

def construct_domain_architecture(pfams, family_id, mgnifam_starts):
    pfams        = ast.literal_eval(pfams)
    fam_names    = []
    start_points = []

    for pfam in pfams:
        fam_names.append(pfam[0])
        start_points.append(pfam[3])

    for mgnifam_start in mgnifam_starts:
        fam_names.append(str(family_id))
        start_points.append(int(mgnifam_start))

    sorted_data             = sorted(zip(start_points, fam_names))
    start_points, fam_names = zip(*sorted_data)
    domain_architecture     = '\t'.join(map(str, fam_names))

    return domain_architecture

def construct_solo_domain_architecture(family_id, mgnifam_starts):
    fam_names           = [str(family_id)] * len(mgnifam_starts)
    domain_architecture = '\t'.join(fam_names)

    return domain_architecture

def count_domain_architectures(file_path, family_id, mgyp_lookup):
    """mgyp_lookup: {mgyp: [mgnifam_start, ...]} for this family only."""
    family_domain_architectures = []

    with open(file_path, 'r') as file:
        for line in file:
            parts           = line.strip().split('\t')
            mgyp            = parts[0]
            mgnifam_starts  = mgyp_lookup.get(mgyp, [])

            if len(parts) == 3: # aka pfams not empty
                pfams               = parts[2]
                domain_architecture = construct_domain_architecture(pfams, family_id, mgnifam_starts)
            else: # no pfams
                domain_architecture = construct_solo_domain_architecture(family_id, mgnifam_starts)

            family_domain_architectures.append(domain_architecture)

    return Counter(family_domain_architectures)

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

    for architecture_text, count in domain_architecture_counts.most_common():
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
    hex_color = hex_color.lstrip('#')

    r = int(hex_color[0:2], 16)
    g = int(hex_color[2:4], 16)
    b = int(hex_color[4:6], 16)

    return (r, g, b)

def calculate_luminosity(rgb):
    def linearize(color):
        c = color / 255.0
        if c <= 0.03928:
            return c / 12.92
        return ((c + 0.055) / 1.055) ** 2.4

    r, g, b = rgb
    luminance = 0.2126 * linearize(r) + 0.7152 * linearize(g) + 0.0722 * linearize(b)

    return luminance

def decide_font_color(hex_color):
    rgb        = hex_to_rgb(hex_color)
    luminosity = calculate_luminosity(rgb)

    if luminosity > 0.2:
        return 'black'
    else:
        return 'white'

def load_pfam_mapping(pfam_mapping_file):
    pfam_mapping = {}
    with open(pfam_mapping_file, 'r') as f:
        for line in f:
            pfam_id, name = line.strip().split('\t', 1)
            pfam_mapping[pfam_id] = name
    return pfam_mapping

def translate_architecture(architecture_json, pfam_mapping):
    translated_top_json = subset_json(architecture_json)

    for architecture_container in translated_top_json:
        for domain in architecture_container['domains']:
            if 'PF' in domain['id']: # pfam
                domain['link']  = f'https://www.ebi.ac.uk/interpro/entry/pfam/{domain["id"]}'
                domain['name']  = pfam_mapping[domain['id']]
                domain['color'] = string_to_hex_color(domain['name'])
            else: # mgnifam
                domain['link'] = f'http://mgnifams-demo.mgnify.org/details/{construct_name(int(domain["id"]))}'
                domain['name'] = f'MGnifam{domain["id"]}'
            domain['font_color'] = decide_font_color(domain['color'])

    return translated_top_json

def write_out(translated_json, out_file):
    translated_json = {"architecture_containers": translated_json}

    with open(out_file, 'w') as f:
        json.dump(translated_json, f, indent=4)

def load_family_data(refined_families_file):
    """Load all family-protein mappings into a nested dict without requiring sorted input.
    Returns {family_id: {mgyp: [mgnifam_start, ...]}}."""
    family_data = defaultdict(lambda: defaultdict(list))

    with open(refined_families_file, 'r') as f:
        for line in f:
            family_id_str, protein_name = line.strip().split('\t')
            family_id     = int(family_id_str)
            mgyp          = extract_mgyp(protein_name)
            mgnifam_start = calculate_mgnifam_start(protein_name)
            family_data[family_id][mgyp].append(mgnifam_start)

    return family_data

def process_family(family_id, mgyp_lookup, query_results_dir, pfam_mapping, output_dir):
    tsv_path = os.path.join(query_results_dir, f"{family_id}.tsv")
    if not os.path.exists(tsv_path):
        return

    domain_architecture_counts = count_domain_architectures(tsv_path, family_id, mgyp_lookup)
    architecture_json           = construct_architecture_json(domain_architecture_counts)
    translated_top_json         = translate_architecture(architecture_json, pfam_mapping)

    out_file = os.path.join(output_dir, f"{family_id}.json")
    write_out(translated_top_json, out_file)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Parse the protein query TSV results into domain architecture JSONs.")
    parser.add_argument("--query_results",    help="Path to the output query results dir generated by the query_mgnprotein_db module")
    parser.add_argument("--pfam_mapping",     help="Path to the biome id to names mapping tsv file")
    parser.add_argument("--refined_families", help="Path to the family-proteins TSV file with two columns")
    parser.add_argument("--output_dir",       help="Path to the output directory with per family JSON domain architectures")

    args = parser.parse_args()

    if not os.path.exists(args.query_results):
        raise FileNotFoundError(f"The folder {args.query_results} does not exist.")

    if not os.path.exists(args.pfam_mapping):
        raise FileNotFoundError(f"The file {args.pfam_mapping} does not exist.")

    if not os.path.exists(args.refined_families):
        raise FileNotFoundError(f"The file {args.refined_families} does not exist.")

    pfam_mapping = load_pfam_mapping(args.pfam_mapping)
    family_data  = load_family_data(args.refined_families)

    os.makedirs(args.output_dir, exist_ok=True)
    for tsv in os.listdir(args.query_results):
        family_id   = int(os.path.splitext(tsv)[0])
        mgyp_lookup = family_data.get(family_id, {})
        process_family(family_id, mgyp_lookup, args.query_results, pfam_mapping, args.output_dir)
