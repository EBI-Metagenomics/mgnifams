import sys
import os

def parse_args():
    global arg_families_dir
        
    if not (len(sys.argv) == 2):
        print("Incorrect number of args.")
        sys.exit(1)

    arg_families_dir = sys.argv[1]

def pool_directory(input_dir, output_filename):
    path_to_folder = os.path.join(arg_families_dir, input_dir)
    output_file    = os.path.join(arg_families_dir, output_filename)    
    
    with open(output_file, 'w') as outfile:
        for filename in os.listdir(path_to_folder):
            filepath = os.path.join(path_to_folder, filename)
            with open(filepath, 'r') as infile:
                base_filename = os.path.splitext(filename)[0]
                for line in infile:
                    value = line.strip()
                    if value:  # Check if the line is not empty
                        outfile.write(f"{base_filename}_{value}\n")

def pool_clusters_directory(input_dir, output_filename):
    path_to_folder = os.path.join(arg_families_dir, input_dir)
    output_file    = os.path.join(arg_families_dir, output_filename)
    
    with open(output_file, 'w') as outfile:
        for filename in os.listdir(path_to_folder):
            filepath = os.path.join(path_to_folder, filename)
            with open(filepath, 'r') as infile:
                content = infile.read().strip()  # Read and strip leading/trailing whitespace
                if content:  # Only write if the file is not empty
                    outfile.write(content + '\n')

def main():
    parse_args()

    pool_directory("refined_families", "refined_families.tsv")
    pool_directory("family_metadata", "family_metadata.csv")
    pool_directory("converged_families", "converged_families.txt")
    pool_clusters_directory("successful_clusters", "successful_clusters.txt")
    pool_clusters_directory("discarded_clusters", "discarded_clusters.txt")
    
if __name__ == "__main__":
    main()