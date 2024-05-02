import sys
import pandas as pd
import time

def parse_args():
    global linclust_input_file

    if not (len(sys.argv) == 2):
        print("Usage: python3 create_clusters_bookkeeping_df.py [Linclust TSV file]")
        sys.exit(1)
    
    linclust_input_file = sys.argv[1]  

def define_globals():
    global log_file, clusters_bookkeeping_df_pkl_file

    log_file                         = "pkl_log.txt"
    clusters_bookkeeping_df_pkl_file = "clusters_bookkeeping_df.pkl"
    
def create_clusters_bookkeeping_df():
    start_time = time.time()
    
    with open(log_file, 'w') as file:
        file.write("Starting creation of bookkeeping pkl file.\n")

    df = pd.read_csv(linclust_input_file, sep='\t', header=None, names=['representative', 'member'])

    with open(log_file, 'a') as file:
        file.write("DF read.\n")

    family_sizes = df.groupby('representative').size()

    with open(log_file, 'a') as file:
        file.write("Sizes calculated.\n")

    clusters_bookkeeping_df = df.set_index('representative').join(family_sizes.rename('size'), on='representative')

    with open(log_file, 'a') as file:
        file.write("Indexed.\n")

    clusters_bookkeeping_df.to_pickle(clusters_bookkeeping_df_pkl_file)

    with open(log_file, 'a') as file:
        file.write("create_clusters_bookkeeping_df: ")
        file.write(str(time.time() - start_time) + "\n")

    return clusters_bookkeeping_df

def main():
    parse_args()
    define_globals()

    create_clusters_bookkeeping_df()

if __name__ == "__main__":
    main()
