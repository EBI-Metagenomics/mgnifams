import sys
import shutil
import pickle
import pandas as pd

import time # benchmarking, TODO remove

def parse_args():
    if not (len(sys.argv) == 2):
        print("Usage: python3 create_clusters_bookkeeping_df.py [Linclust TSV file]")
        sys.exit(1)

    globals().update({
        "linclust_input_file" : sys.argv[1]
    })    

def define_globals():
    globals().update({
        "log_file"                         : "log.txt",
        "clusters_bookkeeping_df_pkl_file" : "clusters_bookkeeping_df.pkl"
    })
    
def create_clusters_bookkeeping_df():
    start_time = time.time()

    df = pd.read_csv(linclust_input_file, sep='\t', header=None, names=['representative', 'member'])
    family_sizes = df.groupby('representative').size()
    clusters_bookkeeping_df = df.set_index('representative').join(family_sizes.rename('size'), on='representative')
    clusters_bookkeeping_df.to_pickle(clusters_bookkeeping_df_pkl_file)

    with open(log_file, 'w') as file:
        file.write("create_clusters_bookkeeping_df: ")
        file.write(str(time.time() - start_time) + "\n")

    return clusters_bookkeeping_df

def main():
    parse_args()
    define_globals()

    clusters_bookkeeping_df = create_clusters_bookkeeping_df()

if __name__ == "__main__":
    main()
