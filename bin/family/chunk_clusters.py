import sys
import os

def parse_args():
    global arg_clusters, arg_checked_clusters, \
        arg_minimum_members, arg_num_cluster_chunks
        
    if not (len(sys.argv) == 5):
        print("Incorrect number of args.")
        sys.exit(1)

    arg_clusters           = sys.argv[1]
    arg_checked_clusters   = sys.argv[2]
    arg_minimum_members    = int(sys.argv[3])
    arg_num_cluster_chunks = int(sys.argv[4])

def main():
    parse_args()

    folder_name = 'cluster_chunks'
    if not os.path.exists(folder_name):
        os.makedirs(folder_name)

    for i in range(1, arg_num_cluster_chunks + 1):
        file_path = os.path.join(folder_name, f'{i}.txt')
        with open(file_path, 'w') as file:
            file.write(f'This is file number {i}\n')

if __name__ == "__main__":
    main()
