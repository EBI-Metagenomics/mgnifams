import argparse
from pyhmmer import easel
from pyhmmer.plan7 import HMMFile, Pipeline

def main(hmm_file_path, sequence_file_path, out_file_path):
    with HMMFile(hmm_file_path) as hmm_file:
        hmm = next(hmm_file)

    alphabet = easel.Alphabet.amino()
    sequence_db = easel.SequenceFile(sequence_file_path, digital=True)

    pipeline = Pipeline(alphabet)
    hits = pipeline.search_hmm(hmm, sequence_db)

    with open(out_file_path, "wb") as f:
        hits.write(f, format="domains")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Search HMMER database.')
    parser.add_argument('hmm_file_path', type=str, help='Path to the HMM file')
    parser.add_argument('sequence_file_path', type=str, help='Path to the sequence file')
    parser.add_argument('out_file_path', type=str, help='Path to output hits')
    args = parser.parse_args()

    main(args.hmm_file_path, args.sequence_file_path, args.out_file_path)
