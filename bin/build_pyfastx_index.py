#!/usr/bin/env python3

import argparse
import pyfastx


def parse_args(args=None):
    parser = argparse.ArgumentParser(description="Create pyfastx index file.")
    parser.add_argument("-f", "--fasta_file", required=True, type=str, help="Path to amino acid sequence FASTA file.")

    return parser.parse_args(args)


def main():
    args = parse_args()
    
    pyfastx.Fasta(args.fasta_file)


if __name__ == "__main__":
    main()
