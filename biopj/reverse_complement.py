#!/usr/bin/env python3

# import modules used here
import sys
import argparse
from Bio import SeqIO


# Gather our code in a main() function
def main():
    parser = argparse.ArgumentParser(description='Reverse complement a FASTA file.', 
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('input_file', metavar='INPUT', type=str,
                                         help='The input file in FASTA format')
    parser.add_argument('output_file', metavar='OUTPUT', type=str,
                                         help='The output file in FASTA format')

    args = parser.parse_args()
    #print(args.accumulate(args.integers))


    with open(args.output_file, 'w') as o_f:
        with open(args.input_file) as f:
            for record in SeqIO.parse(f, "fasta"):
                record.seq = record.seq.reverse_complement()
                SeqIO.write(record, o_f, "fasta")

        f.close()
    o_f.close()

# Standard boilerplate to call the main() function to begin
# the program. Program can be used as module.
if __name__ == '__main__':
    main()
