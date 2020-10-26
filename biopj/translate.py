#!/usr/bin/env python3

# import modules used here
import sys
import argparse
from Bio import SeqIO


# Gather our code in a main() function
def main():
    parser = argparse.ArgumentParser(description='Translate a FASTA file.',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('input_file', metavar='INPUT', type=str,
                        help='The input file in FASTA format')
    parser.add_argument('output_file', metavar='OUTPUT', type=str,
                        help='The output file in FASTA format')
    parser.add_argument('--all', dest='all_reading_frames', action='store_true',
                        help='Translate in all reading frames, default only the first')

    args = parser.parse_args()

    with open(args.output_file, 'w') as o_f:
        with open(args.input_file) as f:
            for record in SeqIO.parse(f, "fasta"):
                sequence = record.seq
                id = record.id

                if args.all_reading_frames:
                    for rf in range(3):
                        sequence_rf = sequence[rf:(len(sequence) - (len(sequence)-rf) % 3)]  # trim to multiple of 3
                        record.seq = sequence_rf.translate()
                        record.id = "{}_RF{}".format(id, rf+1)
                        SeqIO.write(record, o_f, "fasta")
                else:
                    sequence = sequence[0:(len(sequence)-len(sequence) % 3)]  # trim to multiple of 3
                    record.seq = sequence.translate()
                    SeqIO.write(record, o_f, "fasta")

        f.close()
    o_f.close()


# Standard boilerplate to call the main() function to begin
# the program. Program can be used as module.
if __name__ == '__main__':
    main()
