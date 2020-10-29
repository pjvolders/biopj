#!/usr/bin/env python

# import modules used here
import sys
import argparse
from Bio import SeqIO
import re


def find_splits(sequence, stop_codon_char="*") -> list:
    stop_codons = [m.start()+1 for m in re.finditer(re.escape(stop_codon_char), sequence)]
    stop_codons = [0] + stop_codons
    if stop_codons[-1] != len(sequence):
        stop_codons.append(len(sequence))
    return stop_codons


def split_sequence(sequence, stop_codon_char="*", min_length=1):
    stop_codons = find_splits(sequence, stop_codon_char)
    seq_starts = stop_codons[:-1]
    seq_ends = stop_codons[1:]

    for i in range(len(seq_starts)):
        s = seq_starts[i]
        e = seq_ends[i]
        subseq = sequence[s:e]
        if subseq[-1] == stop_codon_char:
            subseq = subseq[:-1]
        if len(subseq) < min_length:
            continue
        yield subseq, s, e


# Gather our code in a main() function
def main():
    parser = argparse.ArgumentParser(description='Split fasta file on stop codons.',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('input_file', metavar='INPUT', type=str,
                        help='The input file')
    parser.add_argument(
        "-l", "--min-length",
        metavar="<min-length>",
        action="store",
        dest="min_length",
        default=7,
        type=int,
        help="minimum peptide length to include",
    )

    args = parser.parse_args()
    # print(args.accumulate(args.integers))

    fasta_sequences = SeqIO.parse(open(args.input_file), 'fasta')

    for fasta in fasta_sequences:
        name, seq = fasta.id, str(fasta.seq)

        for sequence, start, end in split_sequence(seq, min_length=args.min_length):

            fasta_header = ">{}_{}-{}".format(name, start, end)
            print(fasta_header)
            print(sequence)


# Standard boilerplate to call the main() function to begin
# the program. Program can be used as module.
if __name__ == '__main__':
    main()
