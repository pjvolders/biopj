#!/usr/bin/env python3

# import modules used here
import sys
import argparse
import requests, sys
import json
from tqdm import tqdm


def get_seq(chrom, start, end, strand):
    if chrom[0:3] == "chr":
        chrom = chrom[3:]
    if strand == "+":
        strand = 1
    else:
        strand = -1
    server = "https://rest.ensembl.org"
    ext = "/sequence/region/homo_sapiens/{}:{}..{}:{}".format(chrom, start, end, strand)

    r = requests.get(server + ext, headers={"Content-Type": "application/json"})

    if not r.ok:
        r.raise_for_status()
        sys.exit()

    decoded = r.json()

    return decoded['seq']


def get_circ_seq(chrom, bsj_start, bsj_end, strand, span=70):
    width = span // 2
    last_seq = get_seq(chrom, bsj_start, bsj_start + width - 1, strand)
    first_seq = get_seq(chrom, bsj_end - width + 1, bsj_end, strand)
    if strand == "+":
        return first_seq + last_seq
    return last_seq + first_seq


# Gather our code in a main() function
def main():
    parser = argparse.ArgumentParser(
        description='Birthday Python program for Marieke! Obtain circRNA backsplice sequence from Ensembl.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('input_file', metavar='INPUT', type=str,
                        help='The input file in BED format')
    parser.add_argument('-s', dest='span', type=int, default=70,
                        help='The span around the BSJ we want to obtain the sequence')
    parser.add_argument('-f', dest='fasta', action='store_true',
                        help='The output should be in FASTA format')

    args = parser.parse_args()
    # print(args.accumulate(args.integers))

    with open(args.input_file) as f:
        for line in tqdm(f):
            cols = line.rstrip().split("\t")
            chrom, start, end, name, score, strand = cols[0:6]
            # print((chrom, start, end, name, score, strand))
            if args.fasta:
                print(">{}".format(name))
            print(get_circ_seq(chrom, int(start) + 1, int(end), strand, span=args.span))

    f.close()


# Standard boilerplate to call the main() function to begin
# the program. Program can be used as module.
if __name__ == '__main__':
    main()
