#!/usr/bin/env python3

# import modules used here
import argparse
import requests
import sys


def get_pop_freqs(snp_id, population='gnomADe:ALL'):
    server = "https://rest.ensembl.org"
    ext = "/variation/human/{}?pops=1".format(snp_id)

    r = requests.get(server + ext, headers={"Content-Type": "application/json"})

    if not r.ok:
        r.raise_for_status()
        sys.exit()

    decoded = r.json()
    populations = decoded['populations']

    return {p['allele']: p['frequency'] for p in populations if p['population'] == population}

# Gather our code in a main() function
def main():
    parser = argparse.ArgumentParser(
        description='Gets the population frequencies for SNPs based on the rsID',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('snp_id', metavar='SNP_ID', type=str, nargs='+',
                        help='The input file in BED format')
    parser.add_argument('-p', dest='population', type=str, default="gnomADe:ALL",
                        help='The gnomAD population ID')
    parser.add_argument('-t', dest='table', action='store_true',
                        help='The output should be in a tabular format instead of JSON. \nThe collumns for the table '
                             'are NAME, MAJOR_ALLELE, MINOR_ALLELE, MAF')

    args = parser.parse_args()
    if args.table:
        print("NAME\tMAJOR_ALLELE\tMINOR_ALLELE\tMAF")
    for rs_id in args.snp_id:
        try:
            freqs = get_pop_freqs(rs_id, population=args.population)
        except requests.exceptions.HTTPError:
            freqs = {}

        if args.table:
            if len(freqs) == 0:
                print("{}\t-\t-\t-".format(rs_id))
                continue
            freqs = [{'allele': a, 'frequency': freqs[a]} for a in freqs.keys()]
            freqs.sort(key=lambda i: i['frequency'], reverse=True)
            freqs = "{}\t{}\t{}\t{}".format(rs_id, freqs[0]['allele'], freqs[1]['allele'], freqs[1]['frequency'])
        print(freqs)
    # print(args.accumulate(args.integers))



# Standard boilerplate to call the main() function to begin
# the program. Program can be used as module.
if __name__ == '__main__':
    main()
