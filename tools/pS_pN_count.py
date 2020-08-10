#!/usr/bin/env python3
# Chongwei 20200725
# Email: chongwei.bi@kaust.edu.sa

import sys
import argparse
import os
import logging


def get_argparse():
    parser = argparse.ArgumentParser(description='calculate synonymous and nonsynonymous position of input sequence')
    parser.add_argument('-c', '--codon_table', type=validate_file, required=True,
                        help='path/to/codon_table format as ATG M\nATG M\nTAA stop')
    parser.add_argument('-d', '--dna_sequence', type=str,
                        help='DNA sequence in ATGC, should be only the coding sequence with length of multiple of 3')
    parser.add_argument('-n', '--sequence_name', default="unknown", help='gene name of input DNA sequence')
    parser.add_argument('-o', '--out_file', type=str, required=True,
                        help='write results to this file, please note this will add new result to existing file')
    args = parser.parse_args()
    return args


def validate_file(x):
    if not os.path.exists(x):
        raise argparse.ArgumentTypeError("{0} does not exist".format(x))
    return x


def codon_table(args):
    # codon_file = "/Users/bic/Desktop/vertebrate_mitochondrial_code.txt"
    # dna_sequence = "cccgggttt"
    codon_file = args.codon_table
    dna_sequence = args.dna_sequence

    codon_dir = {}
    with open(codon_file, "r") as infile:
        for line in infile:
            base, aa = line.strip().split()
            codon_dir[base.strip()] = aa.strip()

    # check if input sequence is valid
    dna_sequence=dna_sequence.upper()

    if len(dna_sequence) % 3 != 0:
        sys.stderr.write('the length of input DNA sequence should be multiple of 3')
        sys.exit(1)

    for letter in dna_sequence:
        if letter not in "ATGC":
            sys.stderr.write('input DNA sequence contains letter not [ATGC]')
            sys.exit(1)
    return codon_dir, dna_sequence


def calculate_pS_pN(codon_dir, dna_sequence, args):
    # calculate synonymous and nonsynonymous position of input sequence
    # split dna_sequence to 3bp list
    n = 3
    aa_sequence = [dna_sequence[i:i+n] for i in range(0, len(dna_sequence), n)]
    nonsynonymous_site = 0

    for seq_3bp in aa_sequence:
        seq_aa = codon_dir[seq_3bp]
        seq_3bp_list = [letter for letter in seq_3bp]

        for number in [0, 1, 2]:
            base_list = ['A', 'G', 'T', 'C']
            base_list.remove(seq_3bp[number])
            seq_3bp_list_used = list(seq_3bp_list)  # or use seq_3bp_list.copy() to make a new list

            for base in base_list:
                seq_3bp_list_used[number] = base
                current_3bp = "".join(seq_3bp_list_used)
                if codon_dir[current_3bp] != seq_aa:
                    nonsynonymous_site += 0.3333333
                # else:
                #     print(current_3bp)
                #     print(codon_dir[current_3bp])

    logging.info(" ".join([args.sequence_name, str(len(dna_sequence)), str(nonsynonymous_site),
                           str(len(dna_sequence)-nonsynonymous_site), str(dna_sequence)]))


def main():
    args = get_argparse()
    print("\nformat: name length nonsynonymous_site synonymous_site")

    log_file = args.out_file
    logging.basicConfig(level=logging.DEBUG,
                        format='%(message)s',
                        filename=log_file,
                        filemode='a')
    logging.getLogger().addHandler(logging.StreamHandler(sys.stdout))

    codon_dir, dna_sequence = codon_table(args)
    calculate_pS_pN(codon_dir, dna_sequence, args)


if __name__ == '__main__':
    main()
