#!/usr/bin/env python3
# Chongwei 20191027
# bicwei@gmail.com

import argparse
import os
import re
import sys
from operator import itemgetter


def get_argparse():
    parser = argparse.ArgumentParser(description='This is for correcting the position of SNP in vcf file')
    parser.add_argument('-v', '--vcf_file', type=validate_file, required=True, help='path/to/file.vcf')
    parser.add_argument('-c', '--chr_name', type=str, required=True, help='chromosome name')
    parser.add_argument('-p', '--pos_change', type=str, required=True, help='position change')
    parser.add_argument('-b', '--total_base', type=int, help='the length of reference genome.'
                                                             'if value is not None, will reverse the coordinate')
    parser.add_argument('-s', '--save_path', type=str, default="./", help='path/to/save/')
    
    return parser.parse_args()


def validate_file(x):
    if not os.path.exists(x):
        raise argparse.ArgumentTypeError("{0} does not exist".format(x))
    return x


def reverse_coordinate(args):
    vcf = args.vcf_file
    file_name = re.split(r'/', vcf)[-1]
    file_name = re.split(r'.vcf', file_name)[0]
    outvcf = args.save_path + '/' + file_name + '_PosBaseRev.vcf.indel_pos_not_right.tmp'
    total_base = int(args.total_base) + 1
    with open(vcf, 'r') as infile, open(outvcf, 'w') as outfile:
        for line in infile:
            if line.startswith('#'):
                outfile.writelines(line)
            else:
                line = re.split(r'\t', line)
                pos = total_base - int(line[1])
                REF_base = line[3]
                ALT_base = line[4]
                complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
                REF_reverse_complement = "".join(complement.get(base, base) for base in reversed(REF_base))
                ALT_reverse_complement = "".join(complement.get(base, base) for base in reversed(ALT_base))
                outfile.writelines('\t'.join((line[0], str(pos), line[2], str(REF_reverse_complement),
                                              str(ALT_reverse_complement), '\t'.join(line[5:]))))
    return outvcf


def correct_vcf(args):
    chr_name = args.chr_name
    vcf = args.vcf_file
    file_name = re.split(r'/', vcf)[-1]
    file_name = re.split(r'.vcf', file_name)[0]
    outvcf1 = args.save_path + '/' + file_name + '_correct.pos.vcf'
    outvcf2 = args.save_path + '/' + file_name + '_correct.pos.sorted.vcf'
    pos_change_value = int(args.pos_change[1:])
    pos_change_by = args.pos_change[0]
    with open(vcf, 'r') as infile, open(outvcf1, 'w') as outfile1, open(outvcf2, 'w') as outfile2:
        new_line = []
        for line in infile:
            if line.startswith('#'):
                outfile1.writelines(line)
                outfile2.writelines(line)
            else:
                line = re.split(r'\t', line)
                pos = int(line[1])

                if pos_change_by == "+":
                    pos = pos + pos_change_value
                    outfile1.writelines('\t'.join((chr_name, str(pos), '\t'.join(line[2:]))))
                    line = tuple([chr_name, pos] + line[2:])
                    new_line.append(line)

                elif pos_change_by == "-":
                    pos = pos - pos_change_value
                    outfile1.writelines('\t'.join((chr_name, str(pos), '\t'.join(line[2:]))))
                    line = tuple([chr_name, pos] + line[2:])
                    new_line.append(line)

                else:
                    print("\npos change is not +/-\n")

        new_line.sort(key=itemgetter(1), reverse=False)
        for record in new_line:
            outfile2.writelines('\t'.join([record[0], str(record[1])] + list(record[2:])))


def position_main(args):
    if not os.path.isdir(args.save_path):
        os.mkdir(args.save_path)
    print("\n-------------- received arguments -------------")
    print("vcf file  :   %s" % args.vcf_file)
    print("save path :   %s" % args.save_path)
    print("chr name  :   %s" % args.chr_name)
    print("pos change:   %s" % args.pos_change)

    if args.total_base is not None:
        print("\nThe coordinate will be reversed, to make the strand direction same as published reference\n")
        args.vcf_file = reverse_coordinate(args)

    correct_vcf(args)
    print("\nchange VCF file pos and chr done!\n")
    sys.exit(0)


if __name__ == '__main__':
    args = get_argparse()
    position_main(args)
