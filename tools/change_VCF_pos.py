#!/home/bic/anaconda3/envs/bioconda/bin/python
#
# Chongwei 20191027
# Email: chongwei.bi@kaust.edu.sa

import argparse
import os
import re
from operator import itemgetter


def get_argparse():
    parser = argparse.ArgumentParser(description='This is for correcting the position of SNP in vcf file')
    parser.add_argument('-s', '--save_path', type=str, default="./", help='path/to/save/')
    parser.add_argument('-v', '--vcf_file', type=validate_file, required=True, help='path/to/file.vcf')
    parser.add_argument('-c', '--chr_name', type=str, required=True, help='chromosome name')
    parser.add_argument('-p', '--pos_change', type=str, required=True, help='position change')

    return parser.parse_args()


def validate_file(x):
    if not os.path.exists(x):
        raise argparse.ArgumentTypeError("{0} does not exist".format(x))
    return x


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

                if pos_change_by is "+":
                    pos = pos + pos_change_value
                    outfile1.writelines('\t'.join((chr_name, str(pos), '\t'.join(line[2:]))))
                    line = tuple([chr_name, pos] + line[2:])
                    new_line.append(line)

                elif pos_change_by is "-":
                    pos = pos - pos_change_value
                    outfile1.writelines('\t'.join((chr_name, str(pos), '\t'.join(line[2:]))))
                    line = tuple([chr_name, pos] + line[2:])
                    new_line.append(line)

                else:
                    print("\npos change is not +/-\n")

        new_line.sort(key=itemgetter(1), reverse=False)
        for record in new_line:
            outfile2.writelines('\t'.join([record[0], str(record[1])] + list(record[2:])))


def main():
    args = get_argparse()
    if not os.path.isdir(args.save_path):
        os.mkdir(args.save_path)
    print("\n-------------- received arguments -------------")
    print("vcf file  :   %s" % args.vcf_file)
    print("save path :   %s" % args.save_path)
    print("chr name  :   %s" % args.chr_name)
    print("pos change:   %s" % args.pos_change)

    correct_vcf(args)


if __name__ == '__main__':
    main()