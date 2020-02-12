#!/home/bic/anaconda3/envs/bioconda/bin/python
#
# Chongwei 20191029
# Email: chongwei.bi@kaust.edu.sa

import argparse
import os
import re
from operator import itemgetter


def get_argparse():
    parser = argparse.ArgumentParser(description='This is for filtering sv.vcf by the proportion of reads supporting SVs')
    parser.add_argument('-s', '--save_path', type=str, default="./", help='path/to/save/')
    parser.add_argument('-v', '--vcf_file', type=validate_file, required=True, help='path/to/file.vcf')
    parser.add_argument('-f', '--sv_percent', type=float, default="0.5", help='percentage of sv reads')

    return parser.parse_args()


def validate_file(x):
    if not os.path.exists(x):
        raise argparse.ArgumentTypeError("{0} does not exist".format(x))
    return x


def filter_sv(vcf_file, save_path, sv_percent):
    vcf = vcf_file
    file_name = re.split(r'/', vcf)[-1]
    file_name = re.split(r'.vcf', file_name)[0]
    outvcf1 = save_path + '/' + file_name + '.filtered.' + str(sv_percent) + '.vcf'
    outvcf2 = save_path + '/' + file_name + '.filtered.' + str(sv_percent) + '.sorted.vcf'

    with open(vcf, 'r') as infile, open(outvcf1, 'w') as outfile1, open(outvcf2, 'w') as outfile2:
        new_line = []
        for line in infile:
            if line.startswith('#'):
                outfile1.writelines(line)
                outfile2.writelines(line)
            else:
                line = re.split(r'\t', line)
                pos = int(line[1])
                umi_name = line[2]
                read_count = int(re.split(r'_', umi_name)[0])
                info = line[7]
                sv_read_count = int(re.split(r'RE=', info)[1])

                if sv_read_count/read_count >= sv_percent:
                    outfile1.writelines('\t'.join(line))
                    line = tuple([line[0], pos] + line[2:])
                    new_line.append(line)

        new_line.sort(key=itemgetter(1), reverse=False)
        for record in new_line:
            outfile2.writelines('\t'.join([record[0], str(record[1])] + list(record[2:])))


def main():
    args = get_argparse()
    if not os.path.isdir(args.save_path):
        os.mkdir(args.save_path)
    print("\n-------------- received arguments -------------")
    print("vcf file      :   %s" % args.vcf_file)
    print("save path     :   %s" % args.save_path)
    print("SV percentage :   %s" % args.sv_percent)

    filter_sv(args.vcf_file, args.save_path, args.sv_percent)


if __name__ == '__main__':
    main()