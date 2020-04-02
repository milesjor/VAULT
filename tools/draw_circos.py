#!/home/bic/anaconda3/envs/bioconda/bin/python
#
# Chongwei 20190930
# Email: chongwei.bi@kaust.edu.sa

import numpy as np
import pandas as pd
import argparse
import os
import re
import io
import sys
from operator import itemgetter


def get_argparse():
    parser = argparse.ArgumentParser(description='This is for preparing data for circos. For vcf file, it will correct position and remove indel.'
                                                 'For depth file, it will bin the depeth file based on user defined bin size')
    parser.add_argument('-n', '--name', type=str, default="circos", help='prefix of output file')
    parser.add_argument('-d', '--depth_file', type=validate_file, help='depth file from samtools depth')
    parser.add_argument('-s', '--save_path', type=str, required=True, help='path/to/save/')
    parser.add_argument('-v', '--vcf_file', type=validate_file, help='path/to/file.vcf')
    parser.add_argument('-b', '--bin_size', type=int, default="30", help='how many bases per bin')
    parser.add_argument('-g', '--genome', type=str, default='hchrM.F9.UMIs',
                     choices=['mmt_nod_F6_10N_to_C57', 'hchrM.F9.UMIs'], help='the genome used in VAULT')
    parser.add_argument('-c', '--chr_name', type=str, help='chromosome name showed in vcf file')

    args = parser.parse_args()
    if args.genome == 'mmt_nod_F6_10N_to_C57':
        args.chr_name = 'chrMT'
    elif args.genome == 'hchrM.F9.UMIs':
        args.chr_name = 'chrM'

    return args


def validate_file(x):
    if not os.path.exists(x):
        raise argparse.ArgumentTypeError("{0} does not exist".format(x))
    return x


def point2bin(args):
    depth_file1 = pd.read_csv(args.depth_file, header=None, sep='\t')
    bin_size = args.bin_size
    depth = depth_file1[2]
    depth_np = depth.to_numpy()
    length = len(depth_np)
    bin_depth = np.mean(depth_np[:(length//bin_size)*bin_size].reshape(-1, bin_size), axis=1)
    rest = depth_np[(length//bin_size)*bin_size:]

    if rest.size != 0:
        bin_depth = np.append(bin_depth, rest.mean())

    column4 = bin_depth.astype(int)
    column2 = np.arange(1, length, bin_size)
    column3 = np.arange(bin_size, length, bin_size)
    column1 = np.chararray(len(column4), itemsize=3)
    column1[:] = 'mt1'

    if rest.size != 0:
        column3 = np.append(column3, length)

    combine = np.column_stack((column1, column2, column3, column4)).astype(str)
    save_file = args.save_path + '/' + args.name + '_' + str(bin_size) + 'bin.depth.txt'
    np.savetxt(save_file, combine, fmt='%s')


def correct_vcf(args):
    chr_name = args.chr_name
    vcf = args.vcf_file
    outvcf1 = args.save_path + '/' + args.name + '_correct.pos.vcf'
    outvcf2 = args.save_path + '/' + args.name + '_correct.pos.sorted.vcf'
    genome_used = args.genome
    with open(vcf, 'r') as infile, open(outvcf1, 'w') as outfile1, open(outvcf2, 'w') as outfile2:
        new_line = []
        for line in infile:
            if line.startswith('#'):
                outfile1.writelines(line)
                outfile2.writelines(line)
            else:
                line = re.split(r'\t', line)
                pos = int(line[1])
                info = line[7]
                # if not info.startswith('INDEL'):  # below is for mmt_nod_F6_10N.fa to mouse C57 genome
                if genome_used == 'mmt_nod_F6_10N_to_C57':
                    if 41 <= pos <= 1599:
                        pos = pos + 14700
                        outfile1.writelines('\t'.join((chr_name, str(pos), '\t'.join(line[2:]))))
                        line = tuple([chr_name, pos] + line[2:])
                        new_line.append(line)

                    elif 1599 < pos <= 11427:
                        pos = pos - 1599
                        outfile1.writelines('\t'.join((chr_name, str(pos), '\t'.join(line[2:]))))
                        line = tuple([chr_name, pos] + line[2:])
                        new_line.append(line)

                    elif pos > 11427:
                        pos = pos - 1601
                        outfile1.writelines('\t'.join((chr_name, str(pos), '\t'.join(line[2:]))))
                        line = tuple([chr_name, pos] + line[2:])
                        new_line.append(line)

                    else:
                        print("\nWrong position: \n%s\n" % '\t'.join(line))

                elif genome_used == 'hchrM.F9.UMIs':
                    if 38 <= pos <= 13307:
                        pos = pos + 3262
                        outfile1.writelines('\t'.join((chr_name, str(pos), '\t'.join(line[2:]))))
                        line = tuple([chr_name, pos] + line[2:])
                        new_line.append(line)

                    elif 13308 < pos <= 16571:
                        pos = pos - 13307
                        outfile1.writelines('\t'.join((chr_name, str(pos), '\t'.join(line[2:]))))
                        line = tuple([chr_name, pos] + line[2:])
                        new_line.append(line)

                    else:
                        print("\nWrong position: \n%s\n" % '\t'.join(line))

                else:
                    sys.stderr.write('There is no preset parameter for changing chromosome position!')
                    sys.exit(1)

        new_line.sort(key=itemgetter(1), reverse=False)
        for record in new_line:
            outfile2.writelines('\t'.join([record[0], str(record[1])] + list(record[2:])))


def vcf_circos(args):
    vcf = args.save_path + '/' + args.name + '_correct.pos.sorted.vcf'
    out = args.save_path + '/' + args.name + '_correct.pos.snp.txt'

    def read_vcf(path):
        with open(path, 'r') as f:
            lines = [l for l in f if not l.startswith('##')]
        return pd.read_csv(
            io.StringIO(''.join(lines)),
            dtype={'#CHROM': str, 'POS': int, 'ID': str, 'REF': str, 'ALT': str,
                   'QUAL': str, 'FILTER': str, 'INFO': str},
            sep='\t'
        ).rename(columns={'#CHROM': 'CHROM'})

    pd_vcf = read_vcf(vcf)

    # print(pd_vcf['POS'])
    # print(pd_vcf.POS.nunique())
    # print(pd_vcf.POS.count())
    # print(pd_vcf.POS.value_counts())
    # print(pd_vcf.POS.value_counts().reset_index())

    count_pos = pd_vcf.POS.value_counts().reset_index().rename(columns={'POS': 'count', 'index': 'POS'})
    sort_count = count_pos.sort_values('POS').reset_index(drop=True)
    sort_count.insert(0, 'chr', 'mt1')
    snp_df = sort_count[['chr', 'POS', 'POS', 'count']]

    np.savetxt(out, snp_df, fmt='%s')


def main():
    args = get_argparse()
    if not os.path.isdir(args.save_path):
        os.mkdir(args.save_path)
    print("\n-------------- received arguments -------------")
    print("file name:   %s" % args.name)
    print("save path:   %s" % args.save_path)

    if os.path.exists(str(args.depth_file)):
        print("depth file:  %s" % args.depth_file)
        point2bin(args)

    if os.path.exists(str(args.vcf_file)):
        print("vcf file:    %s\n" % args.vcf_file)
        correct_vcf(args)
        vcf_circos(args)


if __name__ == '__main__':
    main()