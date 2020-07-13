#!/usr/bin/env python3
# Chongwei 20190930
# bicwei@gmail.com

import numpy as np
import pandas as pd
import argparse
import os
import re
import io
import sys
import subprocess
from operator import itemgetter


def get_argparse():
    parser = argparse.ArgumentParser(description='This is for preparing data for circos. For vcf file, it will correct position and remove indel.'
                                                 'For depth file, it will bin the depeth file based on user defined bin size')
    parser.add_argument('-n', '--name', type=str, default="circos", help='prefix of output file')
    parser.add_argument('-d', '--depth_file', type=validate_file, help='depth file from samtools depth')
    parser.add_argument('-s', '--save_path', type=str, required=True, help='path/to/save/')
    parser.add_argument('-v', '--vcf_file', type=validate_file, help='path/to/file.vcf')
    parser.add_argument('-b', '--bin_size', type=int, default="30", help='how many bases per bin')
    parser.add_argument('-g', '--genome', type=str,
                        choices=['mmt_nod_F6_10N_to_C57', 'hchrM.F9.UMIs', 'hchrM.F6.UMIs', 'mmt_c57_F6_10N'],
                        help='the genome used in VAULT')
    parser.add_argument('-c', '--chr_name', type=str, help='chromosome name showed in vcf file')
    parser.add_argument('-A', '--keep_1alt', action='store_true',
                        help='leave only 1 Alt in vcf file for helping SNV annotation')
    parser.add_argument('--depth_pos', type=validate_file, help='change the position of depth file')

    args = parser.parse_args()
    if args.genome == 'mmt_nod_F6_10N_to_C57':
        args.chr_name = 'chrMT'
    elif args.genome == 'hchrM.F9.UMIs':
        args.chr_name = 'chrM'
    elif args.genome == 'hchrM.F6.UMIs':
        args.chr_name = 'chrM'
    elif args.genome == 'mmt_c57_F6_10N':
        args.chr_name = 'chrMT'

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
            # info_field = re.split(r'\t', line)[7]
            if line.startswith('#'):
                outfile1.writelines(line)
                outfile2.writelines(line)

            # elif info_field.startswith('INDEL'):  # skip indel in vcf
            #     pass

            else:
                line = re.split(r'\t', line)
                pos = int(line[1])
                # below is for mmt_nod_F6_10N.fa to mouse C57 genome
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
                        print("WARNING!  Wrong position: %s" % '\t'.join(line).strip())

                elif genome_used == 'hchrM.F9.UMIs':
                    if 38 <= pos <= 13307:
                        pos = pos + 3262
                        outfile1.writelines('\t'.join((chr_name, str(pos), '\t'.join(line[2:]))))
                        line = tuple([chr_name, pos] + line[2:])
                        new_line.append(line)

                    elif 13307 < pos <= 16571:
                        pos = pos - 13307
                        outfile1.writelines('\t'.join((chr_name, str(pos), '\t'.join(line[2:]))))
                        line = tuple([chr_name, pos] + line[2:])
                        new_line.append(line)

                    else:
                        print("WARNING!  Wrong position: %s" % '\t'.join(line).strip())

                elif genome_used == 'hchrM.F6.UMIs':
                    if 39 <= pos <= 1748:
                        pos = pos + 14821
                        outfile1.writelines('\t'.join((chr_name, str(pos), '\t'.join(line[2:]))))
                        line = tuple([chr_name, pos] + line[2:])
                        new_line.append(line)

                    elif 1748 < pos <= 16573:
                        pos = pos - 1748
                        outfile1.writelines('\t'.join((chr_name, str(pos), '\t'.join(line[2:]))))
                        line = tuple([chr_name, pos] + line[2:])
                        new_line.append(line)

                    else:
                        print("WARNING!  Wrong position: %s" % '\t'.join(line).strip())

                elif genome_used == 'mmt_c57_F6_10N':
                    if 41 <= pos <= 1599:
                        pos = pos + 14700
                        outfile1.writelines('\t'.join((chr_name, str(pos), '\t'.join(line[2:]))))
                        line = tuple([chr_name, pos] + line[2:])
                        new_line.append(line)

                    elif 1599 < pos <= 16296:
                        pos = pos - 1599
                        outfile1.writelines('\t'.join((chr_name, str(pos), '\t'.join(line[2:]))))
                        line = tuple([chr_name, pos] + line[2:])
                        new_line.append(line)

                    else:
                        print("WARNING!  Wrong position: %s" % '\t'.join(line).strip())

                else:
                    sys.stderr.write('ERROR!  There is no preset parameter for changing chromosome position!')
                    sys.exit(1)

        new_line.sort(key=itemgetter(1), reverse=False)
        for record in new_line:
            outfile2.writelines('\t'.join([record[0], str(record[1])] + list(record[2:])))


def correct_depth(args):
    depth_file = args.depth_pos
    file_name = re.split(r'/', depth_file)[-1]
    file_name = re.split(r'.txt', file_name)[0]
    save_file = args.save_path + '/' + file_name + '.correct.txt'
    save_file2 = args.save_path + '/' + file_name + '.correct.sorted.txt'
    genome_used = args.genome

    with open(depth_file, 'r') as infile, open(save_file, 'w') as outfile1, open(save_file2, 'w') as outfile2:
        new_line = []

        def write_line(line, pos):
            outfile1.writelines('\t'.join((line[0], str(pos), str(line[2]))))
            line = tuple([line[0], pos, line[2]])
            new_line.append(line)

        for line in infile:
            line = re.split(r'\t', line)
            pos = int(line[1])
            # below is for mmt_nod_F6_10N.fa to mouse C57 genome
            if genome_used == 'mmt_nod_F6_10N_to_C57':
                if 41 <= pos <= 1599:
                    pos = pos + 14700
                    write_line(line, pos)

                elif 1599 < pos <= 11427:
                    pos = pos - 1599
                    write_line(line, pos)

                elif pos > 11427:
                    pos = pos - 1601
                    write_line(line, pos)

                else:
                    print("WARNING!  Wrong position: %s" % '\t'.join(line).strip())

            elif genome_used == 'hchrM.F9.UMIs':
                if 38 <= pos <= 13307:
                    pos = pos + 3262
                    write_line(line, pos)

                elif pos > 13307:
                    pos = pos - 13307
                    write_line(line, pos)

                else:
                    print("WARNING!  Wrong position: %s" % '\t'.join(line).strip())

            elif genome_used == 'hchrM.F6.UMIs':
                if 39 <= pos <= 1748:
                    pos = pos + 14821
                    write_line(line, pos)

                elif pos > 1748:
                    pos = pos - 1748
                    write_line(line, pos)

                else:
                    print("WARNING!  Wrong position: %s" % '\t'.join(line).strip())

            elif genome_used == 'mmt_c57_F6_10N':
                if 41 <= pos <= 1599:
                    pos = pos + 14700
                    write_line(line, pos)

                elif pos > 1599:
                    pos = pos - 1599
                    write_line(line, pos)

                else:
                    print("WARNING!  Wrong position: %s" % '\t'.join(line).strip())

            else:
                sys.stderr.write('ERROR!  There is no preset parameter for changing chromosome position!')
                sys.exit(1)

        new_line.sort(key=itemgetter(1), reverse=False)
        for record in new_line:
            outfile2.writelines('\t'.join([record[0], str(record[1]), str(record[2])]))


def vcf_circos(args):
    vcf = args.save_path + '/' + args.name + '_correct.pos.sorted.vcf'
    out = args.save_path + '/' + args.name + '_correct.pos.snv.txt'
    snv_vcf = args.save_path + '/' + args.name + '_correct.pos.sorted.snv.vcf'

    with open(vcf, 'r') as infile, open(snv_vcf, 'w') as outfile:
        for line in infile:
            if line.startswith('#'):
                outfile.writelines(line)
            else:
                info_field = re.split(r'\t', line)[7]
                if not info_field.startswith('INDEL'):
                    outfile.writelines(line)

    def read_vcf(path):
        with open(path, 'r') as f:
            lines = [l for l in f if not l.startswith('##')]
        return pd.read_csv(
            io.StringIO(''.join(lines)),
            dtype={'#CHROM': str, 'POS': int, 'ID': str, 'REF': str, 'ALT': str,
                   'QUAL': str, 'FILTER': str, 'INFO': str},
            sep='\t'
        ).rename(columns={'#CHROM': 'CHROM'})

    pd_vcf = read_vcf(snv_vcf)

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


def count_to_freq(args):
    position = args.save_path + '/all.coverage.3plus.correct.pos.txt'
    circos = args.save_path + '/' + args.name + '_correct.pos.snv.txt'
    corrected_circos = args.save_path + '/' + args.name + '_correct.pos.snv.freq.txt'

    with open(position, 'r') as infile1, open(circos, 'r') as infile2, open(corrected_circos, 'w') as outfile:
        dir = {}
        for line in infile1:
            v, k = line.strip().split(' ')
            dir[k.strip()] = v.strip()
            # print(dir)

        for line2 in infile2:
            line2_split = line2.split(' ')
            count = line2_split[3]
            pos = line2_split[1]
            depth = dir[pos]
            freq = round(int(count) / int(depth), 6)
            new_line2 = line2_split[0] + ' ' + line2_split[1] + ' ' + line2_split[2] + ' ' + str(freq) + '\n'
            outfile.writelines(new_line2)


def snv_load_per_group(args):
    snv_file = args.save_path + '/' + args.name + '_correct.pos.sorted.snv.vcf'
    coverage_file = args.save_path + '/coverage.3plus.txt'
    snv_load_distribution = args.save_path + '/' + args.name + '_snv_load_distribution.txt'

    with open(coverage_file, 'r') as infile1, open(snv_load_distribution, 'w') as outfile:
        dir = {}
        for line in infile1:
            k, v = line.strip().split(' ')
            dir[k.strip()] = v.strip()

        awk = """awk '{print $1,$2}' """
        snv_count = args.save_path + '/snv.count.per.group.txt'
        cmd = """cat {snv_file} | grep -v "^#" | cut -f3 | sort | uniq -c | {awk} > {snv_count} 
                """.format(snv_file=snv_file, awk=awk, snv_count=snv_count)
        subprocess.run(cmd, shell=True, check=True, stdout=subprocess.PIPE)

        with open(snv_count, 'r') as infile2:
            for line2 in infile2:
                line2_split = line2.split(' ')
                group_name = line2_split[1].strip()
                snv_number = line2_split[0].strip()
                coverage = dir[group_name]
                load = round(int(snv_number) / int(coverage), 6)
                new_line = " ".join((group_name, str(snv_number), str(coverage), str(load))) + '\n'
                outfile.writelines(new_line)


def rem_mlt_alt(args):
    vcf = args.vcf_file
    file_name = re.split(r'/', vcf)[-1]
    file_name = re.split(r'.vcf', file_name)[0]
    outvcf = args.save_path + '/' + file_name + '.1alt.vcf'
    with open(vcf, 'r') as infile, open(outvcf, 'w') as outfile:
        for line in infile:
            if line.startswith('#'):
                # pass
                outfile.writelines(line)
            else:
                # oldline = line
                line = re.split(r'\t', line)
                alt = line[4]
                info = line[7]
                last_field = line[9]
                pl = re.split(r':', last_field)[1]
                gt_pl = '1:' + str(re.split(r',', pl)[0]) + ',0\n'

                if ',' in alt:
                    alt = re.split(r',', alt)[0]
                    # print(info)
                    info = re.split(r'AC=|;AN=', info)
                    info_1 = info[0]
                    info_2 = info[2:]
                    ac = info[1].split(',')[0]
                    info = str(info_1) + 'AC=' + str(ac) + ';AN=' + str(''.join(info_2))
                    # print(info)
                    newline = '\t'.join(line[0:4] + [alt] + line[5:7] + [info] + [line[8]] + [gt_pl])
                    # print(oldline)
                    # print(newline)
                    outfile.writelines(newline)

                else:
                    newline = '\t'.join(line[0:9] + [gt_pl])
                    outfile.writelines(newline)


def circos_main(args):
    if not os.path.isdir(args.save_path):
        os.mkdir(args.save_path)
    print("\n-------------- received arguments -------------")
    print("file name:     %s" % args.name)
    print("save path:     %s" % args.save_path)

    if os.path.exists(str(args.depth_file)):
        print("depth file:    %s" % args.depth_file)
        point2bin(args)

    if os.path.exists(str(args.vcf_file)):
        print("vcf file:      %s" % args.vcf_file)
        print("used genome:   %s\n" % args.genome)
        if args.keep_1alt is True:
            rem_mlt_alt(args)
        else:
            correct_vcf(args)
            vcf_circos(args)
            count_to_freq(args)
            snv_load_per_group(args)
    print("\n%s    ------------ Jobs done! ------------\n")
    sys.exit(0)


if __name__ == '__main__':
    args = get_argparse()
    if args.depth_pos is not None:
        correct_depth(args)
    else:
        circos_main(args)
