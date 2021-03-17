import subprocess
import re
import io
import numpy as np
import pandas as pd
import os
import sys


def rem_mlt_alt(save_path):
    vcf = save_path + '/pass_snp_from_perfect_umi.flt.vcf'
    file_name = re.split(r'/', vcf)[-1]
    file_name = re.split(r'.vcf', file_name)[0]
    outvcf = save_path + '/' + file_name + '.1alt.vcf'
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


def count_snv(save_path, name):
    vcf = save_path + '/pass_snp_from_perfect_umi.flt.1alt.vcf'
    snv_vcf = save_path + '/pass_snp_from_perfect_umi.flt.1alt.snv.vcf'
    out = save_path + '/' + name + '_snv.count.txt'

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
    # print(pd_vcf.POS.unique())
    # print(pd_vcf.POS.count())
    # print(pd_vcf.POS.value_counts())
    # print(pd_vcf.POS.value_counts().reset_index())
    # count_pos = pd_vcf.POS.value_counts().reset_index().rename(columns={'POS': 'count', 'index': 'POS'})

    count_pos = pd_vcf.groupby(["POS", "REF", "ALT"]).size().reset_index(name="count")
    sort_count = count_pos.sort_values('POS').reset_index(drop=True)
    sort_count.insert(0, 'chr', 'mt1')
    snp_df = sort_count[['chr', 'POS', 'POS', 'count', 'REF', 'ALT']]
    np.savetxt(out, snp_df, fmt='%s')


def count_to_freq(save_path, name):
    position = save_path + '/all.coverage.3plus.pos.txt'
    snv_count = save_path + '/' + name + '_snv.count.txt'
    snv_freq = save_path + '/' + name + '_snv.freq.txt'

    with open(position, 'r') as infile1, open(snv_count, 'r') as infile2, open(snv_freq, 'w') as outfile:
        dir = {}
        for line in infile1:
            v, k = line.strip().split(' ')
            dir[k.strip()] = v.strip()
            # print(dir)

        for line2 in infile2:
            line2_split = line2.strip().split(' ')
            count = line2_split[3]
            pos = line2_split[1]
            depth = dir[pos]
            freq = round(int(count) / int(depth), 6)
            new_line2 = " ".join([line2_split[0], line2_split[1], line2_split[2], str(freq),
                                  line2_split[3], line2_split[4], line2_split[5]]) + '\n'
            outfile.writelines(new_line2)


def get_vaf(args):
    print("\n-------------- received arguments -------------")
    print("file name:     %s" % args.name)
    print("save path:     %s" % args.save_path)

    infolder = args.save_path
    name = args.name

    # infolder = "/Users/bic/Desktop/codes/github/VAULT_local/example/result/snp"
    # name = "VAFcalculator"

    # run UMI group filter
    if not os.path.exists(infolder + '/pass.group.lst'):
        sys.stderr.write('ERROR! File not exist -> %s' % str(infolder + '/pass.group.lst'))
        sys.stderr.write('ERROR! Please run `vault filter` first')
        sys.exit(1)

    # get the UMI group number that cover every position of the amplicon
    cmd = r"""save=%s
            mkdir ${save}/coverage
            for line in $(cat ${save}/pass.group.lst); do cp ${save}/perfect_umi/${line}/${line}.coverage.txt ${save}/coverage ; done
            for file in ${save}/coverage/*.coverage.txt; do cat $file | awk '{if($3>=3)print $2}' >> ${save}/all.coverage.3plus.pos.txt.tmp ; done
            cat ${save}/all.coverage.3plus.pos.txt.tmp | sort -n | uniq -c > ${save}/all.coverage.3plus.pos.txt
            rm ${save}/all.coverage.3plus.pos.txt.tmp""" % infolder

    subprocess.run(cmd, shell=True, check=True, stdout=subprocess.PIPE)

    # remove multiple allele in VCF files
    rem_mlt_alt(infolder)

    # count the number of UMI groups with the same SNV
    count_snv(infolder, name)

    # get VAF, VAF = (UMI group number with the same SNV) / (UMI groups number covering the SNV position)
    count_to_freq(infolder, name)

    print("\nCalculate VAF done!\n")
    sys.exit(0)
