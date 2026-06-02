#!/usr/bin/env python3
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
    parser.add_argument('-f', '--sv_percent', type=float, default="0.5", help='percentage of sv reads [0.5]')

    return parser.parse_args()


def validate_file(x):
    if not os.path.exists(x):
        raise argparse.ArgumentTypeError("{0} does not exist".format(x))
    return x


def _get_sv_support(info_field):
    """
    兼容 Sniffles1 和 Sniffles2 的支持读数字段解析。
    优先使用 Sniffles2 的 SUPPORT=，如果没有，再尝试旧的 RE=。
    如果都没有，返回 None（后面用一个保守策略处理）。
    """
    # Sniffles2: INFO 里有 SUPPORT=XX
    m = re.search(r'SUPPORT=(\d+)', info_field)
    if m:
        return int(m.group(1))

    # 一些旧 pipeline: INFO 里有 RE=XX
    m = re.search(r'RE=(\d+)', info_field)
    if m:
        return int(m.group(1))

    # 如果都没有，就返回 None，由调用方决定怎么处理
    return None


def filter_sv(vcf_file, save_path, sv_percent):
    vcf = vcf_file
    file_name = re.split(r'/', vcf)[-1]
    file_name = re.split(r'.vcf', file_name)[0]
    outvcf1 = save_path + '/' + file_name + '.filtered.' + str(sv_percent) + '.vcf'
    outvcf2 = save_path + '/' + file_name + '.filtered.' + str(sv_percent) + '.sorted.vcf'

    with open(vcf, 'r') as infile, open(outvcf1, 'w') as outfile1, open(outvcf2, 'w') as outfile2:
        new_line = []
        for raw in infile:
            # header 行原样抄过去
            if raw.startswith('#'):
                outfile1.writelines(raw)
                outfile2.writelines(raw)
                continue

            # 跳过空行或只有空白的行（包括那个奇怪的 perfect_umi/*）
            if not raw.strip():
                continue

            # 按制表符拆分
            line = re.split(r'\t', raw.strip())
            # 正常 VCF 记录至少 8 列 (#CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO...)
            if len(line) < 8:
                # 例如合并时意外写出的 "perfect_umi/*" 只有一列，这里直接跳过
                continue

            try:
                pos = int(line[1])
            except ValueError:
                # POS 不是整数，直接跳过这条
                continue

            umi_name = line[2]
            # UMI 名字类似 "32_5end_AACAGTGCTGCT"，前面的 32 是 read_count
            try:
                read_count = int(re.split(r'_', umi_name)[0])
            except ValueError:
                # 解析不出 read_count，就跳过（很极端的异常情况）
                continue

            info = line[7]
            sv_read_count = _get_sv_support(info)

            # 如果 read_count 异常（0 或负数），这条记录没法用比例过滤
            if read_count <= 0:
                continue

            # two possible branches:
            # 1) 正常解析出 sv_read_count，用比例过滤；
            # 2) 解析失败（sv_read_count is None），为了保守起见，默认留下这条 SV。
            keep = False
            if sv_read_count is None:
                # 保守策略：不知道支持读数，就不强行过滤，直接保留
                keep = True
            else:
                if sv_read_count / read_count >= sv_percent:
                    keep = True

            if keep:
                # 写入未排序版本
                outfile1.writelines('\t'.join(line) + '\n')
                # 记录下来用于排序输出；record[1] 是整数 pos
                record = tuple([line[0], pos] + line[2:])
                new_line.append(record)

        # 根据 POS 升序排序，再写入 sorted 版本
        new_line.sort(key=itemgetter(1), reverse=False)
        for record in new_line:
            # record = (CHROM, pos(int), ID, REF, ALT, QUAL, FILTER, INFO, FORMAT, SAMPLE...)
            out_fields = [record[0], str(record[1])] + list(record[2:])
            outfile2.writelines('\t'.join(out_fields) + '\n')


def main():
    args = get_argparse()
    if not os.path.isdir(args.save_path):
        os.mkdir(args.save_path)
    print("\n-------------- received arguments -------------")
    print("vcf file      :   %s" % args.vcf_file)
    print("save path     :   %s" % args.save_path)
    print("SV percentage :   %s" % args.sv_percent)

    filter_sv(args.vcf_file, args.save_path, args.sv_percent)
    print("\n%s    ------------ Jobs done! ------------\n")


if __name__ == '__main__':
    main()
