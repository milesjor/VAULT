#!/home/bic/anaconda3/envs/bioconda/bin/python
#
# Chongwei 20191203
# Email: chongwei.bi@kaust.edu.sa

import argparse
import os
import subprocess
import re


def get_argparse():
    parser = argparse.ArgumentParser(description='This is for extracting mapped reads')
    parser.add_argument('-s', '--save_path', type=str, help='path/to/save/ default:[dir/to/fastq]')
    parser.add_argument('-r', '--refer', type=validate_file, required=True, help='path/to/ref.fa')
    parser.add_argument('-q', '--fastq', type=validate_file, required=True, help='path/to/reads.fastq')
    parser.add_argument('-a', '--align_mode', type=str, default='map-ont',
                        choices=['sr', 'map-ont', 'map-pb'], help='parameter in alignment, minimap2 -ax [sr|map-ont|map-pb]')
    parser.add_argument('-t', '--thread', type=int, default='56', help='thread/process number, default: [56]')
    parser.add_argument('-p', '--pe_fastq', type=str, help='read2.fastq for illumina pair-end sequencing')

    args = parser.parse_args()
    if args.pe_fastq is not None:
        args.align_mode = 'sr'

    return args


def validate_file(x):
    if not os.path.exists(x):
        raise argparse.ArgumentTypeError("{0} does not exist".format(x))
    return x


def pre_alignment(args, save_path, file_name):
    cmd = """minimap2 -t {thread} -ax {ax} {refer} {fastq} > {sam} 2>>{log}
             samtools sort -@ {thread} {sam} -o {bam} 2>>{log}
             samtools index -@ {thread} {bam}
             echo -e "\n\n\n{name}_samtools_flagstat\t{sam}" >> {log}
             samtools flagstat {sam} >> {log}
             grep -v ^@ {sam} | cut -f 3 | sort | uniq -c | sort -n >> {log}""".format(thread=args.thread,
                                                                                       ax=args.align_mode,
                                                                                       refer=args.refer,
                                                                                       fastq=args.fastq,
                                                                                       name=file_name,
                                                                                       sam=save_path + file_name + '.sam',
                                                                                       bam=save_path + file_name + '.bam',
                                                                                       log=save_path + file_name + '_alignment_summary.log')
    subprocess.run(cmd, shell=True, check=True, stdout=subprocess.PIPE)


def extract_reads_by_name(args, save_path, file_name):
    cmd = """cat {sam} | grep -v ^@ | cut -f1,3 | grep -v "*" | cut -f1 | sort | uniq > {name_lst}
             seqtk subseq {fastq} {name_lst} | seqtk seq -l0 > {mapped_fastq}""".format(fastq=args.fastq,
                                                                                        sam=save_path + file_name + '.sam',
                                                                                        name_lst=save_path + file_name + '.mapped.lst',
                                                                                        mapped_fastq=save_path + file_name + '.mapped.fastq')
    subprocess.run(cmd, shell=True, check=True, stdout=subprocess.PIPE)

    if args.pe_fastq is not None:
        cmd = "seqtk subseq {fastq} {name_lst} | seqtk seq -l0 > {mapped_fastq}".format(fastq=args.pe_fastq,
                                                                                        name_lst=save_path + file_name + '.mapped.lst',
                                                                                        mapped_fastq=save_path + file_name + '.read2.mapped.fastq')
        subprocess.run(cmd, shell=True, check=True, stdout=subprocess.PIPE)


def get_name(args):
    save_path = args.save_path

    split = re.split(r'/', args.fastq)
    file_name = split[-1]
    file_name = re.split(r'.fastq', file_name)[0]

    if not save_path:
        if not split[1:]:
            save_path = './'
        else:
            save_path = "".join(split[0:-1])

    if not os.path.isdir(save_path):
        os.mkdir(save_path)
    save_path = save_path + "/"

    return save_path, file_name


def main():
    args = get_argparse()
    [save_path, file_name] = get_name(args)

    print("\n-------------- received arguments -------------")
    print("fastq file :   %s" % args.fastq)
    if args.pe_fastq is not None:
        print("read2 file :   %s" % args.pe_fastq)
    print("save path  :   %s" % args.save_path)
    print("align mode :   %s" % args.align_mode)
    print("refer file :   %s\n" % args.refer)

    pre_alignment(args, save_path, file_name)
    extract_reads_by_name(args, save_path, file_name)


if __name__ == '__main__':
    main()