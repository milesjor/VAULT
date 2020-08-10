#!/usr/bin/env python3
# Chongwei 20190619
# Email: bicwei@gmail.com

import argparse
import subprocess
import threading
import os
import multiprocessing
import re
import glob
import sys
from datetime import datetime
from _version import vault_version
from tools import extract_mapped_reads, filter_sv, umi_group_filter, call_consensus, change_VCF_pos, draw_circos
from variants_calling import snp_calling, sv_calling


def get_argparse():
    args = argparse.ArgumentParser(description='This is for analyzing UMI labeled reads in IDMseq and IMTseq.',
                                   epilog='Refer github https://github.com/milesjor/VAULT for more detail')
    subparsers = args.add_subparsers(title='subcommands',
                                     description='valid subcommands',
                                     help='additional help',
                                     dest='subcommands_name')
    args.add_argument('-v', '--version', action='version', version=vault_version, help='show the current version')

    req = args.add_argument_group(title='Required options')
    req.add_argument('-u', '--umi_adapter', type=str, help='UMI sequence, automatically detect NNNATGCNNN as UMI')
    req.add_argument('-s', '--save_path', type=str, help='path/to/save/')
    req.add_argument('-r', '--refer', type=validate_file, help='path/to/ref.fa')
    req.add_argument('-q', '--fastq', type=validate_file, help='path/to/reads.fastq, or fastq.gz')

    opt = args.add_argument_group(title='Optional options')
    opt.add_argument('-e', '--error', type=str, default='0.11', help='error tolerate rate in umi analysis '
                                                                     '(depend on read error rate) [0.11]')
    opt.add_argument('-t', '--thread', type=int, default='5', help='parallel process number [5]')
    opt.add_argument('-T', '--threshold', type=int, default='5', help='threshold of read number for snp analysis [5]')
    opt.add_argument('-b', '--bash_thread', type=int, default='1', help='thread for running bash cmd [1]')
    opt.add_argument('-F', '--allele_freq', type=float, default='0.67', help='filter SNVs '
                                                                             'by allelic frequency [0.67]')
    opt.add_argument('-f', '--sv_freq', type=float, default='0.5', help='filter SVs by allelic frequency [0.5]')
    opt.add_argument('-p', '--pe_fastq', type=str, help='read2.fastq for illumina pair-end sequencing')
    opt.add_argument('-a', '--align_mode', type=str, default='map-ont',
                     choices=['sr', 'map-ont', 'map-pb'], help='parameter in alignment, '
                                                               'minimap2 -ax [sr|map-ont|map-pb]')
    opt.add_argument('--unmapped_reads', action='store_true', help='extract mapped reads before UMI analysis')
    opt.add_argument('--group_filter', action='store_true', help='filter out low-confidence UMI groups')

    # subcommands
    consensus = subparsers.add_parser('consensus',
                                      help='get consensus sequence from VAULT result',
                                      description='This is for generating consensus sequence from VAULT result '
                                                  'using canu + medaka. '
                                                  'Just input the path of VAULT result folder '
                                                  '--> the folder containing (/snp /grouped_reads /umi_analysis). ')
    consensus.set_defaults(func=call_consensus.consensus_main)
    consensus.add_argument('-s', '--save_path', type=str, required=True, help='path/to/VAULT_result_folder/')
    consensus.add_argument('-t', '--thread', type=int, default='6', help='parallel worker [6]')
    consensus.add_argument('-T', '--sub_thread', type=int, default='4', help='thread for every parallel worker [4]')
    consensus.add_argument('--threshold', type=int, default='20', help='process UMI group with reads more than [20]')

    position = subparsers.add_parser('position',
                                     help='correct the position and chr_name in vcf file',
                                     description='It will correct the position and chr_name in vcf file '
                                                 'to enable variant annotation. It can also reverse the coordinate'
                                                 ' and DNA base if the used referene in VAULT '
                                                 'is reverse complimentary of the referene genome.')
    position.set_defaults(func=change_VCF_pos.position_main)
    position.add_argument('-v', '--vcf_file', type=validate_file, required=True, help='path/to/file.vcf')
    position.add_argument('-c', '--chr_name', type=str, required=True, help='chromosome name')
    position.add_argument('-p', '--pos_change', type=str, required=True, help='position change, e.g. +12300 or -55789')
    position.add_argument('-b', '--total_base', type=int,
                          help='the length of reference genome used. When provided, it will reverse the '
                               'coordinate(position) and also do reverse complimentary of the DNA bases in vcf file')
    position.add_argument('-s', '--save_path', type=str, default="./", help='path/to/save/')

    circos = subparsers.add_parser('circos',
                                   help='prepare data for circos, used in IMTseq',
                                   description='This is for preparing data for circos. It is used in IMTseq. '
                                               'For vcf file, it will correct position and remove indel. '
                                               'For depth file, it will bin depth based on '
                                               'user defined bin size.')
    circos.set_defaults(func=draw_circos.circos_main)
    circos.add_argument('-n', '--name', type=str, default="circos", help='prefix of output file')
    circos.add_argument('-d', '--depth_file', type=validate_file, help='depth file from samtools depth')
    circos.add_argument('-s', '--save_path', type=str, required=True, help='path/to/save/')
    circos.add_argument('-v', '--vcf_file', type=validate_file, help='path/to/file.vcf')
    circos.add_argument('-b', '--bin_size', type=int, default="30", help='how many bases per bin')
    circos.add_argument('-g', '--genome', type=str,
                        choices=['mmt_nod_F6_10N_to_C57', 'hchrM.F9.UMIs', 'hchrM.F6.UMIs', 'mmt_c57_F6_10N'],
                        help='the genome used in VAULT')
    circos.add_argument('-c', '--chr_name', type=str, help='chromosome name showed in vcf file')
    circos.add_argument('-A', '--keep_1alt', action='store_true',
                        help='leave only 1 Alt in vcf file for helping SNV annotation')

    group_filter = subparsers.add_parser('filter',
                                         help='filter out low-confidence UMI groups',
                                         description='This is for filtering out low-confidence UMI groups after '
                                                     'VAULT analysis. It is the same as <vault --group_filter>')
    group_filter.set_defaults(func=umi_group_filter.filter_main)
    group_filter.add_argument('-s', '--save_path', type=str, required=True,
                              help='path/to/snp (the path to snp folder generated by VAULT)')
    group_filter.add_argument('-F', '--allele_freq', type=float, default='0.67', help='filter SNVs '
                                                                                      'by allelic frequency [0.67]')

    args = args.parse_args()

    if args.subcommands_name in ["consensus", "position", "circos", "filter"]:
        args.func(args)

    else:
        if args.pe_fastq is not None:
            args.align_mode = 'sr'
        if args.umi_adapter is None or args.save_path is None or args.refer is None or args.fastq is None:
            sys.stderr.write("usage: vault [-h] [-v] -u UMI_ADAPTER -s SAVE_PATH -r REFER -q FASTQ\n"
                             "             [-e ERROR] [-t THREAD] [-T THRESHOLD] [-b BASH_THREAD]\n"
                             "             [-F ALLELE_FREQ] [-f SV_FREQ] [-a {sr,map-ont,map-pb}]\n"
                             "             [-p PE_FASTQ] [--unmapped_reads] [--group_filter]\n"
                             "             {consensus} ...\n"
                             "vault: error: the following arguments are required: "
                             "-u/--umi_adapter, -s/--save_path, -r/--refer, -q/--fastq\n")
            sys.exit(1)

    return args


def validate_file(x):
    if not os.path.exists(x):
        raise argparse.ArgumentTypeError("{0} does not exist".format(x))
    return x


def analyze_umi(args):
    # divide UMI adapter to left flank, umi, right flank
    def find_all_index(adapt, search='N'):
        index = [x for x in range(adapt.find(search), len(adapt)) if adapt[x] == search]
        return index

    args.umi_adapter = args.umi_adapter.upper()
    n_index = find_all_index(args.umi_adapter)
    left_flank = args.umi_adapter[0:n_index[0]]
    right_flank = args.umi_adapter[n_index[-1]+1:]
    umi = args.umi_adapter[n_index[0]:n_index[-1]+1]
    print("left_flank:  ", left_flank)
    print("umi:         ", umi)
    print("right_flank:  %s" % right_flank)

    if not os.path.isdir(args.save_path):
        os.mkdir(args.save_path)

    def run_bash_cmd(name, fastq, thread, end):
        cmd = """sh {path}/check_umi.sh -n {name} -l {left_flank} -u {umi} -r {right_flank} \
                 -q {fastq} -s {save} -e {error} -t {thread} -{end}""".format(path=sys.path[0],
                                                                              name=name,
                                                                              left_flank=left_flank,
                                                                              umi=umi,
                                                                              right_flank=right_flank,
                                                                              fastq=fastq,
                                                                              save=args.save_path,
                                                                              error=args.error,
                                                                              thread=thread,
                                                                              end=end)
        subprocess.run(cmd, shell=True, check=True, stdout=subprocess.PIPE)  # print only stderr to screen

    py_threads = []
    if args.pe_fastq is None:
        thread = int(args.thread) // 2
        t1 = threading.Thread(target=run_bash_cmd, args=("umi_analysis", args.fastq, thread, str(5)))
        t2 = threading.Thread(target=run_bash_cmd, args=("umi_analysis", args.fastq, thread, str(3)))
        py_threads.append(t1)
        py_threads.append(t2)
    else:
        thread = int(args.thread) // 4
        t1 = threading.Thread(target=run_bash_cmd, args=("umi_analysis", args.fastq, thread, str(5)))
        t2 = threading.Thread(target=run_bash_cmd, args=("umi_analysis", args.fastq, thread, str(3)))
        t3 = threading.Thread(target=run_bash_cmd, args=("umi_analysis2", args.pe_fastq, thread, str(5)))
        t4 = threading.Thread(target=run_bash_cmd, args=("umi_analysis2", args.pe_fastq, thread, str(3)))
        py_threads.append(t1)
        py_threads.append(t2)
        py_threads.append(t3)
        py_threads.append(t4)

    for t in py_threads:
        t.start()

    for t in py_threads:
        t.join()
    return umi


def group_read_name(seq, sav, count, name_lst):
    outfile = sav + '/reads_with_same_UMIs/' + count + '_' + seq + '.lst'
    with open(outfile, 'w') as of:
        for name in name_lst:
            of.writelines([name, "\n"])


def extract_name_lst(args):
    end5 = args.save_path + '/umi_analysis/5end_UMIs'
    end3 = args.save_path + '/umi_analysis/3end_UMIs'
    if not os.path.isdir(end5 + '/reads_with_same_UMIs'):
        os.mkdir(end5 + '/reads_with_same_UMIs')
    if not os.path.isdir(end3 + '/reads_with_same_UMIs'):
        os.mkdir(end3 + '/reads_with_same_UMIs')

    if args.pe_fastq is None:
        file5 = end5 + '/umi_analysis.07.1read_count.2umi.3read_name.lst'
    else:
        file5 = args.save_path + '/umi_analysis/5end_UMIs/umi_analysis.07.1read_count.2umi.3read_name.pe.lst'

        # # Below is for the specific Illumina data set
        # file5 = end5 + '/umi_analysis.07.1read_count.2umi.3read_name.lst'

    file3 = end3 + '/umi_analysis.07.1read_count.2umi.3read_name.lst'

    # multiprocessing module for extracting read name list with same umi.
    p = multiprocessing.Pool(args.thread)
    with open(file5, "r") as infile5:
        for record in infile5:
            record_splt = record.strip().split(' ')
            count5 = record_splt[0]
            umi = record_splt[1]
            read_name_lst = record_splt[2:]
            p.apply_async(group_read_name, args=(umi, end5, count5, read_name_lst))

    if args.pe_fastq is None:
        with open(file3, "r") as infile3:
            for record in infile3:
                record_splt = record.strip().split(' ')
                count3 = record_splt[0]
                umi = record_splt[1]
                read_name_lst = record_splt[2:]
                p.apply_async(group_read_name, args=(umi, end3, count3, read_name_lst))
    p.close()
    p.join()


def lst2fastq_variant(args, umi_n5):
    umi_regex5 = '^' + umi_n5 + '$'
    umi_regex5 = re.sub('N', '[ATGC]', umi_regex5)

    complement_dict = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    umi_n3 = "".join(complement_dict.get(base, base) for base in reversed(umi_n5))
    umi_regex3 = '^' + umi_n3 + '$'
    umi_regex3 = re.sub('N', '[ATGC]', umi_regex3)

    if args.pe_fastq is None:
        file5 = args.save_path + '/umi_analysis/5end_UMIs/umi_analysis.07.1read_count.2umi.3read_name.lst'
    else:
        file5 = args.save_path + '/umi_analysis/5end_UMIs/umi_analysis.07.1read_count.2umi.3read_name.pe.lst'

        # # Below is for the specific Illumina data set
        # file5 = args.save_path + '/umi_analysis/5end_UMIs/umi_analysis.07.1read_count.2umi.3read_name.lst'

    if not os.path.isdir(args.save_path + '/grouped_reads'):
        os.mkdir(args.save_path + '/grouped_reads')
        os.mkdir(args.save_path + '/grouped_reads/perfect_umi')
        os.mkdir(args.save_path + '/grouped_reads/perfect_umi/5end_and_3end')
        os.mkdir(args.save_path + '/grouped_reads/perfect_umi/only_5end')
        os.mkdir(args.save_path + '/grouped_reads/perfect_umi/only_3end')

    if not os.path.isdir(args.save_path + '/snp'):
        os.mkdir(args.save_path + '/snp')
        os.mkdir(args.save_path + '/snp/perfect_umi')

    # multiprocessing module for combining read name from same 5' and 3' umi, and extract reads, and analyze snp
    p = multiprocessing.Pool(args.thread)
    with open(file5, "r") as infile5:
        for seq in infile5:
            seq_splt = seq.strip().split(' ')
            count5 = seq_splt[0]
            umi = seq_splt[1]
            p.apply_async(use_class_5end, args=(umi, args.save_path, args.fastq, count5, umi_regex5, args.threshold,
                                                args.refer, args.bash_thread, args.align_mode, args.pe_fastq))
    p.close()
    p.join()
    print("%s    5'+3' and 5' end UMI analysis finished" % datetime.now().strftime("%Y-%m-%d %H:%M:%S"))

    # multiprocessing module for processing unused 3' end umi, and extract reads, and analyze snp
    if args.pe_fastq is None:
        p = multiprocessing.Pool(args.thread)
        for file in os.listdir(args.save_path + '/umi_analysis/3end_UMIs/reads_with_same_UMIs/'):
            if file.endswith('.lst'):
                name_lst3 = args.save_path + '/umi_analysis/3end_UMIs/reads_with_same_UMIs/' + file
                umi3 = re.split('[_.]', file)[1]
                count3 = re.split('[_.]', file)[0]
                p.apply_async(use_class_3end, args=(umi3, args.save_path, args.fastq, count3, umi_regex3, args.threshold,
                                                    args.refer, args.bash_thread, args.align_mode, args.pe_fastq, name_lst3))
        p.close()
        p.join()
        print("%s    The rest 3' end UMI analysis finished" % datetime.now().strftime("%Y-%m-%d %H:%M:%S"))


def use_class_5end(a, b, c, d, e, f, g, h, i, j=''):
    process_5end = ProcessUmi(a, b, c, d, e, f, g, h, i, j)
    process_5end.process_5()


def use_class_3end(a, b, c, d, e, f, g, h, i, j='', k=''):
    process_3end = ProcessUmi(a, b, c, d, e, f, g, h, i, j, k)
    process_3end.process_3()


class ProcessUmi:
    umi = ''
    save = ''
    fastq = ''
    count = ''
    umi_regex = ''
    threshold = ''
    refer = ''
    bash = ''

    def __init__(self, umi, save, fastq, count, umi_regex, threshold, refer, bash, ax, read2='', name_lst3=''):
        self.umi = umi
        self.save = save
        self.fastq = fastq
        self.count = count
        self.umi_regex = umi_regex
        self.threshold = threshold
        self.refer = refer
        self.bash = bash
        self.ax = ax
        self.read2 = read2
        self.name_lst3 = name_lst3

    def get_3umi_file_name(self):
        umi5 = self.umi
        name_lst5 = self.save + '/umi_analysis/5end_UMIs/reads_with_same_UMIs/' + self.count + '_' + umi5 + '.lst'
        complement_dict = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
        umi3 = "".join(complement_dict.get(base, base) for base in reversed(umi5))
        name_lst3 = self.save + '/umi_analysis/3end_UMIs/reads_with_same_UMIs/*' + '_' + umi3 + '.lst'
        return [name_lst5, name_lst3, umi3]

    def combine_2_lst(self, name_lst5, name_lst3, umi5, umi3, folder):
        name_lst3 = glob.glob(str(name_lst3))[0]
        name_lst = self.save + '/grouped_reads/' + folder + '/5end_and_3end/' + umi5 + '_' + umi3 + '.lst'
        cmd = 'cat ' + name_lst5 + ' ' + name_lst3 + ' | sort | uniq > ' + name_lst
        subprocess.run(cmd, shell=True, check=True, stdout=subprocess.PIPE)
        read_count = sum(1 for line in open(name_lst, 'r'))
        name_count = self.save + '/grouped_reads/' + folder + '/5end_and_3end/' + str(
            read_count) + '_' + umi5 + '_' + umi3 + '.lst'
        os.rename(name_lst, name_count)
        os.remove(name_lst5)
        os.remove(name_lst3)
        return [read_count, name_count]

    def lst2fastq_snp(self, read_count, umi_name, name_lst, read_folder, snp_folder, end='_'):
        of = self.save + '/grouped_reads/' + read_folder + str(read_count) + '_' + umi_name + '.fastq'
        cmd = """seqtk subseq {fastq} {name_lst} > {outfile}""".format(fastq=self.fastq,
                                                                       name_lst=name_lst,
                                                                       outfile=of)
        subprocess.run(cmd, shell=True, check=True, stdout=subprocess.PIPE)
        fastq_file = of

        # Below is for normal data set
        if self.read2 is not None:
            of2 = self.save + '/grouped_reads/' + read_folder + str(read_count) + '_' + umi_name + '.read2.fastq'
            cmd = """seqtk subseq {fastq} {name_lst} > {outfile}""".format(fastq=self.read2,
                                                                           name_lst=name_lst,
                                                                           outfile=of2)
            subprocess.run(cmd, shell=True, check=True, stdout=subprocess.PIPE)
            fastq_file = of + ' ' + of2

        # # Below is for the specific Illumina data set
        # if self.read2 is not None:
        #     of3 = self.save + '/grouped_reads/' + read_folder + str(read_count) + '_' + umi_name + '.read2.raw.fastq'
        #     of2 = self.save + '/grouped_reads/' + read_folder + str(read_count) + '_' + umi_name + '.read2.fastq'
        #     cmd = """seqtk subseq {fastq} {name_lst} > {outfile}
        #              sh /home/bic/bic/downloads/bbmap/bbduk.sh in={outfile} out={outfile2} ftr=30
        #              """.format(fastq=self.read2,
        #                         name_lst=name_lst,
        #                         outfile=of3,
        #                         outfile2=of2)
        #     subprocess.run(cmd, shell=True, check=True, stdout=subprocess.PIPE)
        #     fastq_file = of + ' ' + of2

        os.remove(name_lst)

        if not os.path.isdir(self.save + snp_folder):
            os.mkdir(self.save + snp_folder)
        save_path = self.save + snp_folder + '/' + str(read_count) + end + umi_name
        os.mkdir(save_path)
        os.mkdir(save_path + '/vcf')
        name = str(read_count) + end + umi_name

        snp_calling(save_path, self.bash, self.ax, self.refer, fastq_file, name)
        if self.ax != 'sr':
            sv_calling(save_path, self.bash, self.refer, fastq_file, name)

    def process_5(self):
        result = self.get_3umi_file_name()
        name_lst5 = result[0]
        name_lst3 = result[1]
        umi5 = self.umi
        umi3 = result[2]

        if re.match(str(self.umi_regex), str(umi5)):
            folder = 'perfect_umi'
            if len(glob.glob(str(name_lst3))) == 1:
                result = self.combine_2_lst(name_lst5, name_lst3, umi5, umi3, folder)
                read_count = result[0]
                name_lst = result[1]
                if int(read_count) >= int(self.threshold):
                    umi_name = umi5 + '_' + umi3
                    read_folder = folder + '/5end_and_3end/'
                    snp_folder = '/snp/' + folder
                    self.lst2fastq_snp(read_count, umi_name, name_lst, read_folder, snp_folder)

            elif len(glob.glob(str(name_lst3))) > 1:
                print("\n!!!!!! multiple match for read_name.lst file: %s!!!!!!\n" % name_lst3)
                if int(self.count) >= int(self.threshold):
                    read_folder = folder + '/only_5end/'
                    snp_folder = '/snp/' + folder
                    self.lst2fastq_snp(self.count, self.umi, name_lst5, read_folder, snp_folder, str('_5end_'))

            else:
                if int(self.count) >= int(self.threshold):
                    read_folder = folder + '/only_5end/'
                    snp_folder = '/snp/' + folder
                    self.lst2fastq_snp(self.count, self.umi, name_lst5, read_folder, snp_folder, str('_5end_'))

    def process_3(self):
        read_count = sum(1 for line in open(self.name_lst3, 'r'))
        if re.match(str(self.umi_regex), str(self.umi)):
            if read_count >= self.threshold:
                read_folder = '/perfect_umi/only_3end/'
                snp_folder = '/snp/perfect_umi'
                self.lst2fastq_snp(read_count, self.umi, self.name_lst3, read_folder, snp_folder, str('_3end_'))


def final_clean_up(args):
    # process SNV and InDel vcf
    snp_regex = args.save_path + '/snp/perfect_umi/*/vcf/*snp.vcf'
    all_snp = args.save_path + '/snp/all_snp_from_perfect_umi.vcf'
    pass_snp = args.save_path + '/snp/pass_snp_from_perfect_umi.vcf'

    all_snp_pcent = args.save_path + '/snp/all_snp_from_perfect_umi.pcent.vcf'
    all_snp_pcent_flt = args.save_path + '/snp/all_snp_from_perfect_umi.pcent.flt.vcf'

    find_vcf_files(snp_regex, all_snp)

    save_dir = args.save_path + '/snp'
    umi_group_filter.add_pcent(save_dir, all_snp_pcent)

    if args.group_filter is True:
        awk_filter = "awk '{if($3>0.35 && $3<0.65)print}' | awk '{if($7 < 0.4) print}'"
        cmd = """cat {save_dir}/umi_group.flt.summary.txt | {filter} > {save_dir}/wrong.group.summary.txt
                 cat {save_dir}/wrong.group.summary.txt | cut -f1 | sort -k1 > {save_dir}/wrong.group.lst
                 ls {save_dir}/perfect_umi | sort | uniq > {save_dir}/all.group.lst
                 comm -23 {save_dir}/all.group.lst {save_dir}/wrong.group.lst > {save_dir}/pass.group.lst
                 
                 split -l 50 {save_dir}/pass.group.lst {save_dir}/pattern-file.split.
                 cat {save_dir}/all_snp_from_perfect_umi.pcent.vcf | grep "^#" > {save_dir}/all_snp_from_perfect_umi.pcent.rem.vcf
                 for CHUNK in {save_dir}/pattern-file.split.* ; do
                    cat {save_dir}/all_snp_from_perfect_umi.pcent.vcf | grep -v "^#" | grep -f "$CHUNK" \
                        >> {save_dir}/all_snp_from_perfect_umi.pcent.rem.vcf
                 done
                 rm {save_dir}/pattern-file.split.* {save_dir}/all.group.lst
                 
                 cat {save_dir}/all_snp_from_perfect_umi.pcent.rem.vcf | \
                    bcftools filter -s LowQual -e '%QUAL<10 || INFO/DP<3 || INFO/IMF<0.5 || INFO/VARP<{alle_freq} || INFO/ALTC<3' \
                    > {save_dir}/all_snp_from_perfect_umi.pcent.rem.flt.vcf
                    
                 cat {save_dir}/all_snp_from_perfect_umi.pcent.rem.flt.vcf | grep -v LowQual \
                    > {save_dir}/pass_snp_from_perfect_umi.flt.vcf
                 find {save_dir}/perfect_umi/ -name "*coverage.3plus.txt" -print0 | xargs -0 cat > {save_dir}/coverage.3plus.txt
                 """.format(save_dir=save_dir,
                            alle_freq=args.allele_freq,
                            filter=awk_filter)

    else:
        cmd = """cat {all_snp_pcent} | bcftools filter -s LowQual -e '%QUAL<10 || INFO/DP<3 || INFO/IMF<0.5 || \
                    INFO/VARP<{alle_freq} || INFO/ALTC<3' > {fltvcf}
                 rm {all_snp_pcent} {save}/snp/umi_group.flt.summary.txt
                 cat {fltvcf} | grep -v LowQual > {pass_snp}
                 find {save}/snp/perfect_umi/ -name "*coverage.3plus.txt" -print0 | xargs -0 cat > {save}/snp/coverage.3plus.txt
                 """.format(all_snp_pcent=all_snp_pcent,
                            fltvcf=all_snp_pcent_flt,
                            save=args.save_path,
                            alle_freq=args.allele_freq,
                            pass_snp=pass_snp)

    subprocess.run(cmd, shell=True, check=True, stdout=subprocess.PIPE)

    # process SV vcf
    if args.align_mode != 'sr':
        sv_regex = args.save_path + '/snp/perfect_umi/*/vcf/*sv.vcf'
        all_sv = args.save_path + '/snp/all_sv_from_perfect_umi.vcf'
        find_vcf_files(sv_regex, all_sv)

        save_path = args.save_path + '/snp/'
        filter_sv.filter_sv(all_sv, save_path, args.sv_freq)


def find_vcf_files(file_regex, allvcf):
    file_list = glob.glob(str(file_regex))
    if len(file_list) >= 1:
        with open(allvcf, 'w') as of, open(file_list[0], 'r') as header:
            for line in header:
                if line.startswith('##'):
                    of.writelines(line)
                elif line.startswith('#'):
                    line = re.split(r'perfect_umi', line)
                    new_line = line[0] + 'perfect_umi/*\n'
                    of.writelines(''.join(new_line))

            for vcf in file_list:
                vcf_file = open(vcf, 'r')
                sample_id = re.split(r'(perfect_umi/|/vcf)', vcf)[2]
                for line in vcf_file:
                    if not line.startswith('#'):
                        line = re.split(r'\t', line)
                        of.writelines('\t'.join((line[0], line[1], sample_id, '\t'.join(line[3:]))))

    else:
        sys.stderr.write("No VCF file detected! Please check %s\n" % file_regex)
        sys.exit(1)


def combine_pe_umi(args):
    if args.pe_fastq is not None:
        umi_csv_5end1 = args.save_path + '/umi_analysis/5end_UMIs/umi_analysis.06.1umi.2read_name.csv'
        umi_csv_5end2 = args.save_path + '/umi_analysis2/5end_UMIs/umi_analysis2.06.1umi.2read_name.csv'
        name_lst5 = args.save_path + '/umi_analysis/5end_UMIs/umi_analysis.07.1read_count.2umi.3read_name.pe.lst'

        prepare_umi_file = """awk -F "," '
                            {
                                k=$2
                                for (i=3;i<=NF;i++)
                                    k=k" "$i
                                if (! a[$1])
                                    a[$1]=k
                                else
                                    a[$1]=a[$1]" "k
                            }
                            END{
                                for (i in a)
                                    print i" "a[i]
                            }' | awk '{print NF-1" "$0}' """

        cmd = """cp {umi_csv_5end1} {umi_csv_5end1}.read1 
                 cp {umi_csv_5end2} {umi_csv_5end1}.read2
                 cat {umi_csv_5end2} >> {umi_csv_5end1} 
                 cat {umi_csv_5end1} | sort | uniq | {prepare_umi_file} > {name_lst5}  
                 """.format(umi_csv_5end2=umi_csv_5end2,
                            umi_csv_5end1=umi_csv_5end1,
                            prepare_umi_file=prepare_umi_file,
                            name_lst5=name_lst5)
        subprocess.run(cmd, shell=True, check=True, stdout=subprocess.PIPE)


def check_tools():
    if call_consensus.check_tool("minimap2") is not True:
        sys.exit("ERROR! Executable minimap2 is not found!\n"
                 "ERROR! Please install minimap2 -> `conda install minimap2=2.11`")
    if call_consensus.check_tool("cutadapt") is not True:
        sys.exit("ERROR! Executable cutadapt is not found!\n"
                 "ERROR! Please install cutadapt -> `conda install cutadapt=2.7`")
    if call_consensus.check_tool("samtools") is not True:
        sys.exit("ERROR! Executable samtools is not found!\n"
                 "ERROR! Please install samtools -> `conda install samtools=1.9`")
    if call_consensus.check_tool("bcftools") is not True:
        sys.exit("ERROR! Executable bcftools is not found!\n"
                 "ERROR! Please install bcftools -> `conda install bcftools=1.9`")
    if call_consensus.check_tool("seqtk") is not True:
        sys.exit("ERROR! Executable seqtk is not found!\n"
                 "ERROR! Please install seqtk -> `conda install seqtk=1.3`")
    if call_consensus.check_tool("sniffles") is not True:
        sys.exit("ERROR! Executable sniffles is not found!\n"
                 "ERROR! Please install sniffles -> `conda install sniffles=1.0.11`")


def main():
    args = get_argparse()
    check_tools()
    print("\n%s    --------- received arguments ---------" % datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
    print("umi adapter:   ", args.umi_adapter)
    print("fastq file:    ", args.fastq)
    if args.pe_fastq is not None:
        print("read2 file:    ", args.pe_fastq)
    print("refer file:    ", args.refer)
    print("error rate:    ", args.error)
    print("align mode:    ", args.align_mode)
    print("save path:      %s\n" % args.save_path)

    if args.unmapped_reads is True:
        [save_path, file_name] = extract_mapped_reads.get_name(args)
        extract_mapped_reads.pre_alignment(args, save_path, file_name)
        extract_mapped_reads.extract_reads_by_name(args, save_path, file_name)
        args.fastq = save_path + file_name + '.mapped.fastq'
        if args.pe_fastq is not None:
            args.pe_fastq = save_path + file_name + '.read2.mapped.fastq'

    print("%s    Start extracting UMIs from reads ...\n" % datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
    umi_n5 = analyze_umi(args)

    # Below is for the specific Illumina data set
    # combine_pe_umi(args)

    combine_pe_umi(args)
    print("\n%s    Extract UMIs from reads finished" % datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
    extract_name_lst(args)
    print("%s    Group reads by UMIs finished" % datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
    # umi_n5 = 'NNNNNNNNNN'
    print("%s    Start individual UMI group analysis ..." % datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
    lst2fastq_variant(args, umi_n5)
    final_clean_up(args)
    print("\n%s    ------------ Jobs done! ------------\n" % datetime.now().strftime("%Y-%m-%d %H:%M:%S"))


if __name__ == '__main__':
    main()
