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
import arguments
from datetime import datetime
from tools import extract_mapped_reads, filter_sv, umi_group_filter, call_consensus, change_VCF_pos, draw_circos
from tools import read_length_filter
from variants_calling import snp_calling, sv_calling
from check_umi2 import run_umi_analysis
import logging


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
    logging.info("left_flank:   %s" % left_flank)
    logging.info("umi:          %s" % umi)
    logging.info("right_flank:  %s" % right_flank)

    if not os.path.isdir(args.save_path):
        os.mkdir(args.save_path)

    def run_bash_cmd(name, fastq, thread, end):
        cmd = """bash {path}/check_umi.sh -n {name} -l {left_flank} -u {umi} -r {right_flank} \
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

        if args.py_only is True:
            t1 = threading.Thread(target=run_umi_analysis, args=("umi_analysis", args.fastq, thread, str(5),
                                                                 "g", args.save_path, args.error, left_flank,
                                                                 umi, right_flank))
            t2 = threading.Thread(target=run_umi_analysis, args=("umi_analysis", args.fastq, thread, str(3),
                                                                 "a", args.save_path, args.error, left_flank,
                                                                 umi, right_flank))
        else:
            t1 = threading.Thread(target=run_bash_cmd, args=("umi_analysis", args.fastq, thread, str(5)))
            t2 = threading.Thread(target=run_bash_cmd, args=("umi_analysis", args.fastq, thread, str(3)))
        py_threads.append(t1)
        py_threads.append(t2)
    else:
        thread = int(args.thread) // 4
        if args.py_only is True:
            t1 = threading.Thread(target=run_umi_analysis, args=("umi_analysis", args.fastq, thread, str(5),
                                                                 "g", args.save_path, args.error, left_flank,
                                                                 umi, right_flank))
            t2 = threading.Thread(target=run_umi_analysis, args=("umi_analysis", args.fastq, thread, str(3),
                                                                 "a", args.save_path, args.error, left_flank,
                                                                 umi, right_flank))
            t3 = threading.Thread(target=run_umi_analysis, args=("umi_analysis2", args.pe_fastq, thread, str(5),
                                                                 "g", args.save_path, args.error, left_flank,
                                                                 umi, right_flank))
            t4 = threading.Thread(target=run_umi_analysis, args=("umi_analysis2", args.pe_fastq, thread, str(3),
                                                                 "a", args.save_path, args.error, left_flank,
                                                                 umi, right_flank))
        else:
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


def extract_umi_reads(args):
    path = args.save_path
    [save_path, file_name] = extract_mapped_reads.get_name(args)
    outfastq = args.save_path + '/' + file_name + '_umi.fastq'
    cmd = """cat {path}/umi_analysis/5end_UMIs/umi_analysis.05.adapter_with_right_flank.fasta | grep "^>" > {path}/umi_analysis/5.read.with.umi.name.lst
            cat {path}/umi_analysis/3end_UMIs/umi_analysis.05.adapter_with_right_flank.fasta | grep "^>" > {path}/umi_analysis/3.read.with.umi.name.lst
            cat {path}/umi_analysis/5.read.with.umi.name.lst {path}/umi_analysis/3.read.with.umi.name.lst | sort | \
            uniq | cut -d ">" -f2 | cut -d " " -f1  > {path}/umi_analysis/5_3.read.with.umi.name.lst
            seqtk subseq {fastq} {path}/umi_analysis/5_3.read.with.umi.name.lst > {outfile}
            """.format(path=path,
                       fastq=args.fastq,
                       outfile=outfastq)
    subprocess.run(cmd, shell=True, check=True, stdout=subprocess.PIPE,
                   stderr=subprocess.STDOUT).stdout.decode('utf-8').strip()

    if args.pe_fastq is not None:
        outfastq_pe = args.save_path + '/' + file_name + '_umi.pe.fastq'
        cmd = """cat {path}/umi_analysis/5end_UMIs/umi_analysis.05.adapter_with_right_flank.fasta | grep "^>" > {path}/umi_analysis/5.read.with.umi.name.lst
                    cat {path}/umi_analysis/3end_UMIs/umi_analysis.05.adapter_with_right_flank.fasta | grep "^>" > {path}/umi_analysis/3.read.with.umi.name.lst
                    cat {path}/umi_analysis/5.read.with.umi.name.lst {path}/umi_analysis/3.read.with.umi.name.lst | \
                    sort | uniq | cut -d ">" -f2 | cut -d " " -f1  > {path}/umi_analysis/5_3.read.with.umi.name.lst
                    seqtk subseq {fastq} {path}/umi_analysis/5_3.read.with.umi.name.lst > {outfile}
                    """.format(path=path,
                               fastq=args.fastq,
                               outfile=outfastq_pe)
        subprocess.run(cmd, shell=True, check=True, stdout=subprocess.PIPE,
                       stderr=subprocess.STDOUT).stdout.decode('utf-8').strip()
    else:
        outfastq_pe = None
    return outfastq, outfastq_pe


def lst2fastq_variant(args, umi_n5):
    umi_regex5 = '^' + umi_n5 + '$'
    umi_regex5 = re.sub('N', '[ATGC]', umi_regex5)
    umi_regex5 = re.sub('B', '[GTC]', umi_regex5)
    umi_regex5 = re.sub('D', '[GAT]', umi_regex5)
    umi_regex5 = re.sub('H', '[ATC]', umi_regex5)
    umi_regex5 = re.sub('K', '[GT]', umi_regex5)
    umi_regex5 = re.sub('M', '[AC]', umi_regex5)
    umi_regex5 = re.sub('R', '[AG]', umi_regex5)
    umi_regex5 = re.sub('S', '[GC]', umi_regex5)
    umi_regex5 = re.sub('V', '[AGC]', umi_regex5)
    umi_regex5 = re.sub('W', '[AT]', umi_regex5)
    umi_regex5 = re.sub('Y', '[CT]', umi_regex5)

    complement_dict = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    umi_n3 = "".join(complement_dict.get(base, base) for base in reversed(umi_n5))
    umi_regex3 = '^' + umi_n3 + '$'
    # umi_regex3 = re.sub('N', '[ATGC]', umi_regex3)
    # umi_regex3 = re.sub('B', '[GTC]', umi_regex3)
    # umi_regex3 = re.sub('D', '[GAT]', umi_regex3)
    # umi_regex3 = re.sub('H', '[ATC]', umi_regex3)
    # umi_regex3 = re.sub('K', '[GT]', umi_regex3)
    # umi_regex3 = re.sub('M', '[AC]', umi_regex3)
    # umi_regex3 = re.sub('R', '[AG]', umi_regex3)
    # umi_regex3 = re.sub('S', '[GC]', umi_regex3)
    # umi_regex3 = re.sub('V', '[AGC]', umi_regex3)
    # umi_regex3 = re.sub('W', '[AT]', umi_regex3)
    # umi_regex3 = re.sub('Y', '[CT]', umi_regex3)

    umi_regex3 = re.sub('N', '[ATGC]', umi_regex3)
    umi_regex3 = re.sub('B', '[CAG]', umi_regex3)
    umi_regex3 = re.sub('D', '[CTA]', umi_regex3)
    umi_regex3 = re.sub('H', '[TAG]', umi_regex3)
    umi_regex3 = re.sub('K', '[CA]', umi_regex3)
    umi_regex3 = re.sub('M', '[TG]', umi_regex3)
    umi_regex3 = re.sub('R', '[TC]', umi_regex3)
    umi_regex3 = re.sub('S', '[CG]', umi_regex3)
    umi_regex3 = re.sub('V', '[TCG]', umi_regex3)
    umi_regex3 = re.sub('W', '[TA]', umi_regex3)
    umi_regex3 = re.sub('Y', '[GA]', umi_regex3)

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

    [umi_fastq, umi_fastq_pe] = extract_umi_reads(args)

    # multiprocessing module for combining read name from same 5' and 3' umi, and extract reads, and analyze snp
    p = multiprocessing.Pool(args.thread)
    with open(file5, "r") as infile5:
        for seq in infile5:
            seq_splt = seq.strip().split(' ')
            count5 = seq_splt[0]
            umi = seq_splt[1]
            # p.apply_async(use_class_5end, args=(umi, args.save_path, args.fastq, count5, umi_regex5, args.threshold,
            #                                     args.refer, args.bash_thread, args.align_mode, args.pe_fastq))
            p.apply_async(use_class_5end, args=(umi, args.save_path, umi_fastq, count5, umi_regex5, args.threshold,
                                                args.refer, args.bash_thread, args.align_mode, umi_fastq_pe))
    p.close()
    p.join()
    logging.info("%s    5'+3' and 5' end UMI analysis finished" % datetime.now().strftime("%Y-%m-%d %H:%M:%S"))

    # multiprocessing module for processing unused 3' end umi, and extract reads, and analyze snp
    if args.pe_fastq is None:
        p = multiprocessing.Pool(args.thread)
        for file in os.listdir(args.save_path + '/umi_analysis/3end_UMIs/reads_with_same_UMIs/'):
            if file.endswith('.lst'):
                name_lst3 = args.save_path + '/umi_analysis/3end_UMIs/reads_with_same_UMIs/' + file
                umi3 = re.split('[_.]', file)[1]
                count3 = re.split('[_.]', file)[0]
                # p.apply_async(use_class_3end, args=(umi3, args.save_path, args.fastq, count3, umi_regex3,
                #                                     args.threshold, args.refer, args.bash_thread,
                #                                     args.align_mode, args.pe_fastq, name_lst3))
                p.apply_async(use_class_3end, args=(umi3, args.save_path, umi_fastq, count3, umi_regex3,
                                                    args.threshold, args.refer, args.bash_thread,
                                                    args.align_mode, umi_fastq_pe, name_lst3))
        p.close()
        p.join()
        logging.info("%s    The rest 3' end UMI analysis finished" % datetime.now().strftime("%Y-%m-%d %H:%M:%S"))


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
        a = subprocess.run(cmd, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        if a.stdout != b'':
            logging.info("stdout from " + umi_name + " :\n" + a.stdout.decode('utf-8'))
        if a.stderr != b'':
            logging.info("stderr from " + umi_name + " :\n" + a.stderr.decode('utf-8'))

        fastq_file = of

        # Below is for normal data set
        if self.read2 is not None:
            of2 = self.save + '/grouped_reads/' + read_folder + str(read_count) + '_' + umi_name + '.read2.fastq'
            cmd = """seqtk subseq {fastq} {name_lst} > {outfile}""".format(fastq=self.read2,
                                                                           name_lst=name_lst,
                                                                           outfile=of2)
            a = subprocess.run(cmd, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            if a.stdout != b'':
                logging.info("stdout from " + umi_name + " :\n" + a.stdout.decode('utf-8'))
            if a.stderr != b'':
                logging.info("stderr from " + umi_name + " :\n" + a.stderr.decode('utf-8'))
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
                logging.info("\n!!!!!! multiple match for read_name.lst file: %s!!!!!!\n" % name_lst3)
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
    args = arguments.get_argparse()
    check_tools()

    if not os.path.isdir(args.save_path):
        os.mkdir(args.save_path)

    ctime = datetime.now().strftime("%Y%m%d_%H.%M.%S")
    log_file = args.save_path + '/' + ctime + '_vault.log'
    logging.basicConfig(level=logging.DEBUG,
                        format='%(message)s',
                        filename=log_file,
                        filemode='w')
    logging.getLogger().addHandler(logging.StreamHandler(sys.stdout))

    logging.info(" ".join(sys.argv))
    logging.info("\n%s    --------- received arguments ---------" % datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
    logging.info("umi adapter:    %s" % args.umi_adapter)
    logging.info("fastq file:     %s" % args.fastq)
    if args.pe_fastq is not None:
        logging.info("read2 file:    %s" % args.pe_fastq)
    logging.info("refer file:     %s" % args.refer)
    logging.info("error rate:     %s" % args.error)
    logging.info("align mode:     %s" % args.align_mode)
    logging.info("thread:         %s" % args.thread)
    logging.info("save path:      %s\n" % args.save_path)

    fastq_gzip = ""
    pe_fastq_gzip = ""
    tmp_fastq = ""
    tmp_pe_fastq = ""
    if args.fastq.split(".")[-1] == "gz":
        cmd = """gunzip -k {fastq}""".format(fastq=args.fastq)
        subprocess.run(cmd, shell=True, check=True, stdout=subprocess.PIPE)
        args.fastq = ".".join(args.fastq.split(".")[:-1])
        tmp_fastq = args.fastq
        fastq_gzip = True

    if args.pe_fastq is not None:
        if args.pe_fastq.split(".")[-1] == "gz":
            cmd = """gunzip -k {fastq}""".format(fastq=args.pe_fastq)
            subprocess.run(cmd, shell=True, check=True, stdout=subprocess.PIPE)
            args.pe_fastq = ".".join(args.pe_fastq.split(".")[:-1])
            tmp_pe_fastq = args.pe_fastq
            pe_fastq_gzip = True

    if args.minlength != 0 or args.maxlength != 0:
        if args.pe_fastq is not None:
            logging.info("ERROR! Pair-end fastq file provided: %s" % args.refer)
            logging.info("ERROR! NO NEED for read length filter for short reads!")
            logging.info("ERROR! Please run without --minlength and --maxlength")
            sys.exit(1)

        if args.maxlength == 0:
            logging.info("%s    Start filtering reads by length >= %d ..."
                  % (datetime.now().strftime("%Y-%m-%d %H:%M:%S"), args.minlength))
        else:
            logging.info("%s    Start filtering reads by length %d - %d ..."
                  % (datetime.now().strftime("%Y-%m-%d %H:%M:%S"), args.minlength, args.maxlength))

        [save_path, file_name] = extract_mapped_reads.get_name(args)
        args.fastq = read_length_filter.filter_length(args.fastq, args.minlength, args.maxlength, save_path, file_name)

    if args.unmapped_reads is True:
        logging.info("%s    Start extracting mapped reads ..." % datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
        [save_path, file_name] = extract_mapped_reads.get_name(args)
        extract_mapped_reads.pre_alignment(args, save_path, file_name)
        extract_mapped_reads.extract_reads_by_name(args, save_path, file_name)
        args.fastq = save_path + file_name + '.mapped.fastq'
        if args.pe_fastq is not None:
            args.pe_fastq = save_path + file_name + '.read2.mapped.fastq'

    logging.info("%s    Start extracting UMIs from reads ...\n" % datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
    umi_n5 = analyze_umi(args)

    # Below is for the specific Illumina data set
    # combine_pe_umi(args)

    combine_pe_umi(args)
    logging.info("\n%s    Extract UMIs from reads finished" % datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
    extract_name_lst(args)
    logging.info("%s    Group reads by UMIs finished" % datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
    # umi_n5 = 'NNNNNNNNNN'
    logging.info("%s    Start individual UMI group analysis ..." % datetime.now().strftime("%Y-%m-%d %H:%M:%S"))

    lst2fastq_variant(args, umi_n5)
    final_clean_up(args)

    if fastq_gzip is True:
        os.remove(tmp_fastq)

    if pe_fastq_gzip is True:
        os.remove(tmp_pe_fastq)

    logging.info("\n%s    ------------ Jobs done! ------------\n" % datetime.now().strftime("%Y-%m-%d %H:%M:%S"))


if __name__ == '__main__':
    main()
