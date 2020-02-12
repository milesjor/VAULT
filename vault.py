#!/home/bic/anaconda3/envs/bioconda/bin/python
# This is for analyzing UMI labeled reads
# Chongwei 20190619
# Email: chongwei.bi@kaust.edu.sa

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
from tools import extract_mapped_reads, filter_sv
from variants_calling import snp_calling, sv_calling


def get_argparse():
    args = argparse.ArgumentParser(description='This is a python script for analyzing UMI labeled reads.'
                                               'require: cutadapt, minimap2, seqtk, samtools, sniffles',
                                   epilog='Refer github https://github.com/milesjor/VAULT for more detail')
    req = args.add_argument_group(title='Required options')
    req.add_argument('-u', '--umi_adapter', type=str, required=True, help='UMI sequence, automatically detect NNNATGCNNN as UMI')
    req.add_argument('-s', '--save_path', type=str, required=True, help='path/to/save/')
    req.add_argument('-r', '--refer', type=validate_file, required=True, help='path/to/ref.fa')
    req.add_argument('-q', '--fastq', type=validate_file, required=True, help='path/to/reads.fastq, or fastq.gz')

    opt = args.add_argument_group(title='Optional options')
    opt.add_argument('-e', '--error', type=str, default='0.11', help='error tolerate rate for umi analysis (read error rate) [0.11]')
    opt.add_argument('-t', '--thread', type=int, default='5', help='thread/process number [5]')
    opt.add_argument('-T', '--threshold', type=int, default='5', help='Threshold of read number for snp analysis [5]')
    opt.add_argument('-b', '--bash_thread', type=int, default='1', help='Thread for running bash cmd [1]')
    opt.add_argument('-p', '--pe_fastq', type=str, help='read2.fastq for illumina pair-end sequencing')
    opt.add_argument('-a', '--align_mode', type=str, default='map-ont',
                     choices=['sr', 'map-ont', 'map-pb'], help='parameter in alignment, minimap2 -ax [sr|map-ont|map-pb]')
    opt.add_argument('--unmapped_reads', action='store_true', help='extract mapped reads before UMI analysis')
    opt.add_argument('-v', '--version', action='version', version=vault_version, help='show the current version')

    args = args.parse_args()
    if args.pe_fastq is not None:
        args.align_mode = 'sr'

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

    n_index = find_all_index(args.umi_adapter)
    left_flank = args.umi_adapter[0:n_index[0]]
    right_flank = args.umi_adapter[n_index[-1]+1:]
    umi = args.umi_adapter[n_index[0]:n_index[-1]+1]
    print("\nleft_flank:  ", left_flank)
    print("umi:         ", umi)
    print("right_flank:  %s\n" % right_flank)

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
        subprocess.run(cmd, shell=True, check=True, stdout=subprocess.PIPE)

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

    print("%s    Extract read name list finished" % datetime.now().strftime("%Y-%m-%d %H:%M:%S"))


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
    print("%s    5 end snp analysis finished" % datetime.now().strftime("%Y-%m-%d %H:%M:%S"))

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
        print("%s    3 end snp analysis finished" % datetime.now().strftime("%Y-%m-%d %H:%M:%S"))


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
        if self.read2 is not None:
            of2 = self.save + '/grouped_reads/' + read_folder + str(read_count) + '_' + umi_name + '.read2.fastq'
            cmd = """seqtk subseq {fastq} {name_lst} > {outfile}""".format(fastq=self.read2,
                                                                           name_lst=name_lst,
                                                                           outfile=of2)
            subprocess.run(cmd, shell=True, check=True, stdout=subprocess.PIPE)
            fastq_file = of + ' ' + of2
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
    snp_regex = args.save_path + '/snp/perfect_umi/*/vcf/*snp.vcf'
    all_snp = args.save_path + '/snp/all_snp_from_perfect_umi.vcf'
    pass_snp = args.save_path + '/snp/pass_snp_from_perfect_umi.vcf'
    find_vcf_files(snp_regex, all_snp, pass_snp)

    if args.align_mode != 'sr':
        sv_regex = args.save_path + '/snp/perfect_umi/*/vcf/*sv.vcf'
        all_sv = args.save_path + '/snp/all_sv_from_perfect_umi.vcf'
        find_vcf_files(sv_regex, all_sv)

        flt_sv_regex = args.save_path + '/snp/perfect_umi/*/vcf/*sv.flt.vcf'
        flt_sv = args.save_path + '/snp/flt_sv_from_perfect_umi.vcf'
        find_vcf_files(flt_sv_regex, flt_sv)
        save_path = args.save_path + '/snp/'
        filter_sv.filter_sv(all_sv, save_path, float('0.5'))


def find_vcf_files(file_regex, allvcf, passvcf='no'):
    file_list = glob.glob(str(file_regex))
    if len(file_list) >= 1:
        with open(allvcf, 'w') as of, open(file_list[0], 'r') as header:
            for line in header:
                if line.startswith('#'):
                    of.writelines(line)

            for vcf in file_list:
                vcf_file = open(vcf, 'r')
                sample_id = re.split(r'(perfect_umi/|/vcf)', vcf)[2]
                for line in vcf_file:
                    if not line.startswith('#'):
                        line = re.split(r'\t', line)
                        of.writelines('\t'.join((line[0], line[1], sample_id, '\t'.join(line[3:]))))

        if passvcf != 'no':
            with open(passvcf, 'w') as ps_out:
                with open(allvcf, 'r') as infile:
                    for line in infile:
                        if line.startswith('#'):
                            ps_out.writelines(line)
                        elif re.split(r'\t', line)[6] == 'PASS':
                            ps_out.writelines(line)

    else:
        print("No VCF file detected! Please check %s\n" % file_regex)


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


def main():
    args = get_argparse()
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

    umi_n5 = analyze_umi(args)
    combine_pe_umi(args)
    print("\n%s    UMI analysis finished" % datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
    extract_name_lst(args)
    # umi_n5 = 'NNNNNNNNNN'
    lst2fastq_variant(args, umi_n5)
    final_clean_up(args)
    print("\n%s    ------------ Jobs done! ------------\n" % datetime.now().strftime("%Y-%m-%d %H:%M:%S"))


if __name__ == '__main__':
    main()