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
import shlex
import shutil
import arguments
import time
from datetime import datetime
from tools import extract_mapped_reads, filter_sv, umi_group_filter, call_consensus, change_VCF_pos, draw_circos
from tools import read_length_filter
from tools import cutadapt27, sniffles1
from variants_calling import snv_calling, sv_calling
import logging
import gzip


def _is_group_only(args) -> bool:
    """Return True if running in 'grouping only' mode.

    Backward compatible with the old flag name (--only_group_reads).
    """
    return bool(getattr(args, "only_group_reads", False) or getattr(args, "only_group_fastq", False))


def _command_output(cmd, timeout=10):
    try:
        result = subprocess.run(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            check=False,
            timeout=timeout,
        )
    except Exception as e:
        return "ERROR: %s" % e

    output = "\n".join(
        line.strip()
        for line in (result.stdout + "\n" + result.stderr).splitlines()
        if line.strip()
    )
    if output:
        return output
    return "no version output (exit code %s)" % result.returncode


def _first_version_line(output, needles=None):
    lines = [line.strip() for line in output.splitlines() if line.strip()]
    if not lines:
        return output

    if needles:
        for line in lines:
            lower = line.lower()
            if any(needle in lower for needle in needles):
                return line

    return lines[0]


def _tool_path(tool):
    if not tool:
        return "not used"
    return shutil.which(tool) or tool


def _tool_version(tool, version_args, needles=None):
    if not tool:
        return "not used"

    output = _command_output([tool] + list(version_args))
    return _first_version_line(output, needles)


def log_software_versions(args, cutadapt_bin, sniffles_bin=None):
    logging.info("%s    --------- software versions ---------"
                 % datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
    logging.info("VAULT:        %s", getattr(arguments, "vault_version", "unknown"))
    logging.info("python:       %s (%s)", sys.version.split()[0], sys.executable)
    logging.info("cutadapt:     %s (%s)",
                 _tool_version(cutadapt_bin, ["--version"]),
                 _tool_path(cutadapt_bin))
    logging.info("seqtk:        %s (%s)",
                 _tool_version("seqtk", [], ["version:"]),
                 _tool_path("seqtk"))

    if _is_group_only(args):
        logging.info("")
        return

    logging.info("minimap2:     %s (%s)",
                 _tool_version("minimap2", ["--version"]),
                 _tool_path("minimap2"))
    logging.info("samtools:     %s (%s)",
                 _tool_version("samtools", ["--version"], ["samtools"]),
                 _tool_path("samtools"))
    logging.info("bcftools:     %s (%s)",
                 _tool_version("bcftools", ["--version"], ["bcftools"]),
                 _tool_path("bcftools"))

    if sniffles_bin:
        sniffles_args = (
            ["--version"]
            if getattr(args, "sv_caller", "sniffles1") == "sniffles2"
            else ["-h"]
        )
        logging.info("sniffles:     %s (%s)",
                     _tool_version(sniffles_bin, sniffles_args, ["version"]),
                     _tool_path(sniffles_bin))
    logging.info("")


def validate_file(x):
    if not os.path.exists(x):
        raise argparse.ArgumentTypeError("{0} does not exist".format(x))
    return x


def find_all_index(adapt, search='N'):
    index = [x for x in range(adapt.find(search), len(adapt)) if adapt[x] == search]
    return index


def get_base_name(fq_path):
    """
    Return base name without directory and compressed/FASTQ extension.
    Examples:
        /path/nanopore_reads.fastq.gz -> nanopore_reads
        /path/nanopore_reads.300-20000.fastq.gz -> nanopore_reads.300-20000
        /path/nanopore_reads.fastq -> nanopore_reads
    """
    base = os.path.basename(fq_path)
    if base.endswith('.gz'):
        base = base[:-3]
    if base.endswith('.fastq') or base.endswith('.fq'):
        base = base.rsplit('.', 1)[0]
    return base


def filter_fastq_by_length_gz(fastq_in, min_len, max_len, save_path, base_name):
    """
    Filter reads by length and write gzipped FASTQ.
    This replaces read_length_filter.filter_length() in the main pipeline.
    """
    out_path = os.path.join(save_path, base_name + '.fastq.gz')

    if fastq_in.endswith('.gz'):
        fin = gzip.open(fastq_in, 'rt')
    else:
        fin = open(fastq_in, 'r')

    with fin, gzip.open(out_path, 'wt') as fout:
        while True:
            header = fin.readline()
            if not header:
                break
            seq = fin.readline()
            plus = fin.readline()
            qual = fin.readline()
            if not qual:
                break

            L = len(seq.strip())
            if (min_len == 0 or L >= min_len) and (max_len == 0 or L <= max_len):
                fout.write(header)
                fout.write(seq)
                fout.write(plus)
                fout.write(qual)

    return out_path


def analyze_umi(args):
    """
    Run UMI extraction via check_umi.sh.
    This step splits the adapter into left flank / UMI / right flank,
    then runs 5'-UMI and 3'-UMI analysis in parallel.
    """
    # split UMI adapter into left flank, UMI, right flank
    args.umi_adapter = args.umi_adapter.upper()
    n_index = find_all_index(args.umi_adapter)
    left_flank = args.umi_adapter[0:n_index[0]]
    right_flank = args.umi_adapter[n_index[-1] + 1:]
    umi = args.umi_adapter[n_index[0]:n_index[-1] + 1]
    logging.info("left_flank:   %s", left_flank)
    logging.info("umi:          %s", umi)
    logging.info("right_flank:  %s", right_flank)

    if not os.path.isdir(args.save_path):
        os.mkdir(args.save_path)

    def run_bash_cmd(name, fastq, thread, end):
        cmd = """bash {path}/check_umi.sh -n {name} -l {left_flank} -u {umi} -r {right_flank} \
                 -q {fastq} -s {save} -e {error} -t {thread} -{end}""".format(
            path=sys.path[0],
            name=name,
            left_flank=left_flank,
            umi=umi,
            right_flank=right_flank,
            fastq=fastq,
            save=args.save_path,
            error=args.error,
            thread=thread,
            end=end
        )
        # only stderr is printed to screen by the shell script itself
        subprocess.run(cmd, shell=True, check=True, stdout=subprocess.PIPE)

    umi_threads = []
    if args.pe_fastq is None:
        thread = int(args.thread) // 2
        umi_threads.extend([
            threading.Thread(
                target=run_bash_cmd,
                args=("umi_analysis", args.fastq, thread, str(5))
            ),
            threading.Thread(
                target=run_bash_cmd,
                args=("umi_analysis", args.fastq, thread, str(3))
            ),
        ])
    else:
        thread = int(args.thread) // 4
        for name, fastq, end in [
            ("umi_analysis", args.fastq, str(5)),
            ("umi_analysis", args.fastq, str(3)),
            ("umi_analysis2", args.pe_fastq, str(5)),
            ("umi_analysis2", args.pe_fastq, str(3)),
        ]:
            umi_threads.append(
                threading.Thread(
                    target=run_bash_cmd,
                    args=(name, fastq, thread, end)
                )
            )

    for t in umi_threads:
        t.start()

    for t in umi_threads:
        t.join()


def group_read_name(seq, sav, count, name_lst):
    outfile = sav + '/reads_with_same_UMIs/' + count + '_' + seq + '.lst'
    with open(outfile, 'w') as of:
        for name in name_lst:
            of.writelines([name, "\n"])


def _group_read_name_from_task(task):
    seq, sav, count, name_lst, label = task
    group_read_name(seq, sav, count, name_lst)
    return label


def _use_class_5end_from_task(task):
    (
        umi,
        save_path,
        umi_fastq,
        count,
        umi_regex,
        threshold,
        refer,
        bash_thread,
        align_mode,
        umi_fastq_pe,
        label,
    ) = task
    use_class_5end(
        umi,
        save_path,
        umi_fastq,
        count,
        umi_regex,
        threshold,
        refer,
        bash_thread,
        align_mode,
        umi_fastq_pe,
    )
    return label


def _use_class_3end_from_task(task):
    (
        umi,
        save_path,
        umi_fastq,
        count,
        umi_regex,
        threshold,
        refer,
        bash_thread,
        align_mode,
        umi_fastq_pe,
        name_lst3,
        label,
    ) = task
    use_class_3end(
        umi,
        save_path,
        umi_fastq,
        count,
        umi_regex,
        threshold,
        refer,
        bash_thread,
        align_mode,
        umi_fastq_pe,
        name_lst3,
    )
    return label


def _run_limited_process_tasks(task_iterable, worker_func, max_workers, task_name):
    """
    Run process-pool tasks with a bounded submission window.

    At most max_workers tasks are submitted at one time. When one task
    finishes, the next task is submitted.
    """
    worker_count = max(1, int(max_workers))
    task_iter = iter(task_iterable)
    completed = 0
    submitted = 0
    pending = []

    logging.info(
        "%s    Run %s with up to %d worker process(es).",
        datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        task_name,
        worker_count,
    )

    def submit_next(executor):
        nonlocal submitted
        try:
            task = next(task_iter)
        except StopIteration:
            return False
        label = task[-1]
        result = executor.apply_async(worker_func, args=(task,))
        pending.append((result, label))
        submitted += 1
        return True

    pool = multiprocessing.Pool(worker_count)
    try:
        for _ in range(worker_count):
            if not submit_next(pool):
                break

        if not pending:
            logging.info("No tasks found for %s", task_name)
            pool.close()
            pool.join()
            return

        while pending:
            progressed = False
            for result, label in list(pending):
                if not result.ready():
                    continue
                pending.remove((result, label))
                try:
                    result.get()
                except Exception:
                    logging.exception(
                        "%s failed for task: %s after %d completed task(s).",
                        task_name,
                        label,
                        completed,
                    )
                    raise

                completed += 1
                submit_next(pool)
                progressed = True

            if not progressed:
                time.sleep(0.1)

        pool.close()
        pool.join()
    except Exception:
        pool.terminate()
        pool.join()
        raise


def extract_name_lst(args):
    """
    Split the large "UMI -> read names" list into multiple .lst files
    under:
      save_path/umi_analysis/5end_UMIs/reads_with_same_UMIs/
      save_path/umi_analysis/3end_UMIs/reads_with_same_UMIs/
    using a bounded process pool.
    """
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

    def iter_read_name_tasks():
        with open(file5, "r") as infile5:
            for record in infile5:
                record_splt = record.strip().split(' ')
                count5 = record_splt[0]
                umi = record_splt[1]
                read_name_lst = record_splt[2:]
                yield (umi, end5, count5, read_name_lst, count5 + "_" + umi)

        if args.pe_fastq is None:
            with open(file3, "r") as infile3:
                for record in infile3:
                    record_splt = record.strip().split(' ')
                    count3 = record_splt[0]
                    umi = record_splt[1]
                    read_name_lst = record_splt[2:]
                    yield (umi, end3, count3, read_name_lst, count3 + "_" + umi)

    _run_limited_process_tasks(
        iter_read_name_tasks(),
        _group_read_name_from_task,
        args.thread,
        "read-name-list writing",
    )


def extract_umi_reads(args):
    """
    Collect read IDs that contain UMI (5' or 3') and use seqtk subseq
    to extract their reads into <basename>_umi.fastq.gz (and .pe.fastq.gz if PE).
    """
    path = args.save_path
    base_name = get_base_name(args.fastq)
    outfastq = os.path.join(args.save_path, base_name + '_umi.fastq.gz')

    # collect read IDs with UMI, and call seqtk once
    five_fa = os.path.join(path, "umi_analysis", "5end_UMIs", "umi_analysis.05.adapter_with_right_flank.fasta")
    three_fa = os.path.join(path, "umi_analysis", "3end_UMIs", "umi_analysis.05.adapter_with_right_flank.fasta")
    five_name_file = os.path.join(path, "umi_analysis", "5.read.with.umi.name.lst")
    three_name_file = os.path.join(path, "umi_analysis", "3.read.with.umi.name.lst")
    both_name_file = os.path.join(path, "umi_analysis", "5_3.read.with.umi.name.lst")

    def collect_headers(fasta_path, header_path):
        """
        Collect read IDs from FASTA headers and write raw headers to header_path.
        Returns a list of read IDs.
        """
        names = []
        if not os.path.exists(fasta_path):
            return names
        with open(fasta_path, 'r') as fin, open(header_path, 'w') as fout:
            for line in fin:
                if line.startswith('>'):
                    fout.write(line)
                    header = line[1:].strip()
                    read_id = header.split()[0]
                    names.append(read_id)
        return names

    names5 = collect_headers(five_fa, five_name_file)
    names3 = collect_headers(three_fa, three_name_file)

    all_names = sorted(set(names5) | set(names3))
    with open(both_name_file, 'w') as fh:
        for rid in all_names:
            fh.write(rid + "\n")

    cmd = (
        "seqtk subseq "
        f"{shlex.quote(args.fastq)} {shlex.quote(both_name_file)} "
        f"| gzip -c > {shlex.quote(outfastq)}"
    )
    subprocess.run(
        cmd,
        shell=True,
        check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT
    ).stdout.decode('utf-8').strip()

    if args.pe_fastq is not None:
        outfastq_pe = args.save_path + '/' + base_name + '_umi.pe.fastq.gz'
        cmd_pe = (
            "seqtk subseq "
            f"{shlex.quote(args.pe_fastq)} {shlex.quote(both_name_file)} "
            f"| gzip -c > {shlex.quote(outfastq_pe)}"
        )
        subprocess.run(
            cmd_pe,
            shell=True,
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT
        ).stdout.decode('utf-8').strip()
    else:
        outfastq_pe = None

    return outfastq, outfastq_pe


def lst2fastq(args, umi_n5):
    """
    Step1: combine read-name lists by UMI, and generate per-UMI FASTQ files
           under grouped_reads/perfect_umi/*.fastq.gz.
    Variant calling itself is handled later in ProcessUmi.run_variant_on_grouped_fastq().
    """
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

    os.makedirs(os.path.join(args.save_path, 'grouped_reads'), exist_ok=True)
    os.makedirs(os.path.join(args.save_path, 'grouped_reads', 'perfect_umi'), exist_ok=True)
    os.makedirs(
        os.path.join(args.save_path, 'grouped_reads', 'perfect_umi_unanalyzed'),
        exist_ok=True,
    )

    # per_umi_process is only needed for downstream per-UMI variant calling/aggregation.
    # In grouping-only mode we skip creating it to keep the result folder clean.
    if not _is_group_only(args):
        os.makedirs(os.path.join(args.save_path, 'per_umi_process', 'perfect_umi'), exist_ok=True)

    umi_fastq, umi_fastq_pe = extract_umi_reads(args)

    logging.info("%s    Start generating per-UMI FASTQ files (5'+3' and 5'-only groups) ..."
                 % datetime.now().strftime("%Y-%m-%d %H:%M:%S"))

    def iter_5end_fastq_tasks():
        with open(file5, "r") as infile5:
            for seq in infile5:
                seq_splt = seq.strip().split(' ')
                count5 = seq_splt[0]
                umi = seq_splt[1]
                yield (
                    umi,
                    args.save_path,
                    umi_fastq,
                    count5,
                    umi_regex5,
                    args.threshold,
                    args.refer,
                    args.bash_thread,
                    args.align_mode,
                    umi_fastq_pe,
                    count5 + "_" + umi,
                )

    _run_limited_process_tasks(
        iter_5end_fastq_tasks(),
        _use_class_5end_from_task,
        args.thread,
        "5end per-UMI FASTQ generation",
    )
    logging.info("%s    5'+3' and 5' end per-UMI FASTQ generation finished"
                 % datetime.now().strftime("%Y-%m-%d %H:%M:%S"))

    # multiprocessing: generate per-UMI FASTQ for remaining 3' only UMIs
    if args.pe_fastq is None:
        logging.info("%s    Start generating per-UMI FASTQ files for remaining 3' UMIs ..."
                     % datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
        def iter_3end_fastq_tasks():
            end3_dir = args.save_path + '/umi_analysis/3end_UMIs/reads_with_same_UMIs/'
            for file in os.listdir(end3_dir):
                if file.endswith('.lst'):
                    name_lst3 = end3_dir + file
                    umi3 = re.split('[_.]', file)[1]
                    count3 = re.split('[_.]', file)[0]
                    yield (
                        umi3,
                        args.save_path,
                        umi_fastq,
                        count3,
                        umi_regex3,
                        args.threshold,
                        args.refer,
                        args.bash_thread,
                        args.align_mode,
                        umi_fastq_pe,
                        name_lst3,
                        count3 + "_" + umi3,
                    )

        _run_limited_process_tasks(
            iter_3end_fastq_tasks(),
            _use_class_3end_from_task,
            args.thread,
            "3end per-UMI FASTQ generation",
        )
        logging.info("%s    Remaining 3' end per-UMI FASTQ generation finished"
                     % datetime.now().strftime("%Y-%m-%d %H:%M:%S"))


def use_class_5end(a, b, c, d, e, f, g, h, i, j=''):
    process_5end = ProcessUmi(a, b, c, d, e, f, g, h, i, j)
    process_5end.process_5()


def use_class_3end(a, b, c, d, e, f, g, h, i, j='', k=''):
    process_3end = ProcessUmi(a, b, c, d, e, f, g, h, i, j, k)
    process_3end.process_3()


class ProcessUmi:
    """
    Class that handles per-UMI grouping logic and (optionally) per-UMI variant calling.

    In the main VAULT run, this class is used in two phases:

      1) write_per_umi_fastq / process_5 / process_3:
         generate per-UMI FASTQ files in grouped_reads/perfect_umi

      2) run_variant_on_grouped_fastq (static):
         run snv_calling / sv_calling on all grouped FASTQ files in parallel
    """
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
        """
        Combine 5' and 3' read-name list files for a paired UMI, write to grouped_reads/<folder>,
        and rename to use the total read count as prefix.
        """
        name_lst3 = glob.glob(str(name_lst3))[0]

        name_lst = self.save + '/grouped_reads/' + folder + '/' + umi5 + '_' + umi3 + '.lst'

        cmd = 'cat ' + name_lst5 + ' ' + name_lst3 + ' | sort | uniq > ' + name_lst
        subprocess.run(cmd, shell=True, check=True, stdout=subprocess.PIPE)

        read_count = sum(1 for line in open(name_lst, 'r'))

        name_count = self.save + '/grouped_reads/' + folder + '/' + str(
            read_count) + '_' + umi5 + '_' + umi3 + '.lst'

        os.rename(name_lst, name_count)
        os.remove(name_lst5)
        os.remove(name_lst3)
        return [read_count, name_count]

    def keep_unanalyzed_lst(self, name_lst):
        """
        Keep read-name lists that did not pass the analysis threshold outside
        grouped_reads/perfect_umi, so that perfect_umi contains FASTQ files only.
        """
        if not name_lst or not os.path.exists(name_lst):
            return None

        out_dir = os.path.join(
            self.save,
            "grouped_reads",
            "perfect_umi_unanalyzed",
        )
        os.makedirs(out_dir, exist_ok=True)

        out_path = os.path.join(out_dir, os.path.basename(name_lst))
        if os.path.abspath(name_lst) != os.path.abspath(out_path):
            if os.path.exists(out_path):
                os.remove(out_path)
            os.rename(name_lst, out_path)

        return out_path

    def write_per_umi_fastq(self, read_count, umi_name, name_lst, end="_"):
        """
        Step1: given a read-name list and UMI name, generate per-UMI FASTQ (and optional read2).
        The files are written to:
          save/grouped_reads/perfect_umi/
        with names:
          <read_count><end><umi_name>.fastq.gz
          <read_count><end><umi_name>.read2.fastq.gz  (if read2 is available)
        """
        base = f"{read_count}{end}{umi_name}"

        # read1
        fastq1 = os.path.join(self.save, "grouped_reads", "perfect_umi", base + ".fastq.gz")
        cmd1 = (
            "seqtk subseq "
            f"{shlex.quote(self.fastq)} {shlex.quote(name_lst)} "
            f"| gzip -c > {shlex.quote(fastq1)}"
        )
        a1 = subprocess.run(
            cmd1,
            shell=True,
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        if a1.stdout:
            logging.info("stdout from %s:\n%s", umi_name, a1.stdout.decode("utf-8"))
        if a1.stderr:
            logging.info("stderr from %s:\n%s", umi_name, a1.stderr.decode("utf-8"))

        fastq2 = None
        if self.read2:
            fastq2 = os.path.join(self.save, "grouped_reads", "perfect_umi", base + ".read2.fastq.gz")
            cmd2 = (
                "seqtk subseq "
                f"{shlex.quote(self.read2)} {shlex.quote(name_lst)} "
                f"| gzip -c > {shlex.quote(fastq2)}"
            )
            a2 = subprocess.run(
                cmd2,
                shell=True,
                check=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
            )
            if a2.stdout:
                logging.info("stdout(read2) from %s:\n%s", umi_name, a2.stdout.decode("utf-8"))
            if a2.stderr:
                logging.info("stderr(read2) from %s:\n%s", umi_name, a2.stderr.decode("utf-8"))

        os.remove(name_lst)

        return fastq1, fastq2, base

    def process_5(self):
        """
        Handle one 5'-UMI entry: try to find its matching 3'-UMI, merge read lists if possible,
        and generate per-UMI FASTQ files for:
          - 5'+3' UMI pairs
          - 5'-only UMI groups (if no unique 3'-UMI match)
        """
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
                    self.write_per_umi_fastq(read_count, umi_name, name_lst, end="_")
                else:
                    self.keep_unanalyzed_lst(name_lst)

            elif len(glob.glob(str(name_lst3))) > 1:
                logging.info("\n!!!!!! multiple match for read_name.lst file: %s!!!!!!\n", name_lst3)
                if int(self.count) >= int(self.threshold):
                    self.write_per_umi_fastq(self.count, self.umi, name_lst5, end="_5end_")
                else:
                    self.keep_unanalyzed_lst(name_lst5)

            else:
                if int(self.count) >= int(self.threshold):
                    self.write_per_umi_fastq(self.count, self.umi, name_lst5, end="_5end_")
                else:
                    self.keep_unanalyzed_lst(name_lst5)

    def process_3(self):
        """
        Handle one 3'-UMI entry that has not been paired with any 5'-UMI (3'-only UMI groups).
        """
        read_count = sum(1 for line in open(self.name_lst3, 'r'))
        if re.match(str(self.umi_regex), str(self.umi)):
            if read_count >= self.threshold:
                self.write_per_umi_fastq(read_count, self.umi, self.name_lst3, end="_3end_")
            else:
                self.keep_unanalyzed_lst(self.name_lst3)

    # ========= Step2: per-UMI variant calling (parallel) =========
    @staticmethod
    def _variant_worker(
        save_path,
        bash_thread,
        align_mode,
        refer,
        fastq1,
        fastq2,
        base,
        sv_caller,
        sniffles_bin,
        sniffles_min_support,
    ):
        """
        Worker that actually runs snv_calling / sv_calling in multiprocessing.Pool.
        Each worker operates on a single UMI group (one per-UMI FASTQ).
        """
        vcf_dir = os.path.join(save_path, "vcf")
        os.makedirs(vcf_dir, exist_ok=True)

        fastq_file = fastq1 if not fastq2 else f"{fastq1} {fastq2}"

        snv_calling(save_path, bash_thread, align_mode, refer, fastq_file, base)
        if align_mode != "sr":
            sv_calling(
                save_path,
                bash_thread,
                refer,
                fastq_file,
                base,
                sv_caller=sv_caller,
                sniffles_bin=sniffles_bin,
                min_support=sniffles_min_support,
            )

    @staticmethod
    def _variant_worker_from_task(task):
        """
        Run one queued per-UMI variant-calling task and return its group name.

        This wrapper lets Pool.imap_unordered() report completed groups as they
        finish, instead of queuing every task with apply_async() and checking
        errors only after all workers have joined.
        """
        ProcessUmi._variant_worker(*task)
        return task[6]

    @staticmethod
    def run_variant_on_grouped_fastq(args):
        """
        Step2: run per-UMI variant calling on existing per-UMI FASTQ files.

        Input FASTQ:
            grouped_reads/perfect_umi/*.fastq.gz

        For each group <base>, the worker will output:
            per_umi_process/perfect_umi/<base>/vcf/<base>.snv.vcf
            per_umi_process/perfect_umi/<base>/vcf/<base>.sv.vcf (if not SR mode)

        This function uses multiprocessing.Pool with args.thread workers.
        """
        root = os.path.join(args.save_path, "grouped_reads", "perfect_umi")
        if not os.path.isdir(root):
            logging.error("Grouped FASTQ folder not found: %s", root)
            return

        if args.align_mode != "sr" and not getattr(args, "sniffles_bin", None):
            try:
                args.sniffles_bin = sniffles1.resolve_sniffles(
                    getattr(args, "sniffles_path", None),
                    getattr(args, "sv_caller", "sniffles1"),
                )
            except FileNotFoundError as e:
                logging.error(str(e))
                raise

        sv_caller = getattr(args, "sv_caller", "sniffles1")
        sniffles_min_support = getattr(args, "sniffles_min_support", 3)

        perumi_root = os.path.join(args.save_path, "per_umi_process", "perfect_umi")
        os.makedirs(perumi_root, exist_ok=True)

        tasks = []
        for fname in sorted(os.listdir(root)):
            if not fname.endswith(".fastq.gz"):
                continue
            if fname.endswith(".read2.fastq.gz"):
                # read2 is handled together with read1
                continue

            base = fname[:-len(".fastq.gz")]
            fastq1 = os.path.join(root, fname)
            fastq2_path = os.path.join(root, base + ".read2.fastq.gz")
            fastq2 = fastq2_path if os.path.exists(fastq2_path) else None

            save_path = os.path.join(perumi_root, base)
            vcf_dir = os.path.join(save_path, "vcf")
            snv_vcf = os.path.join(vcf_dir, f"{base}.snv.vcf")

            if os.path.exists(snv_vcf):
                logging.info("Skip per-UMI variant calling for %s (VCF already exists)", base)
                continue

            tasks.append((
                save_path,
                args.bash_thread,
                args.align_mode,
                args.refer,
                fastq1,
                fastq2,
                base,
                sv_caller,
                getattr(args, "sniffles_bin", None),
                sniffles_min_support,
            ))

        if not tasks:
            logging.info("No UMI groups found for variant calling under %s", root)
            return

        logging.info("%s    Start per-UMI variant calling on %d grouped FASTQ files ..."
                     % (datetime.now().strftime("%Y-%m-%d %H:%M:%S"), len(tasks)))

        worker_count = max(1, min(int(args.thread), len(tasks)))
        logging.info(
            "%s    Run per-UMI variant calling with %d worker process(es); "
            "each worker uses %d bash thread(s)."
            % (
                datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                worker_count,
                int(args.bash_thread),
            )
        )

        completed = 0
        try:
            with multiprocessing.Pool(worker_count) as pool:
                for base in pool.imap_unordered(
                    ProcessUmi._variant_worker_from_task,
                    tasks,
                    chunksize=1,
                ):
                    completed += 1
                    if not getattr(args, "silence", False):
                        logging.info(
                            "Finished per-UMI variant calling for %s "
                            "(%d/%d)",
                            base,
                            completed,
                            len(tasks),
                        )
        except Exception:
            logging.exception(
                "Per-UMI variant calling failed after %d/%d completed group(s).",
                completed,
                len(tasks),
            )
            raise

        logging.info("%s    Per-UMI variant calling finished"
                     % datetime.now().strftime("%Y-%m-%d %H:%M:%S"))


def final_clean_up(args):
    """
    Aggregate per-UMI VCFs, apply UMI group filtering (optional),
    and produce final SNV/Indel and SV call sets.
    """
    logging.info("%s    Start aggregating per-UMI VCFs and applying final filters ..."
                 % datetime.now().strftime("%Y-%m-%d %H:%M:%S"))

    # process SNV and Indel VCF
    snv_regex = args.save_path + '/per_umi_process/perfect_umi/*/vcf/*snv.vcf'
    all_snv = args.save_path + '/per_umi_process/all_snv_from_perfect_umi.vcf'
    pass_snv = args.save_path + '/per_umi_process/pass_snv_from_perfect_umi.vcf'

    all_snv_pcent = args.save_path + '/per_umi_process/all_snv_from_perfect_umi.pcent.vcf'
    all_snv_pcent_flt = args.save_path + '/per_umi_process/all_snv_from_perfect_umi.pcent.flt.vcf'

    find_vcf_files(snv_regex, all_snv)

    save_dir = args.save_path + '/per_umi_process'
    umi_group_filter.add_pcent(save_dir, all_snv_pcent)

    if args.group_filter is True:
        awk_filter = "awk '{if($3>0.35 && $3<0.65)print}' | awk '{if($7 < 0.4) print}'"
        cmd = """cat {save_dir}/umi_group.flt.summary.txt | {filter} > {save_dir}/wrong.group.summary.txt
                 cat {save_dir}/wrong.group.summary.txt | cut -f1 | sort -k1 > {save_dir}/wrong.group.lst
                 ls {save_dir}/perfect_umi | sort | uniq > {save_dir}/all.group.lst
                 comm -23 {save_dir}/all.group.lst {save_dir}/wrong.group.lst > {save_dir}/pass.group.lst

                 split -l 50 {save_dir}/pass.group.lst {save_dir}/pattern-file.split.
                 cat {save_dir}/all_snv_from_perfect_umi.pcent.vcf | grep "^#" > {save_dir}/all_snv_from_perfect_umi.pcent.rem.vcf
                 for CHUNK in {save_dir}/pattern-file.split.* ; do
                    cat {save_dir}/all_snv_from_perfect_umi.pcent.vcf | grep -v "^#" | grep -f "$CHUNK" \
                        >> {save_dir}/all_snv_from_perfect_umi.pcent.rem.vcf
                 done
                 rm {save_dir}/pattern-file.split.* {save_dir}/all.group.lst

                 cat {save_dir}/all_snv_from_perfect_umi.pcent.rem.vcf | \
                    bcftools filter -s LowQual -e 'QUAL<10 || INFO/DP<3 || INFO/IMF<0.5 || INFO/VARP<{alle_freq} || INFO/ALTC<3' \
                    > {save_dir}/all_snv_from_perfect_umi.pcent.rem.flt.vcf

                 cat {save_dir}/all_snv_from_perfect_umi.pcent.rem.flt.vcf | grep -v LowQual \
                    > {save_dir}/pass_snv_from_perfect_umi.flt.vcf
                 find {save_dir}/perfect_umi/ -name "*coverage.3plus.txt" -print0 | xargs -0 cat > {save_dir}/coverage.3plus.txt
                 """.format(save_dir=save_dir,
                            alle_freq=args.allele_freq,
                            filter=awk_filter)

    else:
        cmd = """cat {all_snv_pcent} | bcftools filter -s LowQual -e 'QUAL<10 || INFO/DP<3 || INFO/IMF<0.5 || \
                    INFO/VARP<{alle_freq} || INFO/ALTC<3' > {fltvcf}
                 rm {all_snv_pcent} {save}/per_umi_process/umi_group.flt.summary.txt
                 cat {fltvcf} | grep -v LowQual > {pass_snv}
                 find {save}/per_umi_process/perfect_umi/ -name "*coverage.3plus.txt" -print0 | xargs -0 cat > {save}/per_umi_process/coverage.3plus.txt
                 """.format(all_snv_pcent=all_snv_pcent,
                            fltvcf=all_snv_pcent_flt,
                            save=args.save_path,
                            alle_freq=args.allele_freq,
                            pass_snv=pass_snv)

    subprocess.run(cmd, shell=True, check=True, stdout=subprocess.PIPE)

    # process SV VCF
    if args.align_mode != 'sr':
        sv_regex = args.save_path + '/per_umi_process/perfect_umi/*/vcf/*sv.vcf'
        all_sv = args.save_path + '/per_umi_process/all_sv_from_perfect_umi.vcf'
        find_vcf_files(sv_regex, all_sv)

        save_path = args.save_path + '/per_umi_process/'
        filter_sv.filter_sv(all_sv, save_path, args.sv_freq)

    logging.info("%s    Final VCF aggregation and filtering finished"
                 % datetime.now().strftime("%Y-%m-%d %H:%M:%S"))


def find_vcf_files(file_regex, allvcf):
    """
    Concatenate per-UMI VCFs that match file_regex into a single VCF (allvcf),
    with an extra "sample_id" column derived from the directory name.
    """
    file_list = glob.glob(str(file_regex))
    if len(file_list) >= 1:
        with open(allvcf, 'w') as of, open(file_list[0], 'r') as header:
            for line in header:
                if line.startswith('##'):
                    if line.startswith("##INFO=<ID=MQ,") and "Type=Integer" in line:
                        line = line.replace("Type=Integer", "Type=Float", 1)
                    of.writelines(line)
                elif line.startswith('#'):
                    of.writelines(line)

            for vcf in file_list:
                vcf_file = open(vcf, 'r')
                sample_id = re.split(r'(perfect_umi/|/vcf)', vcf)[2]
                for line in vcf_file:
                    if not line.startswith('#'):
                        line = re.split(r'\t', line.strip())
                        of.writelines('\t'.join((line[0], line[1], sample_id, '\t'.join(line[3:]))) + '\n')

    else:
        sys.stderr.write("No VCF file detected! Please check %s\n" % file_regex)
        sys.exit(1)


def combine_pe_umi(args):
    """
    For PE data, merge 5'-UMI read-name CSVs from read1 and read2 runs into a single
    umi_analysis.07.1read_count.2umi.3read_name.pe.lst file.
    """
    if args.pe_fastq is not None:
        umi_csv_5end1 = args.save_path + '/umi_analysis/5end_UMIs/umi_analysis.06.1umi.2read_name.csv'
        umi_csv_5end2 = args.save_path + '/umi_analysis2/5end_UMIs/umi_analysis2.06.1umi.2read_name.csv'
        name_lst5 = args.save_path + '/umi_analysis/5end_UMIs/umi_analysis.07.1read_count.2umi.3read_name.pe.lst'

        umi_to_reads = {}

        for csv_path in (umi_csv_5end1, umi_csv_5end2):
            if not os.path.exists(csv_path):
                continue
            with open(csv_path, 'r') as fh:
                for line in fh:
                    line = line.strip()
                    if not line:
                        continue
                    if '\t' in line:
                        umi, read_id = line.split('\t', 1)
                        reads = [read_id]
                    else:
                        parts = line.split(',')
                        umi = parts[0]
                        reads = parts[1:]
                    if umi not in umi_to_reads:
                        umi_to_reads[umi] = set()
                    umi_to_reads[umi].update(reads)

        with open(name_lst5, 'w') as out:
            for umi in sorted(umi_to_reads.keys()):
                reads = sorted(umi_to_reads[umi])
                count = len(reads)
                out.write(f"{count} {umi} {' '.join(reads)}\n")


def auto_summarize_after_main(args, original_fastq):
    """
    Run `vault summarize` automatically after a successful full VAULT run.

    Use the original input FASTQ instead of args.fastq, because args.fastq may
    be replaced by length-filtered or mapped-read intermediates during analysis.
    """
    if _is_group_only(args):
        return

    vault_script = os.path.join(os.path.dirname(os.path.abspath(__file__)), "vault")
    cmd = [
        sys.executable,
        vault_script,
        "summarize",
        "-q",
        original_fastq,
        "-s",
        args.save_path,
        "-r",
        args.refer,
        "-T",
        "0.01",
    ]
    if getattr(args, "unmapped_reads", False):
        cmd.append("--unmapped_reads")

    cmd_text = " ".join(shlex.quote(x) for x in cmd)
    logging.info("%s    Start automatic VAULT summarize ..."
                 % datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
    logging.info("auto summarize command: %s", cmd_text)
    subprocess.run(cmd, check=True)
    logging.info("%s    Automatic VAULT summarize finished"
                 % datetime.now().strftime("%Y-%m-%d %H:%M:%S"))


def check_tools(args, cutadapt_bin, sniffles_bin=None):
    """Sanity check external dependencies.

    VAULT has two very different execution paths:
      1) Full pipeline (alignment + per-UMI variant calling)
      2) Grouping-only (--only_group_reads): extract UMIs and write per-UMI FASTQ files

    In grouping-only mode we intentionally *do not* require alignment/variant tools
    (minimap2/samtools/bcftools/sniffles), because they are never invoked.
    """

    # cutadapt (adapter / UMI extraction)
    if call_consensus.check_tool(cutadapt_bin) is not True:
        logging.warning('ERROR: "cutadapt" is not available / not executable: %s', cutadapt_bin)
        sys.exit(1)

    # seqtk (FASTQ sub-selection to build per-UMI FASTQs)
    if call_consensus.check_tool("seqtk") is not True:
        logging.warning('ERROR: "seqtk" is not available.')
        sys.exit(1)

    # Grouping-only ends here.
    if _is_group_only(args):
        return

    # Alignment + variant calling dependencies
    if call_consensus.check_tool("minimap2") is not True:
        logging.warning('ERROR: "minimap2" is not available.')
        sys.exit(1)

    if call_consensus.check_tool("samtools") is not True:
        logging.warning('ERROR: "samtools" is not available.')
        sys.exit(1)

    if call_consensus.check_tool("bcftools") is not True:
        logging.warning('ERROR: "bcftools" is not available.')
        sys.exit(1)

    # SV caller is only required for long-read modes
    if getattr(args, "align_mode", None) in ("map-ont", "map-pb"):
        if call_consensus.check_tool(sniffles_bin or "sniffles") is not True:
            logging.warning('ERROR: "sniffles" is not available: %s', sniffles_bin)
            sys.exit(1)


def main():
    args = arguments.get_argparse()
    original_fastq = args.fastq

    # Resolve cutadapt path (possibly version 2.7); exit if not found.
    try:
        cutadapt_bin = cutadapt27.resolve_cutadapt(
            getattr(args, "cutadapt_path", None)
        )
    except FileNotFoundError as e:
        sys.stderr.write(str(e) + "\n")
        sys.exit(1)

    sniffles_bin = None
    if (
        getattr(args, "align_mode", None) in ("map-ont", "map-pb")
        and not _is_group_only(args)
    ):
        try:
            sniffles_bin = sniffles1.resolve_sniffles(
                getattr(args, "sniffles_path", None),
                getattr(args, "sv_caller", "sniffles1"),
            )
        except FileNotFoundError as e:
            sys.stderr.write(str(e) + "\n")
            sys.exit(1)

    if not os.path.isdir(args.save_path):
        os.makedirs(args.save_path, exist_ok=True)

    ctime = datetime.now().strftime("%Y%m%d_%H.%M.%S")
    log_file = args.save_path + '/' + ctime + '_vault.log'
    logging.basicConfig(level=logging.DEBUG,
                        format='%(message)s',
                        filename=log_file,
                        filemode='w')
    logging.getLogger().addHandler(logging.StreamHandler(sys.stdout))

    logging.info(" ".join(sys.argv))
    logging.info("\n%s    --------- received arguments ---------" % datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
    logging.info("umi adapter:    %s", args.umi_adapter)
    logging.info("fastq file:     %s", args.fastq)
    if args.pe_fastq is not None:
        logging.info("read2 file:    %s", args.pe_fastq)
    logging.info("refer file:     %s", args.refer)
    logging.info("error rate:     %s", args.error)
    logging.info("align mode:     %s", args.align_mode)
    logging.info("thread:         %s", args.thread)
    logging.info("UMI threshold:  %s", args.threshold)
    logging.info("save path:      %s", args.save_path)
    logging.info("cutadapt bin:   %s", cutadapt_bin)
    if sniffles_bin:
        logging.info("sv caller:      %s", getattr(args, "sv_caller", "sniffles1"))
        logging.info("sniffles bin:   %s", sniffles_bin)
        logging.info("sniffles min support: %s",
                     getattr(args, "sniffles_min_support", 3))
    logging.info("")
    log_software_versions(args, cutadapt_bin, sniffles_bin)

    # export cutadapt path to environment, used by check_umi.sh
    os.environ["VAULT_CUTADAPT_BIN"] = cutadapt_bin
    if sniffles_bin:
        args.sniffles_bin = sniffles_bin
        os.environ["VAULT_SNIFFLES_BIN"] = sniffles_bin

    check_tools(args, cutadapt_bin, sniffles_bin)

    # ====== Step0: read length filter (optional) ======
    if args.minlength != 0 or args.maxlength != 0:
        if args.pe_fastq is not None:
            logging.info("ERROR! Pair-end fastq file provided: %s", args.refer)
            logging.info("ERROR! NO NEED for read length filter for short reads!")
            logging.info("ERROR! Please run without --minlength and --maxlength")
            sys.exit(1)

        if args.maxlength == 0:
            logging.info("%s    Start filtering reads by length >= %d ..."
                         % (datetime.now().strftime("%Y-%m-%d %H:%M:%S"), args.minlength))
        else:
            logging.info("%s    Start filtering reads by length %d - %d ..."
                         % (datetime.now().strftime("%Y-%m-%d %H:%M:%S"), args.minlength, args.maxlength))

        base0 = get_base_name(args.fastq)
        if args.maxlength == 0:
            label = f"{args.minlength}+"
        else:
            label = f"{args.minlength}-{args.maxlength}"
        filtered_basename = f"{base0}.{label}"
        logging.info("Filtered fastq basename will be %s", filtered_basename)

        args.fastq = filter_fastq_by_length_gz(
            args.fastq,
            args.minlength,
            args.maxlength,
            args.save_path,
            filtered_basename
        )

    # ====== Step1: extract mapped reads (optional) ======
    if args.unmapped_reads is True:
        logging.info("%s    Start extracting mapped reads ..."
                     % datetime.now().strftime("%Y-%m-%d %H:%M:%S"))

        save_path = args.save_path + '/'
        file_name = get_base_name(args.fastq)

        extract_mapped_reads.pre_alignment(args, save_path, file_name)
        extract_mapped_reads.extract_reads_by_name(args, save_path, file_name)

        mapped_fastq = save_path + file_name + '.mapped.fastq'
        if os.path.exists(mapped_fastq):
            subprocess.run(["gzip", "-f", mapped_fastq], check=True,
                           stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            args.fastq = mapped_fastq + '.gz'
        else:
            args.fastq = mapped_fastq

        if args.pe_fastq is not None:
            mapped_fastq_pe = save_path + file_name + '.read2.mapped.fastq'
            if os.path.exists(mapped_fastq_pe):
                subprocess.run(["gzip", "-f", mapped_fastq_pe], check=True,
                               stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                args.pe_fastq = mapped_fastq_pe + '.gz'
            else:
                args.pe_fastq = mapped_fastq_pe

    # ====== Step2: UMI extraction ======
    logging.info("%s    Start extracting UMIs from reads ..."
                 % datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
    analyze_umi(args)

    combine_pe_umi(args)
    logging.info("\n%s    Extract UMIs from reads finished"
                 % datetime.now().strftime("%Y-%m-%d %H:%M:%S"))

    # ====== Step3: group read names by UMI ======
    logging.info("%s    Start grouping reads by UMIs ..."
                 % datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
    extract_name_lst(args)
    logging.info("%s    Group reads by UMIs finished"
                 % datetime.now().strftime("%Y-%m-%d %H:%M:%S"))

    # ====== Step4: generate per-UMI FASTQ files ======
    logging.info("%s    Start generating per-UMI FASTQ groups ..."
                 % datetime.now().strftime("%Y-%m-%d %H:%M:%S"))

    args.umi_adapter = args.umi_adapter.upper()
    n_index = find_all_index(args.umi_adapter)
    umi_n5 = args.umi_adapter[n_index[0]:n_index[-1] + 1]

    lst2fastq(args, umi_n5)

    logging.info(
        "%s    Per-UMI FASTQ grouping finished "
        "(FASTQ: grouped_reads/perfect_umi/*.fastq.gz; "
        "unanalyzed LST: grouped_reads/perfect_umi_unanalyzed/*.lst)"
        % datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    )

    # ====== Step5: per-UMI variant calling + final aggregation ======
    if _is_group_only(args):
        logging.info("%s    --only_group_reads is set: skip per-UMI variant calling and final aggregation."
                     % datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
        logging.info("\n%s    ------------ Jobs done (grouping only)! ------------\n"
                     % datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
    else:
        logging.info("%s    Start per-UMI variant calling ..."
                     % datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
        ProcessUmi.run_variant_on_grouped_fastq(args)

        final_clean_up(args)

        logging.info("\n%s    ------------ Jobs done! ------------\n"
                     % datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
        auto_summarize_after_main(args, original_fastq)


if __name__ == '__main__':
    main()
