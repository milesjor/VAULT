#!/usr/bin/env python3

import subprocess
import numpy as np
import pandas as pd
import re
import io
import sys
import logging
from datetime import datetime
import os
import glob
import gzip


def check_read_number(checked_fastq):
    paths = sorted(glob.glob(checked_fastq)) or [checked_fastq]
    total = 0

    for path in paths:
        opener = gzip.open if path.endswith(".gz") else open
        name = path[:-3] if path.endswith(".gz") else path

        if name.endswith((".fastq", ".fq")):
            lines = 0
            with opener(path, "rt") as infile:
                for _ in infile:
                    lines += 1
            total += lines // 4

        elif name.endswith((".fasta", ".fa")):
            with opener(path, "rt") as infile:
                total += sum(1 for line in infile if line.startswith(">"))

    return str(total)


def read_group_list(group_file):
    with open(group_file, "r") as infile:
        return [line.strip() for line in infile if line.strip()]


def format_number(value):
    if isinstance(value, float):
        text = ("%.6f" % value).rstrip("0").rstrip(".")
        return text if text else "0"
    return str(value)


def summarize_numeric_values(values):
    if not values:
        raise RuntimeError(
            "No passed UMI coverage values were found; "
            "check per_umi_process/pass.group.lst and coverage files."
        )

    values = sorted(values)
    count = len(values)
    total = sum(values)
    avg = total / count
    if count % 2 == 1:
        median = values[count // 2]
    else:
        median = (values[count // 2 - 1] + values[count // 2]) / 2

    return ",".join([
        format_number(avg),
        format_number(median),
        format_number(values[0]),
        format_number(values[-1]),
    ])


def summarize_snv_load_values(path, refer_seq_length, min_coverage=1):
    load_file = os.path.join(
        path,
        "per_umi_process",
        "summary",
        "pass.1molecule_2snv_3coverage_4divide.txt",
    )
    values = []
    with open(load_file, "r") as infile:
        for line in infile:
            fields = line.strip().split()
            if len(fields) < 4:
                continue
            try:
                coverage = int(float(fields[2]))
                load = float(fields[3])
            except ValueError:
                continue
            if coverage >= min_coverage:
                values.append(load * float(refer_seq_length))
    return summarize_numeric_values(values)


def coverage_file_for_group(path, group):
    return os.path.join(
        path,
        "per_umi_process",
        "perfect_umi",
        group,
        group + ".coverage.txt",
    )


def count_covered_bases(coverage_file, min_depth=3):
    count = 0
    with open(coverage_file, "r") as infile:
        for line in infile:
            fields = line.strip().split()
            if len(fields) < 3:
                continue
            try:
                depth = int(float(fields[2]))
            except ValueError:
                continue
            if depth >= min_depth:
                count += 1
    return count


def write_pass_position_coverage(path, pass_groups):
    pos_counts = {}
    summary_dir = os.path.join(path, "per_umi_process", "summary")
    pass_cov_dir = os.path.join(summary_dir, "pass_coverage")
    os.makedirs(pass_cov_dir, exist_ok=True)

    for old_file in glob.glob(os.path.join(pass_cov_dir, "*.coverage.txt")):
        os.remove(old_file)

    for group in pass_groups:
        coverage_file = coverage_file_for_group(path, group)
        if not os.path.exists(coverage_file):
            logging.warning("Coverage file not found for group %s: %s", group, coverage_file)
            continue

        copied_file = os.path.join(pass_cov_dir, os.path.basename(coverage_file))
        with open(coverage_file, "r") as infile, open(copied_file, "w") as copied:
            for line in infile:
                copied.write(line)
                fields = line.strip().split()
                if len(fields) < 3:
                    continue
                try:
                    pos = int(fields[1])
                    depth = int(float(fields[2]))
                except ValueError:
                    continue
                if depth >= 3:
                    pos_counts[pos] = pos_counts.get(pos, 0) + 1

    out_file = os.path.join(summary_dir, "pass.coverage.3plus.pos.molecule.txt")
    with open(out_file, "w") as outfile:
        for pos in sorted(pos_counts):
            outfile.write("%s %s\n" % (pos_counts[pos], pos))


def write_coverage_length_summary(path, refer_seq_length):
    summary_dir = os.path.join(path, "per_umi_process", "summary")
    all_cov_dir = os.path.join(summary_dir, "all_coverage")
    os.makedirs(all_cov_dir, exist_ok=True)

    for old_file in glob.glob(os.path.join(all_cov_dir, "*.coverage.txt")):
        os.remove(old_file)

    all_groups = [
        name for name in sorted(os.listdir(os.path.join(path, "per_umi_process", "perfect_umi")))
        if os.path.isdir(os.path.join(path, "per_umi_process", "perfect_umi", name))
    ]
    pass_groups = set(read_group_list(os.path.join(path, "per_umi_process", "pass.group.lst")))

    all_lengths = {}
    for group in all_groups:
        coverage_file = coverage_file_for_group(path, group)
        if not os.path.exists(coverage_file):
            logging.warning("Coverage file not found for group %s: %s", group, coverage_file)
            continue

        copied_file = os.path.join(all_cov_dir, os.path.basename(coverage_file))
        with open(coverage_file, "r") as infile, open(copied_file, "w") as copied:
            copied.writelines(infile.readlines())
        all_lengths[group] = count_covered_bases(coverage_file)

    all_length_file = os.path.join(summary_dir, "all.coverage.3plus.length.txt")
    pass_length_file = os.path.join(summary_dir, "pass.coverage.3plus.length.txt")
    pass_95p_file = os.path.join(summary_dir, "pass.coverage.3plus.length.95p.txt")
    pass_95p_group_file = os.path.join(summary_dir, "pass.95p.group.lst")

    with open(all_length_file, "w") as outfile:
        for group in sorted(all_lengths):
            outfile.write("%s %s\n" % (group, all_lengths[group]))

    pass_values = []
    pass_95p_groups = []
    threshold_95p = float(refer_seq_length) * 0.95
    with open(pass_length_file, "w") as outfile, open(pass_95p_file, "w") as out95:
        for group in sorted(pass_groups):
            if group not in all_lengths:
                logging.warning("Passed group has no coverage length: %s", group)
                continue
            value = all_lengths[group]
            pass_values.append(value)
            line = "%s %s\n" % (group, value)
            outfile.write(line)
            if value >= threshold_95p:
                out95.write(line)
                pass_95p_groups.append(group)

    with open(pass_95p_group_file, "w") as outfile:
        for group in sorted(pass_95p_groups):
            outfile.write(group + "\n")

    write_95p_snv_vcf(path, set(pass_95p_groups))

    return summarize_numeric_values(pass_values), len(pass_95p_groups)


def write_95p_snv_vcf(path, pass_95p_groups):
    input_vcf = os.path.join(path, "per_umi_process", "pass_snv_from_perfect_umi.flt.vcf")
    output_vcf = os.path.join(
        path,
        "per_umi_process",
        "summary",
        "pass_snv_from_perfect_umi.flt.95p.vcf",
    )

    with open(input_vcf, "r") as infile, open(output_vcf, "w") as outfile:
        for line in infile:
            if line.startswith("#"):
                outfile.write(line)
                continue
            fields = line.split("\t")
            if len(fields) >= 3 and fields[2] in pass_95p_groups:
                outfile.write(line)


def snv_load_per_group(path):
    snv_file = path + '/per_umi_process/pass_snv_from_perfect_umi.flt.vcf'
    coverage_file = path + '/per_umi_process/summary/all.coverage.3plus.length.txt'
    snv_load_distribution = path + '/per_umi_process/summary/pass.1molecule_2snv_3coverage_4divide.txt'

    with open(coverage_file, 'r') as infile1, open(snv_load_distribution, 'w') as outfile:
        dir = {}
        for line in infile1:
            k, v = line.strip().split(' ')
            dir[k.strip()] = v.strip()

        awk = """awk '{print $1,$2}' """
        snv_count = path + '/per_umi_process/summary/pass.snv.count.per.group.txt'
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


def vcf_2_circos_like_file(path, refer_name):
    snv_file = path + '/per_umi_process/pass_snv_from_perfect_umi.flt.vcf'
    circos_like_file = path + '/per_umi_process/summary/pass.snv.2pos.4count.txt'
    snv_vcf_no_indel = path + '/per_umi_process/summary/pass_snv_from_perfect_umi.flt.snvONLY.vcf'

    if not os.path.exists(snv_file):
        logging.warning("SNV VCF file not found: %s", snv_file)
        return

    # 先过滤掉 INDEL，只保留 SNV
    with open(snv_file, 'r') as infile, open(snv_vcf_no_indel, 'w') as outfile:
        for line in infile:
            if line.startswith('#'):
                outfile.writelines(line)
            else:
                info_field = re.split(r'\t', line)[7]
                if not info_field.startswith('INDEL'):
                    outfile.writelines(line)

    col_names = [
        "CHROM",
        "POS",
        "ID",
        "REF",
        "ALT",
        "QUAL",
        "FILTER",
        "INFO",
        "FORMAT",
        "SAMPLE"
    ]

    df = pd.read_table(
        snv_vcf_no_indel,
        comment="#",
        header=None,
        sep="\t",
        names=col_names,
        dtype={"CHROM": str, "POS": int}
    )

    if df.empty:
        logging.warning("No SNVs in file: %s", snv_vcf_no_indel)
        return

    count_pos = (
        df.groupby(["CHROM", "POS"])
        .size()
        .reset_index(name="count")
    )

    sort_count = count_pos.sort_values(["CHROM", "POS"]).reset_index(drop=True)
    sort_count.insert(0, 'chr', refer_name)
    snv_df = sort_count[['chr', 'POS', 'POS', 'count']]

    np.savetxt(circos_like_file, snv_df, fmt='%s')


def count_to_freq(path):
    molecule_coverage_file = path + '/per_umi_process/summary/pass.coverage.3plus.pos.molecule.txt'
    circos_like_file = path + '/per_umi_process/summary/pass.snv.2pos.4count.txt'
    circos_like_file_vaf = path + '/per_umi_process/summary/pass.snv.2pos.4count.5VAF.txt'

    with open(molecule_coverage_file, 'r') as infile1, open(circos_like_file, 'r') as infile2, \
            open(circos_like_file_vaf, 'w') as outfile:
        dir = {}
        for line in infile1:
            v, k = line.strip().split(' ')
            dir[k.strip()] = v.strip()

        for line2 in infile2:
            line2_split = line2.split(' ')
            count = line2_split[3]
            pos = line2_split[1]
            depth = dir[pos]
            freq = round(int(count) / int(depth), 6)
            new_line2 = line2_split[0] + ' ' + line2_split[1] + ' ' + line2_split[2] + ' ' + line2_split[3].strip() + \
                        ' ' + str(freq) + '\n'
            outfile.writelines(new_line2)


def read_usage(args):
    path = args.save_path
    raw_fastq = args.fastq
    refer_seq = args.refer
    somatic_snv_vaf_threshold = args.somatic_VAF
    unmapped_reads = args.unmapped_reads
    sv_file_name = "all_sv_from_perfect_umi.filtered.[0-9].[0-9].vcf"

    raw_read_number = ""
    used_read_number = ""
    reads_with_umi = ""
    detected_molecule_number = ""
    detected_passed_molecule_number = ""
    refer_seq_length = ""
    covered_region_of_molecule = ""
    p95_coverage_molecule = ""

    molecule_with_snv = ""
    total_snv_number = ""
    unique_snv_number = ""
    snv_number_per_molecule = ""
    total_somatic_snv_number = ""
    unique_somatic_snv_number = ""
    somatic_snv_load_per_mbp = ""

    molecule_with_sv = ""
    total_sv_number = ""
    unique_sv_number = ""
    molecule_with_deletion = ""
    total_deletion = ""
    molecule_with_insertion = ""
    total_insertion = ""
    molecule_with_inversion = ""
    total_inversion = ""
    molecule_with_duplication = ""
    total_duplication = ""

    raw_read_number = check_read_number(raw_fastq)
    logging.info("raw_read_number is: %s" % raw_read_number)

    # 新实现：兼容 .mapped.fastq 和 .mapped.fastq.gz，并在缺失时安全返回 0
    if unmapped_reads is True:
        pattern_plain = os.path.join(path, '*.mapped.fastq')
        pattern_gz = os.path.join(path, '*.mapped.fastq.gz')

        if glob.glob(pattern_plain):
            used_read_number = check_read_number(pattern_plain)
        elif glob.glob(pattern_gz):
            used_read_number = check_read_number(pattern_gz)
        else:
            logging.warning("No *.mapped.fastq[.gz] found under %s", path)
            used_read_number = "0"
    else:
        used_read_number = "N/A"

    logging.info("used_read_number is: %s" % used_read_number)

    cmd = """cat {path}/umi_analysis/5end_UMIs/umi_analysis.05.adapter_with_right_flank.fasta | grep "^>" > {path}/umi_analysis/5.read.with.umi.name.lst
            cat {path}/umi_analysis/3end_UMIs/umi_analysis.05.adapter_with_right_flank.fasta | grep "^>" > {path}/umi_analysis/3.read.with.umi.name.lst
            cat {path}/umi_analysis/5.read.with.umi.name.lst {path}/umi_analysis/3.read.with.umi.name.lst | sort | uniq | wc -l
            """.format(path=path)
    reads_with_umi = subprocess.run(cmd, shell=True, check=True, stdout=subprocess.PIPE,
                                    stderr=subprocess.STDOUT).stdout.decode('utf-8').strip()
    logging.info("reads_with_umi is: %s" % reads_with_umi)

    cmd = """ls {path}/per_umi_process/perfect_umi | wc -l""".format(path=path)
    detected_molecule_number = subprocess.run(cmd, shell=True, check=True, stdout=subprocess.PIPE,
                                              stderr=subprocess.STDOUT).stdout.decode('utf-8').strip()
    logging.info("detected_molecule_number is: %s" % detected_molecule_number)

    cmd = """wc -l {path}/per_umi_process/pass.group.lst """.format(path=path)
    detected_passed_molecule_number = subprocess.run(cmd, shell=True, check=True, stdout=subprocess.PIPE,
                                                     stderr=subprocess.STDOUT).stdout.decode('utf-8').strip().split(
        " ")[0]
    logging.info("detected_passed_molecule_number is: %s" % detected_passed_molecule_number)

    pass_groups = read_group_list(os.path.join(path, "per_umi_process", "pass.group.lst"))
    write_pass_position_coverage(path, pass_groups)

    awk_length = """awk 'BEGIN{n=0} {if($1 !~ /^>/)n+=length($1)} END{print n}' """
    cmd = """cat {refer_seq} | {awk_length}""".format(refer_seq=refer_seq,
                                                      awk_length=awk_length)
    refer_seq_length = subprocess.run(cmd, shell=True, check=True, stdout=subprocess.PIPE,
                                      stderr=subprocess.STDOUT).stdout.decode('utf-8').strip()
    logging.info("refer_seq_length is: %s" % refer_seq_length)

    summary_dir = os.path.join(path, "per_umi_process", "summary")
    for fname in [
        "all.coverage.3plus.length.txt",
        "pass.coverage.3plus.length.txt",
        "pass.coverage.3plus.length.95p.txt",
        "pass.95p.group.lst",
        "pass_snv_from_perfect_umi.flt.95p.vcf",
    ]:
        fpath = os.path.join(summary_dir, fname)
        if os.path.exists(fpath):
            os.remove(fpath)

    covered_region_of_molecule, p95_coverage_molecule = write_coverage_length_summary(
        path,
        int(refer_seq_length),
    )
    logging.info("covered_region_of_molecule(avg,median,min,max) is: %s" % covered_region_of_molecule)

    logging.info("p95_coverage_molecule is: %s" % p95_coverage_molecule)

    cmd = """cat {path}/per_umi_process/pass_snv_from_perfect_umi.flt.vcf | grep -v "^#" | grep -v "INDEL" | cut -f3 | sort | uniq | wc -l""".format(
        path=path)
    molecule_with_snv = subprocess.run(cmd, shell=True, check=True, stdout=subprocess.PIPE,
                                       stderr=subprocess.STDOUT).stdout.decode('utf-8').strip()
    logging.info("molecule_with_snv is: %s" % molecule_with_snv)

    cmd = """cat {path}/per_umi_process/pass_snv_from_perfect_umi.flt.vcf | grep -v "^#" | grep -v "INDEL" | wc -l""".format(
        path=path)
    total_snv_number = subprocess.run(cmd, shell=True, check=True, stdout=subprocess.PIPE,
                                      stderr=subprocess.STDOUT).stdout.decode('utf-8').strip()
    logging.info("total_snv_number is: %s" % total_snv_number)

    cmd = """cat {path}/per_umi_process/pass_snv_from_perfect_umi.flt.vcf | grep -v "^#" | grep -v "INDEL" | cut -f2,4,5 | sort | uniq | wc -l""".format(
        path=path)
    unique_snv_number = subprocess.run(cmd, shell=True, check=True, stdout=subprocess.PIPE,
                                       stderr=subprocess.STDOUT).stdout.decode('utf-8').strip()
    logging.info("unique_snv_number is: %s" % unique_snv_number)

    snv_load_per_group(path)
    snv_number_per_molecule = summarize_snv_load_values(path, refer_seq_length)
    logging.info("normalized_snv_number_per_molecule(avg,median,min,max) is: %s" % snv_number_per_molecule)

    refer_name = ""
    with open(refer_seq, 'r') as infile:
        for line in infile:
            if line.startswith(">"):
                refer_name = line.split(">")[1].split()[0]

    vcf_2_circos_like_file(path, refer_name)
    count_to_freq(path)

    awk7 = """awk '{if($5<'$threshold')print $2}' """
    awk8 = """awk '{if($2=='$pos')print}' """
    cmd = """threshold={somatic_threshold}
        somatic_snv_pos=$(cat {path}/per_umi_process/summary/pass.snv.2pos.4count.5VAF.txt | {awk7} )
        cat {path}/per_umi_process/summary/pass_snv_from_perfect_umi.flt.snvONLY.vcf | grep "^#" > \
            {path}/per_umi_process/summary/pass_snv_from_perfect_umi.flt.snvONLY.somatic.vcf
        for pos in $somatic_snv_pos; do cat {path}/per_umi_process/summary/pass_snv_from_perfect_umi.flt.snvONLY.vcf | \
            grep -v "^#" | {awk8} >> {path}/per_umi_process/summary/pass_snv_from_perfect_umi.flt.snvONLY.somatic.vcf ; done
        cat {path}/per_umi_process/summary/pass_snv_from_perfect_umi.flt.snvONLY.somatic.vcf | grep -v "^#" | cut -f 2,4,5 | sort \
            | uniq -c > {path}/per_umi_process/summary/pass.somatic.1count.2snv.txt
        unique_somatic_snv_number=0
        total_somatic_snv_number=0
        unique_somatic_snv_number=$(cat {path}/per_umi_process/summary/pass_snv_from_perfect_umi.flt.snvONLY.somatic.vcf | grep -v "^#" | cut -f2,4,5 | sort | uniq | wc -l)
        total_somatic_snv_number=$(cat {path}/per_umi_process/summary/pass_snv_from_perfect_umi.flt.snvONLY.somatic.vcf | grep -v "^#" | wc -l)
        echo $total_somatic_snv_number $unique_somatic_snv_number
        """.format(somatic_threshold=somatic_snv_vaf_threshold,
                   path=path,
                   awk7=awk7,
                   awk8=awk8)
    total_somatic_snv_number, unique_somatic_snv_number = subprocess.run(cmd, shell=True, check=True,
                                                                         stdout=subprocess.PIPE,
                                                                         stderr=subprocess.STDOUT).stdout.decode(
        'utf-8').strip().split()
    logging.info("total_somatic_snv_number is: %s" % total_somatic_snv_number)
    logging.info("unique_somatic_snv_number is: %s" % unique_somatic_snv_number)

    awk9 = """awk 'BEGIN{n=0}{n+=$2}END{print n}'  """
    awk10 = """awk '{print $1/$2*1000000}' """
    cmd = """total_somatic_snv=$(cat {path}/per_umi_process/summary/pass_snv_from_perfect_umi.flt.snvONLY.somatic.vcf | grep -v "^#" | wc -l)
            total_coverage=$(cat {path}/per_umi_process/summary/pass.coverage.3plus.length.txt | {awk9})
            snv_load=$(echo $total_somatic_snv $total_coverage | {awk10})
            echo $snv_load""".format(path=path, awk9=awk9, awk10=awk10)
    somatic_snv_load_per_mbp = subprocess.run(cmd, shell=True, check=True, stdout=subprocess.PIPE,
                                              stderr=subprocess.STDOUT).stdout.decode('utf-8').strip()
    logging.info("somatic_snv_load_per_Mbp is: %s" % somatic_snv_load_per_mbp)

    # ===== 下面是 SV 统计部分 =====
    awk11 = """awk '{print $1,$2,$3,$4}'"""
    awk12 = """awk -F "\t|SVTYPE=|;SUPTYPE=|;SVLEN=|;STRANDS" '{print $2,$9,$11}'  """
    awk13 = """awk -F "\t|SVTYPE=|;SUPTYPE=|;SVLEN=|;STRANDS" '{print $3,$2,$9,$11}'  """
    awk14 = """awk '{print $1}' | sort | uniq | wc -l """
    sv_pattern = os.path.join(path, "per_umi_process", sv_file_name)
    sv_files = glob.glob(sv_pattern)

    if len(sv_files) == 0:
        logging.warning("No SV VCF files found matching pattern: %s", sv_pattern)
        molecule_with_sv = total_sv_number = unique_sv_number = "0"
        molecule_with_deletion = total_deletion = "0"
        molecule_with_insertion = total_insertion = "0"
        molecule_with_inversion = total_inversion = "0"
        molecule_with_duplication = total_duplication = "0"
    else:
        cmd = """
                molecule_with_sv=0
                total_sv_number=0
                unique_sv_number=0
                molecule_with_deletion=0
                total_deletion=0
                molecule_with_insertion=0
                total_insertion=0
                molecule_with_inversion=0
                total_inversion=0
                molecule_with_duplication=0
                total_duplication=0

                molecule_with_sv=$(cat {path}/per_umi_process/{sv_file_name} | grep -v "^#" | cut -f3 \
                    | sort | uniq | wc -l)
                total_sv_number=$(cat {path}/per_umi_process/{sv_file_name} | grep -v "^#" | wc -l)
                unique_sv_number=$(cat {path}/per_umi_process/{sv_file_name} | grep -v "^#" | \
                    {awk12} | sort | uniq | wc -l)
                cat {path}/per_umi_process/{sv_file_name} | grep -v "^#" | \
                    {awk12} | sort | uniq -c | sort -k3 | {awk11} > {path}/per_umi_process/summary/all_sv.1count.2pos.3type.4length.txt
                molecule_with_deletion=$(cat {path}/per_umi_process/{sv_file_name} | grep -v "^#" | {awk13} | grep "DEL" | {awk14})
                total_deletion=$(cat {path}/per_umi_process/{sv_file_name} | grep -v "^#" | {awk13} | grep "DEL" | wc -l)

                molecule_with_insertion=$(cat {path}/per_umi_process/{sv_file_name} | grep -v "^#" | {awk13} | grep "INS" | {awk14})
                total_insertion=$(cat {path}/per_umi_process/{sv_file_name} | grep -v "^#" | {awk13} | grep "INS" | wc -l)

                molecule_with_inversion=$(cat {path}/per_umi_process/{sv_file_name} | grep -v "^#" | {awk13} | grep "INV" | {awk14})
                total_inversion=$(cat {path}/per_umi_process/{sv_file_name} | grep -v "^#" | {awk13} | grep "INV" | wc -l)

                molecule_with_duplication=$(cat {path}/per_umi_process/{sv_file_name} | grep -v "^#" | {awk13} | grep "DUP" | {awk14})
                total_duplication=$(cat {path}/per_umi_process/{sv_file_name} | grep -v "^#" | {awk13} | grep "DUP" | wc -l)
                echo $molecule_with_sv $total_sv_number $unique_sv_number $molecule_with_deletion $total_deletion \
                    $molecule_with_insertion $total_insertion $molecule_with_inversion $total_inversion \
                    $molecule_with_duplication $total_duplication
                """.format(path=path, awk11=awk11, awk12=awk12, awk13=awk13, awk14=awk14,
                           pass_group_number=detected_passed_molecule_number,
                           sv_file_name=sv_file_name)

        sv_out = subprocess.run(
            cmd,
            shell=True,
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT
        ).stdout.decode('utf-8').strip()

        # 只取最后一行（前面可能有提示信息）
        last_line = sv_out.splitlines()[-1] if sv_out else ""
        parts = last_line.split()

        if len(parts) != 11:
            logging.warning("Unexpected SV summary output (expect 11 fields, got %d): %s", len(parts), last_line)
            molecule_with_sv = total_sv_number = unique_sv_number = "0"
            molecule_with_deletion = total_deletion = "0"
            molecule_with_insertion = total_insertion = "0"
            molecule_with_inversion = total_inversion = "0"
            molecule_with_duplication = total_duplication = "0"
        else:
            molecule_with_sv, total_sv_number, unique_sv_number, molecule_with_deletion, total_deletion, \
                molecule_with_insertion, total_insertion, molecule_with_inversion, total_inversion, \
                molecule_with_duplication, total_duplication = parts

    # 打印 SV 统计信息，注意避免除以 0
    if int(detected_passed_molecule_number) > 0:
        logging.info("molecule_with_sv is: %s/%.2f%%" % (
            molecule_with_sv,
            100 * float(int(molecule_with_sv) / int(detected_passed_molecule_number))
        ))
        logging.info("molecule_with_deletion is: %s/%.2f%%" % (
            molecule_with_deletion,
            100 * float(int(molecule_with_deletion) / int(detected_passed_molecule_number))
        ))
        logging.info("molecule_with_insertion is: %s/%.2f%%" % (
            molecule_with_insertion,
            100 * float(int(molecule_with_insertion) / int(detected_passed_molecule_number))
        ))
        logging.info("molecule_with_inversion is: %s/%.2f%%" % (
            molecule_with_inversion,
            100 * float(int(molecule_with_inversion) / int(detected_passed_molecule_number))
        ))
        logging.info("molecule_with_duplication is: %s/%.2f%%" % (
            molecule_with_duplication,
            100 * float(int(molecule_with_duplication) / int(detected_passed_molecule_number))
        ))
    else:
        logging.info("molecule_with_sv is: %s/NA" % molecule_with_sv)
        logging.info("molecule_with_deletion is: %s/NA" % molecule_with_deletion)
        logging.info("molecule_with_insertion is: %s/NA" % molecule_with_insertion)
        logging.info("molecule_with_inversion is: %s/NA" % molecule_with_inversion)
        logging.info("molecule_with_duplication is: %s/NA" % molecule_with_duplication)

    logging.info("total_sv_number is: %s" % total_sv_number)
    logging.info("unique_sv_number is: %s" % unique_sv_number)
    logging.info("total_deletion is: %s" % total_deletion)
    logging.info("total_insertion is: %s" % total_insertion)
    logging.info("total_inversion is: %s" % total_inversion)
    logging.info("total_duplication is: %s" % total_duplication)


def summarize_main(args):
    per_umi_dir = os.path.join(args.save_path, "per_umi_process")
    if not os.path.isdir(per_umi_dir):
        raise FileNotFoundError(
            "VAULT per_umi_process folder not found: %s. "
            "Please check -s/--save_path points to a finished VAULT result folder."
            % per_umi_dir
        )

    os.makedirs(os.path.join(per_umi_dir, "summary"), exist_ok=True)

    ctime = datetime.now().strftime("%Y%m%d_%H.%M.%S")
    log_file = args.save_path + '/per_umi_process/summary/' + ctime + '_vault_summary.txt'
    logging.basicConfig(level=logging.DEBUG,
                        format='%(message)s',
                        filename=log_file,
                        filemode='w')
    logging.getLogger().addHandler(logging.StreamHandler(sys.stdout))

    logging.info(" ".join(sys.argv))
    logging.info("\n-------------- received arguments -------------")
    logging.info("-vault path:             %s" % args.save_path)
    logging.info("-raw fastq file:         %s" % args.fastq)
    logging.info("-reference sequence:     %s" % args.refer)
    logging.info("-somatic VAF threshold:  %s" % args.somatic_VAF)
    logging.info("-unmapped reads?:        %s" % args.unmapped_reads)
    logging.info("-result saved in:        %s" % (args.save_path + "/per_umi_process/summary/\n"))

    read_usage(args)

    logging.info("\n------------ Jobs done! ------------\n")

    # 把 summarize 的输出内容 append 到主 vault log 末尾
    try:
        main_logs = glob.glob(os.path.join(args.save_path, "*_vault.log"))
        if main_logs:
            # 按文件名排序（YYYYmmdd_HH.MM.SS_vault.log），最后一个就是最新的一次 vault run
            main_logs_sorted = sorted(main_logs)
            main_log = main_logs_sorted[-1]

            with open(log_file, "r") as sf, open(main_log, "a") as mf:
                mf.write("\n\n===== VAULT summarize output (%s) =====\n" % ctime)
                mf.writelines(sf.readlines())
        else:
            logging.warning("No main vault log (*_vault.log) found under %s, skip appending summarize output.",
                            args.save_path)
    except Exception as e:
        logging.warning("Failed to append summarize output to main VAULT log: %s", e)

    sys.exit(0)
