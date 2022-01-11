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


def check_read_number(checked_fastq):
    read_number = ""
    if checked_fastq.split(".")[-1] == "gz":
        if checked_fastq.split(".")[-2] == "fastq" or checked_fastq.split(".")[-2] == "fq":
            cmd = """gunzip -c {fastq} | grep -e "^@" | wc -l""".format(fastq=checked_fastq)
            read_number = subprocess.run(cmd, shell=True, check=True, stdout=subprocess.PIPE,
                                             stderr=subprocess.STDOUT).stdout.decode('utf-8').strip()

        elif checked_fastq.split(".")[-2] == "fasta" or checked_fastq.split(".")[-2] == "fa":
            cmd = """gunzip -c {fastq} | grep -e "^>" | wc -l""".format(fastq=checked_fastq)
            read_number = subprocess.run(cmd, shell=True, check=True, stdout=subprocess.PIPE,
                                             stderr=subprocess.STDOUT).stdout.decode('utf-8').strip()
    else:
        if checked_fastq.split(".")[-1] == "fastq" or checked_fastq.split(".")[-1] == "fq":
            cmd = """cat {fastq} | grep -e "^@" | wc -l""".format(fastq=checked_fastq)
            read_number = subprocess.run(cmd, shell=True, check=True, stdout=subprocess.PIPE,
                                             stderr=subprocess.STDOUT).stdout.decode('utf-8').strip()

        elif checked_fastq.split(".")[-1] == "fasta" or checked_fastq.split(".")[-1] == "fa":
            cmd = """cat {fastq} | grep -e "^>" | wc -l""".format(fastq=checked_fastq)
            read_number = subprocess.run(cmd, shell=True, check=True, stdout=subprocess.PIPE,
                                             stderr=subprocess.STDOUT).stdout.decode('utf-8').strip()
    return read_number


def snv_load_per_group(path):
    snv_file = path + '/snp/pass_snp_from_perfect_umi.flt.vcf'
    coverage_file = path + '/snp/summary/all.coverage.3plus.length.txt'
    snv_load_distribution = path + '/snp/summary/pass.1molecule_2snv_3coverage_4divide.txt'

    with open(coverage_file, 'r') as infile1, open(snv_load_distribution, 'w') as outfile:
        dir = {}
        for line in infile1:
            k, v = line.strip().split(' ')
            dir[k.strip()] = v.strip()

        awk = """awk '{print $1,$2}' """
        snv_count = path + '/snp/summary/pass.snv.count.per.group.txt'
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
    snv_file = path + '/snp/pass_snp_from_perfect_umi.flt.vcf'
    circos_like_file = path + '/snp/summary/pass.snv.2pos.4count.txt'
    snv_vcf_no_indel = path + '/snp/summary/pass_snp_from_perfect_umi.flt.snvONLY.vcf'

    with open(snv_file, 'r') as infile, open(snv_vcf_no_indel, 'w') as outfile:
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

    pd_vcf = read_vcf(snv_vcf_no_indel)

    # print(pd_vcf['POS'])
    # print(pd_vcf.POS.nunique())
    # print(pd_vcf.POS.count())
    # print(pd_vcf.POS.value_counts())
    # print(pd_vcf.POS.value_counts().reset_index())

    count_pos = pd_vcf.POS.value_counts().reset_index().rename(columns={'POS': 'count', 'index': 'POS'})
    sort_count = count_pos.sort_values('POS').reset_index(drop=True)
    sort_count.insert(0, 'chr', refer_name)
    snp_df = sort_count[['chr', 'POS', 'POS', 'count']]

    np.savetxt(circos_like_file, snp_df, fmt='%s')


def count_to_freq(path):
    molecule_coverage_file = path + '/snp/summary/pass.coverage.3plus.pos.molecule.txt'
    circos_like_file = path + '/snp/summary/pass.snv.2pos.4count.txt'
    circos_like_file_vaf = path + '/snp/summary/pass.snv.2pos.4count.5VAF.txt'

    with open(molecule_coverage_file, 'r') as infile1, open(circos_like_file, 'r') as infile2, \
            open(circos_like_file_vaf, 'w') as outfile:
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
            new_line2 = line2_split[0] + ' ' + line2_split[1] + ' ' + line2_split[2] + ' ' + line2_split[3].strip() + \
                        ' ' + str(freq) + '\n'
            outfile.writelines(new_line2)


def read_usage(args):
    # input
    # path = "/Users/bic/Desktop/codes/github/VAULT_local/example/result/"
    # raw_fastq = "/Users/bic/Desktop/codes/github/VAULT_local/example/nanopore_reads.fastq.gz"
    # refer_seq = "/Users/bic/Desktop/codes/github/VAULT_local/example/reference.fa"
    # somatic_snv_vaf_threshold = "0.2"
    # unmapped_reads = True

    path = args.save_path
    raw_fastq = args.fastq
    refer_seq = args.refer
    somatic_snv_vaf_threshold = args.somatic_VAF
    unmapped_reads = args.unmapped_reads
    sv_file_name = "all_sv_from_perfect_umi.filtered.[0-9]*\.[0-9].vcf"

    # TBD
    raw_read_number = ""
    used_read_number = ""  # mapped and length filter passed reads
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

    if unmapped_reads is True:
        mapped_fastq = path + '/*.mapped.fastq'
        used_read_number = check_read_number(mapped_fastq)
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

    cmd = """ls {path}/snp/perfect_umi | wc -l""".format(path=path)
    detected_molecule_number = subprocess.run(cmd, shell=True, check=True, stdout=subprocess.PIPE,
                                                     stderr=subprocess.STDOUT).stdout.decode('utf-8').strip()
    logging.info("detected_molecule_number is: %s" % detected_molecule_number)

    cmd = """wc -l {path}/snp/pass.group.lst """.format(path=path)
    detected_passed_molecule_number = subprocess.run(cmd, shell=True, check=True, stdout=subprocess.PIPE,
                                                     stderr=subprocess.STDOUT).stdout.decode('utf-8').strip().split(" ")[0]
    logging.info("detected_passed_molecule_number is: %s" % detected_passed_molecule_number)

    # summary of covered_region_of_molecule is for passed group only
    awk = """awk '{if($3>=3)print $2}' """
    cmd = """# mkdir {path}/snp/summary
            mkdir {path}/snp/summary/pass_coverage
            for line in $(cat {path}/snp/pass.group.lst); \
                do cp {path}/snp/perfect_umi/$line/$line.coverage.txt {path}/snp/summary/pass_coverage ; done
            for file in {path}/snp/summary/pass_coverage/*.coverage.txt; \
                do cat $file | {awk} >> {path}/snp/summary/pass.coverage.3plus.pos.txt.tmp ; done
            cat {path}/snp/summary/pass.coverage.3plus.pos.txt.tmp | sort -n | uniq -c > {path}/snp/summary/pass.coverage.3plus.pos.molecule.txt
            rm {path}/snp/summary/pass.coverage.3plus.pos.txt.tmp
            """.format(path=path,
                       awk=awk)
    subprocess.run(cmd, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)

    awk_length = """awk 'BEGIN{n=0} {if($1 !~ /^>/)n+=length($1)} END{print n}' """
    cmd = """cat {refer_seq} | {awk_length}""".format(refer_seq=refer_seq,
                                                      awk_length=awk_length)
    refer_seq_length = subprocess.run(cmd, shell=True, check=True, stdout=subprocess.PIPE,
                                      stderr=subprocess.STDOUT).stdout.decode('utf-8').strip()
    logging.info("refer_seq_length is: %s" % refer_seq_length)

    awk2 = """awk 'BEGIN {c = 0;sum = 0;} $1 ~ /^(\-)?[0-9]*(\.[0-9]*)?$/ {a[c++] = $1; sum += $1;} \
                END {ave = sum / c;if( (c % 2) == 1 ) {median = a[ int(c/2) ];} else {median = ( a[c/2] + a[c/2-1] ) / 2;} \
                OFS="\t"; print sum, c, ave, median, a[0], a[c-1];}'"""
    awk3 = """awk -F '/' '{print $NF}' | awk -F '.' '{print $1}' """
    awk4 = """awk 'BEGIN{n=0}{if($3>=3)n+=1}END{print "'$name'",n}' """
    cmd = """# mkdir {path}/snp/summary
            mkdir {path}/snp/summary/all_coverage
            ls {path}/snp/perfect_umi/ > {path}/snp/summary/all.group.lst
            for line in $(cat {path}/snp/summary/all.group.lst); do cp {path}/snp/perfect_umi/$line/$line.coverage.txt \
                {path}/snp/summary/all_coverage ; done
            for file in {path}/snp/summary/all_coverage/*.coverage.txt; do name=$(echo $file | {awk3}) ; \
                cat $file | {awk4} >> {path}/snp/summary/all.coverage.3plus.length.txt ; done
            # less {path}/snp/summary/coverage.3plus.txt | cut -d " " -f2 | sort -n | {awk2}
            split -l 100 {path}/snp/pass.group.lst {path}/snp/summary/pattern-file.split.
            for CHUNK in {path}/snp/summary/pattern-file.split.* ; do cat {path}/snp/summary/all.coverage.3plus.length.txt | \
                grep -f "$CHUNK"  >> {path}/snp/summary/pass.coverage.3plus.length.txt ; done
            rm {path}/snp/summary/pattern-file.split.*
            less {path}/snp/summary/pass.coverage.3plus.length.txt | cut -d " " -f2 | sort -n | {awk2}
            """.format(path=path,
                       awk2=awk2,
                       awk3=awk3,
                       awk4=awk4)
    covered_region_of_molecule = ",".join(subprocess.run(cmd, shell=True, check=True, stdout=subprocess.PIPE,
                                                         stderr=subprocess.STDOUT).stdout.decode('utf-8').strip().split("\t")[2:6])
    logging.info("covered_region_of_molecule(avg,median,min,max) is: %s" % covered_region_of_molecule)

    awk5 = """awk '{if($2>='$refer_length')print}' """
    cmd = """refer_length={refer_length}
            cat {path}/snp/summary/pass.coverage.3plus.length.txt | {awk5} > {path}/snp/summary/pass.coverage.3plus.length.95p.txt
            cat {path}/snp/summary/pass.coverage.3plus.length.95p.txt | cut -d " " -f1 > {path}/snp/summary/pass.95p.group.lst
            split -l 100 {path}/snp/summary/pass.95p.group.lst {path}/snp/summary/pattern-file.split.
            cat {path}/snp/pass_snp_from_perfect_umi.flt.vcf | grep "^#" > {path}/snp/summary/pass_snp_from_perfect_umi.flt.95p.vcf
            for CHUNK in {path}/snp/summary/pattern-file.split.* ; do cat {path}/snp/pass_snp_from_perfect_umi.flt.vcf \
                | grep -v "^#" | grep -f "$CHUNK"  >> {path}/snp/summary/pass_snp_from_perfect_umi.flt.95p.vcf ; done
            rm {path}/snp/summary/pattern-file.split.*
            wc -l {path}/snp/summary/pass.95p.group.lst""".format(path=path,
                                                                  awk5=awk5,
                                                                  refer_length=int(refer_seq_length)*0.95)
    subprocess.run(cmd, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)

    cmd = """wc -l {path}/snp/summary/pass.95p.group.lst""".format(path=path)
    p95_coverage_molecule = subprocess.run(cmd, shell=True, check=True, stdout=subprocess.PIPE,
                                           stderr=subprocess.STDOUT).stdout.decode('utf-8').strip().split(" ")[0]

    logging.info("p95_coverage_molecule is: %s" % p95_coverage_molecule)

    cmd = """cat {path}/snp/pass_snp_from_perfect_umi.flt.vcf | grep -v "^#" | grep -v "INDEL" | cut -f3 | sort | uniq | wc -l""".format(path=path)
    molecule_with_snv = subprocess.run(cmd, shell=True, check=True, stdout=subprocess.PIPE,
                                      stderr=subprocess.STDOUT).stdout.decode('utf-8').strip()
    logging.info("molecule_with_snv is: %s" % molecule_with_snv)

    cmd = """cat {path}/snp/pass_snp_from_perfect_umi.flt.vcf | grep -v "^#" | grep -v "INDEL" | wc -l""".format(path=path)
    total_snv_number = subprocess.run(cmd, shell=True, check=True, stdout=subprocess.PIPE,
                                       stderr=subprocess.STDOUT).stdout.decode('utf-8').strip()
    logging.info("total_snv_number is: %s" % total_snv_number)

    cmd = """cat {path}/snp/pass_snp_from_perfect_umi.flt.vcf | grep -v "^#" | grep -v "INDEL" | cut -f2,4,5 | sort | uniq | wc -l""".format(path=path)
    unique_snv_number = subprocess.run(cmd, shell=True, check=True, stdout=subprocess.PIPE,
                                stderr=subprocess.STDOUT).stdout.decode('utf-8').strip()
    logging.info("unique_snv_number is: %s" % unique_snv_number)

    snv_load_per_group(path)
    awk6 = """awk '{if($3>='$length_filter')print $4*'$refer_length'}' """
    cmd = """length_filter={length_filter}
            refer_length={refer_length}
            cat {path}/snp/summary/pass.1molecule_2snv_3coverage_4divide.txt | {awk6} | sort -n | {awk2}
            """.format(length_filter=1,
                       refer_length=refer_seq_length,
                       path=path,
                       awk6=awk6,
                       awk2=awk2)
    snv_number_per_molecule = ",".join(subprocess.run(cmd, shell=True, check=True, stdout=subprocess.PIPE,
                                                      stderr=subprocess.STDOUT).stdout.decode('utf-8').strip().split("\t")[2:6])
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
        somatic_snv_pos=$(cat {path}/snp/summary/pass.snv.2pos.4count.5VAF.txt | {awk7} )
        cat {path}/snp/summary/pass_snp_from_perfect_umi.flt.snvONLY.vcf | grep "^#" > \
            {path}/snp/summary/pass_snp_from_perfect_umi.flt.snvONLY.somatic.vcf
        for pos in $somatic_snv_pos; do cat {path}/snp/summary/pass_snp_from_perfect_umi.flt.snvONLY.vcf | \
            grep -v "^#" | {awk8} >> {path}/snp/summary/pass_snp_from_perfect_umi.flt.snvONLY.somatic.vcf ; done
        cat {path}/snp/summary/pass_snp_from_perfect_umi.flt.snvONLY.somatic.vcf | grep -v "^#" | cut -f 2,4,5 | sort \
            | uniq -c > {path}/snp/summary/pass.somatic.1count.2snv.txt
        unique_somatic_snv_number=0
        total_somatic_snv_number=0
        unique_somatic_snv_number=$(cat {path}/snp/summary/pass_snp_from_perfect_umi.flt.snvONLY.somatic.vcf | grep -v "^#" | cut -f2,4,5 | sort | uniq | wc -l)
        total_somatic_snv_number=$(cat {path}/snp/summary/pass_snp_from_perfect_umi.flt.snvONLY.somatic.vcf | grep -v "^#" | wc -l)
        echo $total_somatic_snv_number $unique_somatic_snv_number
        """.format(somatic_threshold=somatic_snv_vaf_threshold,
                   path=path,
                   awk7=awk7,
                   awk8=awk8)
    total_somatic_snv_number, unique_somatic_snv_number = subprocess.run(cmd, shell=True, check=True, stdout=subprocess.PIPE,
                                                 stderr=subprocess.STDOUT).stdout.decode('utf-8').strip().split()
    logging.info("total_somatic_snv_number is: %s" % total_somatic_snv_number)
    logging.info("unique_somatic_snv_number is: %s" % unique_somatic_snv_number)

    awk9 = """awk 'BEGIN{n=0}{n+=$2}END{print n}'  """
    awk10 = """awk '{print $1/$2*1000000}' """
    cmd = """total_somatic_snv=$(cat {path}/snp/summary/pass_snp_from_perfect_umi.flt.snvONLY.somatic.vcf | grep -v "^#" | wc -l)
            total_coverage=$(cat {path}/snp/summary/pass.coverage.3plus.length.txt | {awk9})
            snv_load=$(echo $total_somatic_snv $total_coverage | {awk10})
            echo $snv_load""".format(path=path, awk9=awk9, awk10=awk10)
    somatic_snv_load_per_mbp = subprocess.run(cmd, shell=True, check=True, stdout=subprocess.PIPE,
                                              stderr=subprocess.STDOUT).stdout.decode('utf-8').strip()
    logging.info("somatic_snv_load_per_Mbp is: %s" % somatic_snv_load_per_mbp)

    awk11 = """awk '{print $1,$2,$3,$4}'"""
    awk12 = """awk -F "\t|SVTYPE=|;SUPTYPE=|;SVLEN=|;STRANDS" '{print $2,$9,$11}'  """
    awk13 = """awk -F "\t|SVTYPE=|;SUPTYPE=|;SVLEN=|;STRANDS" '{print $3,$2,$9,$11}'  """
    awk14 = """awk '{print $1}' | sort | uniq | wc -l """
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
            
            molecule_with_sv=$(cat {path}/snp/{sv_file_name} | grep -v "^#" | cut -f3 \
                | sort | uniq | wc -l)
            total_sv_number=$(cat {path}/snp/{sv_file_name} | grep -v "^#" | wc -l)
            unique_sv_number=$(cat {path}/snp/{sv_file_name} | grep -v "^#" | \
                {awk12} | sort | uniq | wc -l)
            cat {path}/snp/{sv_file_name} | grep -v "^#" | \
                {awk12} | sort | uniq -c | sort -k3 | {awk11} > {path}/snp/summary/all_sv.1count.2pos.3type.4length.txt
            molecule_with_deletion=$(cat {path}/snp/{sv_file_name} | grep -v "^#" | {awk13} | grep "DEL" | {awk14})
            total_deletion=$(cat {path}/snp/{sv_file_name} | grep -v "^#" | {awk13} | grep "DEL" | wc -l)
            
            molecule_with_insertion=$(cat {path}/snp/{sv_file_name} | grep -v "^#" | {awk13} | grep "INS" | {awk14})
            total_insertion=$(cat {path}/snp/{sv_file_name} | grep -v "^#" | {awk13} | grep "INS" | wc -l)
            
            molecule_with_inversion=$(cat {path}/snp/{sv_file_name} | grep -v "^#" | {awk13} | grep "INV" | {awk14})
            total_inversion=$(cat {path}/snp/{sv_file_name} | grep -v "^#" | {awk13} | grep "INV" | wc -l)
            
            molecule_with_duplication=$(cat {path}/snp/{sv_file_name} | grep -v "^#" | {awk13} | grep "DUP" | {awk14})
            total_duplication=$(cat {path}/snp/{sv_file_name} | grep -v "^#" | {awk13} | grep "DUP" | wc -l)
            echo $molecule_with_sv $total_sv_number $unique_sv_number $molecule_with_deletion $total_deletion \
                $molecule_with_insertion $total_insertion $molecule_with_inversion $total_inversion \
                $molecule_with_duplication $total_duplication
            """.format(path=path, awk11=awk11, awk12=awk12, awk13=awk13, awk14=awk14,
                       pass_group_number=detected_passed_molecule_number,
                       sv_file_name=sv_file_name)

    molecule_with_sv, total_sv_number, unique_sv_number, molecule_with_deletion, total_deletion, \
    molecule_with_insertion, total_insertion, molecule_with_inversion, total_inversion, \
    molecule_with_duplication, total_duplication = subprocess.run(cmd, shell=True, check=True, stdout=subprocess.PIPE,
                                                                  stderr=subprocess.STDOUT).stdout.decode('utf-8').strip().split()
    logging.info("molecule_with_sv is: %s/%.2f%%" % (molecule_with_sv, 100*float(int(molecule_with_sv)/int(detected_passed_molecule_number))))
    logging.info("total_sv_number is: %s" % total_sv_number)
    logging.info("unique_sv_number is: %s" % unique_sv_number)
    logging.info("molecule_with_deletion is: %s/%.2f%%" % (molecule_with_deletion, 100*float(int(molecule_with_deletion)/int(detected_passed_molecule_number))))
    logging.info("total_deletion is: %s" % total_deletion)
    logging.info("molecule_with_insertion is: %s/%.2f%%" % (molecule_with_insertion, 100*float(int(molecule_with_insertion)/int(detected_passed_molecule_number))))
    logging.info("total_insertion is: %s" % total_insertion)
    logging.info("molecule_with_inversion is: %s/%.2f%%" % (molecule_with_inversion, 100*float(int(molecule_with_inversion)/int(detected_passed_molecule_number))))
    logging.info("total_inversion is: %s" % total_inversion)
    logging.info("molecule_with_duplication is: %s/%.2f%%" % (molecule_with_duplication, 100*float(int(molecule_with_duplication)/int(detected_passed_molecule_number))))
    logging.info("total_duplication is: %s" % total_duplication)


def summarize_main(args):
    if not os.path.isdir(args.save_path + '/snp/summary/'):
        os.mkdir(args.save_path + '/snp/summary/')

    ctime = datetime.now().strftime("%Y%m%d_%H.%M.%S")
    log_file = args.save_path + '/snp/summary/' + ctime + '_vault_summary.txt'
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
    logging.info("-result saved in:        %s" % args.save_path + "/snp/summary/\n")

    read_usage(args)

    logging.info("\n------------ Jobs done! ------------\n")
    sys.exit(0)
