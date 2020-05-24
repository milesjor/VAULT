import re
import statistics
import subprocess


def add_pcent(save_path, outvcf):
    vcf = save_path + '/all_snp_from_perfect_umi.vcf'
    outsummary = save_path + '/umi_group.flt.summary.txt'

    with open(vcf, 'r') as infile, open(outvcf, 'w') as outfile1, open(outsummary, 'w') as outfile2:
        process_umi = ''
        record_count = 0
        outfile2.writelines('\t'.join(['process_umi', 'record_count', 'var_read_pcet', 'total_read_count',
                                       'reads_combine', 'average_var_pcent', 'var_stdev', 'average_dp4_pcent',
                                       'var_percent_list', 'dp4_percent_list\n']))
        for line in infile:
            if line.startswith('##'):
                outfile1.writelines(line)
            elif line.startswith('#'):
                outfile1.writelines("""##INFO=<ID=VARP,Number=1,Type=Float,Description="Variant Percentage: Variant allele number divided by total allele number">\n""")
                outfile1.writelines("""##INFO=<ID=DP4P,Number=1,Type=Float,Description="DP4 Percentage: Dp4 read number divided by UMI group read number">\n""")
                outfile1.writelines(line)

            else:
                line = re.split(r'\t', line)
                umi_name = line[2]
                read_count = int(re.split(r'_', umi_name)[0])
                info = line[7]
                dp4 = re.split(r'DP4=|;MQ=', info)[1]
                dp4 = re.split(r',', dp4)
                var_reads = int(dp4[2]) + int(dp4[3])
                ref_reads = int(dp4[0]) + int(dp4[1])
                total_coverage = var_reads + ref_reads
                var_percent = round(var_reads/total_coverage, 2)
                dp4_percent = round(total_coverage/read_count, 2)
                new_info = info + ';VARP=' + str(var_percent) + ';DP4P=' + str(dp4_percent)
                revised_line = tuple(line[0:7] + [new_info] + line[8:])
                outfile1.writelines('\t'.join(revised_line))

                # if var_reads > 1 and var_percent > 0.2 and not info.startswith('INDEL'):
                # var_reads > 1 is for excluding very rare variants, var_percent < 0.9 is for excluding species specific variants
                if var_reads > 1 and var_percent < 0.9 and not info.startswith('INDEL'):
                    if umi_name == process_umi:
                        record_count += 1
                        var_percent_list = var_percent_list + [var_percent]
                        dp4_percent_list = dp4_percent_list + [dp4_percent]
                        total_var_reads += var_reads
                        total_ref_reads += ref_reads

                    elif process_umi is '':
                        process_umi = umi_name
                        record_count += 1
                        var_percent_list = [var_percent]
                        dp4_percent_list = [dp4_percent]
                        total_var_reads = var_reads
                        total_ref_reads = ref_reads

                    else:
                        average_var_pcent = round(sum(var_percent_list) / len(var_percent_list), 2)
                        average_dp4_pcent = round(sum(dp4_percent_list) / len(dp4_percent_list), 2)
                        total_read_count = total_var_reads + total_ref_reads
                        var_read_pcet = round(total_var_reads / total_read_count, 2)
                        reads_combine = str(total_ref_reads) + ',' + str(total_var_reads)

                        if len(var_percent_list) > 1:
                            var_stdev = round(statistics.stdev(var_percent_list), 2)
                        else:
                            var_stdev = 0.0

                        result = [process_umi, record_count, var_read_pcet, total_read_count, reads_combine,
                                  average_var_pcent, str(var_stdev), average_dp4_pcent] + \
                                 [','.join(str(e) for e in var_percent_list)] + [','.join(str(e) for e in dp4_percent_list)]
                        outfile2.writelines('\t'.join(str(e) for e in result) + '\n')
                        process_umi = umi_name
                        record_count = 1
                        var_percent_list = [var_percent]
                        dp4_percent_list = [dp4_percent]
                        total_var_reads = var_reads
                        total_ref_reads = ref_reads


if __name__ == '__main__':
    save_path = './'
    outvcf = './all_snp_from_perfect_umi.pcent.vcf'
    add_pcent(save_path, outvcf)

    cmd = """cat umi_group.flt.summary.txt | awk '{if($3>0.35 && $3<0.65)print}' | awk '{if($7 < 0.4) print}' > wrong.group.summary.txt
             cat wrong.group.summary.txt | cut -f1 | sort -k1 > wrong.group.lst
             ls ./perfect_umi | sort | uniq > all.group.lst
             comm -23 all.group.lst wrong.group.lst > pass.group.lst

             split -l 100 pass.group.lst pattern-file.split.
             cat all_snp_from_perfect_umi.pcent.vcf | grep "^#" > all_snp_from_perfect_umi.pcent.rem.vcf
             for CHUNK in pattern-file.split.* ; do
                cat all_snp_from_perfect_umi.pcent.vcf | grep -v "^#" | grep -f "$CHUNK"  >> all_snp_from_perfect_umi.pcent.rem.vcf
             done
             rm pattern-file.split.* all.group.lst

             cat all_snp_from_perfect_umi.pcent.rem.vcf | \
                bcftools filter -s LowQual -e '%QUAL<10 || INFO/DP<3 || INFO/IMF<0.5 || INFO/VARP<0.67' \
                > all_snp_from_perfect_umi.pcent.rem.flt.vcf

             cat all_snp_from_perfect_umi.pcent.rem.flt.vcf | grep -v LowQual \
                > pass_snp_from_perfect_umi.flt.vcf"""

    subprocess.run(cmd, shell=True, check=True, stdout=subprocess.PIPE)
