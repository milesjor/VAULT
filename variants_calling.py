import subprocess
import sys
import logging

"""
bcftools mpileup options
-A, --count-orphans
Do not skip anomalous read pairs in variant calling.
-q, -min-MQ INT
Minimum mapping quality for an alignment to be used [0]
-Q, --min-BQ INT
Minimum base quality for a base to be considered [13]
-x, --ignore-overlaps
Disable read-pair overlap detection.
"""


def snp_calling(save, thread, ax, refer, fastq, name):
    count_number = """awk 'BEGIN{n=0}{if($3>=3)n+=1}END{print n}' """
    awk_identity = """awk 'BEGIN{a=0;b=0}{a+=$10/$11;b++}END{print a/b}' """
    cmd = """# alignment
             minimap2 -t {thread} -a -x {ax} -Y --MD {refer} {fastq} > {save}/{name}.sv.sam 2>>{save}/{name}_alignment_summary.log
             # ngmlr -t {thread} -i 0.8 -x ont -r {refer} -q {fastq} -o {save}/{name}.sv.sam 2>>{save}/{name}_alignment_summary.log
             # minimap2 -t {thread} -a -x asm10 -Y --MD {refer} {fastq} > {save}/{name}.sv.sam 2>>{save}/{name}_alignment_summary.log
             
             # alignment identity
             # minimap2 -t {thread} -x asm10 -Y --MD {refer} {fastq} > {save}/{name}.sv.paf 2>>{save}/{name}_alignment_summary.log
             # cat {save}/{name}.sv.paf | {awk_identity} > {save}/{name}.sv.align_identity.txt
             
             echo "====== Alignment summary ======" >> {save}/{name}_alignment_summary.log
             grep -v "^@" {save}/{name}.sv.sam | cut -f 3 | sort | uniq -c | sort -n >> {save}/{name}_alignment_summary.log
             echo "===============================" >> {save}/{name}_alignment_summary.log
             samtools sort -@ {thread} {save}/{name}.sv.sam -o {save}/{name}.sv.sorted.bam >> {save}/{name}_alignment_summary.log
             samtools index -@ {thread} {save}/{name}.sv.sorted.bam
             # SNP calling
             bcftools mpileup --threads {thread} -Q 10 -Ou -A -f {refer} {save}/{name}.sv.sorted.bam 2>>{save}/{name}_alignment_summary.log | \
             bcftools call --threads {thread} -Ou -mv 2>>{save}/{name}_alignment_summary.log | \
             bcftools norm --threads {thread} -Ou -f {refer} 2>>{save}/{name}_alignment_summary.log | \
             bcftools annotate --set-id ${name} -Ov 1> {save}/vcf/{name}.snp.vcf
             # bcftools filter -s LowQual -e '%QUAL<10 || INFO/DP<5 || INFO/IMF<0.5' > {save}/vcf/{name}.snp.vcf
             samtools depth -d 0 {save}/{name}.sv.sorted.bam > {save}/{name}.coverage.txt
             cat {save}/{name}.coverage.txt | {count_number} | xargs echo {name} > {save}/{name}.coverage.3plus.txt
             """.format(save=save,
                        thread=thread,
                        ax=ax,
                        refer=refer,
                        fastq=fastq,
                        name=name,
                        count_number=count_number,
                        awk_identity=awk_identity)

    a = subprocess.run(cmd, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if a.stdout != b'':
        logging.info("stdout from " + save + " :\n" + a.stdout.decode('utf-8'))
    if a.stderr != b'':
        logging.info("stderr from " + save + " :\n" + a.stderr.decode('utf-8'))


def sv_calling(save, thread, refer, fastq, name, map_mode='-x ont'):
    cmd = """# SV calling
             # ngmlr -t {thread} -r {refer} -q {fastq} -o {save}/{name}.sv.sam  {map_mode} 2>>{save}/{name}_alignment_summary.log
             # echo "====== Alignment summary ======" >> {save}/{name}_alignment_summary.log
             # grep -v "^@" {save}/{name}.sv.sam | cut -f 3 | sort | uniq -c | sort -n >> {save}/{name}_alignment_summary.log
             # echo "===============================" >> {save}/{name}_alignment_summary.log
             # samtools sort -@ {thread} {save}/{name}.sv.sam -o {save}/{name}.sv.sorted.bam
             # samtools index -@ {thread} {save}/{name}.sv.sorted.bam
             # sniffles -t 5 -s 4 -r 2000 -q 20 -d 500 -l 30 -m {save}/{name}.sv.sorted.bam -v {save}/vcf/{name}.sv.vcf
             
             sniffles -t 3 -s 3 -r 1500 -q 20 -d 500 -l 30 -m {save}/{name}.sv.sorted.bam -v {save}/vcf/{name}.sv.vcf >> {save}/{name}_alignment_summary.log
             
             # bcftools view -i '(SVTYPE = "DUP" || SVTYPE = "INS" || SVTYPE = "DEL") && ABS(SVLEN) > 49' {save}/vcf/{name}.sv.vcf \
             # > {save}/vcf/{name}.sv.flt.vcf
             """.format(save=save,
                        thread=thread,
                        refer=refer,
                        fastq=fastq,
                        name=name,
                        map_mode=map_mode)

    a = subprocess.run(cmd, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if a.stdout != b'':
        logging.info("stdout from " + save + " :\n" + a.stdout.decode('utf-8'))
    if a.stderr != b'':
        logging.info("stderr from " + save + " :\n" + a.stderr.decode('utf-8'))


if __name__ == '__main__':
    sys.stderr.write('ERROR! This is a module. See commands with `vault -h`')
    sys.exit(1)
