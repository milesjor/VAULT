import subprocess
import sys
import logging
import os
import shlex

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


def finalize_snv_vcf(vcf_file, variant_id=None):
    """Fix bcftools INFO/MQ type and optionally set variant IDs."""
    if not os.path.exists(vcf_file):
        return

    tmp_file = vcf_file + ".tmp"

    with open(vcf_file, "r") as infile, open(tmp_file, "w") as outfile:
        for line in infile:
            if (
                line.startswith("##INFO=<ID=MQ,")
                and "Type=Integer" in line
            ):
                line = line.replace("Type=Integer", "Type=Float", 1)

            elif variant_id and not line.startswith("#"):
                fields = line.rstrip("\n").split("\t")
                if len(fields) >= 3:
                    fields[2] = variant_id
                    line = "\t".join(fields) + "\n"

            outfile.write(line)

    os.replace(tmp_file, vcf_file)


def run_logged_command(cmd, save):
    try:
        result = subprocess.run(
            cmd,
            shell=True,
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            universal_newlines=True,
        )
    except subprocess.CalledProcessError as e:
        if e.stdout:
            logging.error("stdout from " + save + " :\n" + e.stdout)
        if e.stderr:
            logging.error("stderr from " + save + " :\n" + e.stderr)
        logging.error(
            "Command failed in %s with exit code %s",
            save,
            e.returncode,
        )
        raise

    if result.stdout:
        logging.info("stdout from " + save + " :\n" + result.stdout)
    if result.stderr:
        logging.info("stderr from " + save + " :\n" + result.stderr)

    return result


def snv_calling(save, thread, ax, refer, fastq, name):
    count_number = """awk 'BEGIN{n=0}{if($3>=3)n+=1}END{print n}' """
    awk_identity = """awk 'BEGIN{a=0;b=0}{a+=$10/$11;b++}END{print a/b}' """
    raw_vcf = "{save}/vcf/{name}.snv.raw.vcf".format(
        save=save,
        name=name,
    )
    final_vcf = "{save}/vcf/{name}.snv.vcf".format(
        save=save,
        name=name,
    )

    align_cmd = """# alignment
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
             """.format(save=save,
                        thread=thread,
                        ax=ax,
                        refer=refer,
                        fastq=fastq,
                        name=name,
                        awk_identity=awk_identity)

    call_cmd = """# SNV calling
             bcftools mpileup --threads {thread} -Q 10 -Ou -A -f {refer} {save}/{name}.sv.sorted.bam 2>>{save}/{name}_alignment_summary.log | \
             bcftools call --threads {thread} -Ov -mv 2>>{save}/{name}_alignment_summary.log > {raw_vcf}
             """.format(save=save,
                        thread=thread,
                        refer=refer,
                        name=name,
                        raw_vcf=raw_vcf)

    norm_cmd = """# SNV normalization
             bcftools norm --threads {thread} -Ov -f {refer} {raw_vcf} 2>>{save}/{name}_alignment_summary.log > {final_vcf}
             samtools depth -d 0 {save}/{name}.sv.sorted.bam > {save}/{name}.coverage.txt
             cat {save}/{name}.coverage.txt | {count_number} | xargs echo {name} > {save}/{name}.coverage.3plus.txt
             rm -f {raw_vcf}
             """.format(save=save,
                        thread=thread,
                        refer=refer,
                        name=name,
                        raw_vcf=raw_vcf,
                        final_vcf=final_vcf,
                        count_number=count_number)

    for cmd in (align_cmd, call_cmd):
        run_logged_command(cmd, save)

    finalize_snv_vcf(raw_vcf)

    run_logged_command(norm_cmd, save)
    finalize_snv_vcf(final_vcf, name)


def sv_calling(
    save,
    thread,
    refer,
    fastq,
    name,
    map_mode='-x ont',
    sv_caller='sniffles1',
    sniffles_bin='sniffles',
    min_support=3,
):
    bam = shlex.quote(os.path.join(save, name + ".sv.sorted.bam"))
    vcf = shlex.quote(os.path.join(save, "vcf", name + ".sv.vcf"))
    log = shlex.quote(os.path.join(save, name + "_alignment_summary.log"))
    sniffles_cmd = shlex.quote(str(sniffles_bin or "sniffles"))

    if sv_caller == "sniffles2":
        call_cmd = (
            f"{sniffles_cmd} --input {bam} "
            f"--vcf {vcf} "
            f"--threads {thread} "
            f"--minsupport {min_support} >> {log} 2>&1"
        )
    else:
        call_cmd = (
            f"{sniffles_cmd} -t {thread} "
            f"-s {min_support} "
            f"-r 100 -n -1 -q 20 -d 1000 -l 30 "
            f"--skip_parameter_estimation "
            f"-m {bam} "
            f"-v {vcf} >> {log} 2>&1"
        )

    cmd = """# SV calling
             # ngmlr -t {thread} -r {refer} -q {fastq} -o {save}/{name}.sv.sam  {map_mode} 2>>{save}/{name}_alignment_summary.log
             # echo "====== Alignment summary ======" >> {save}/{name}_alignment_summary.log
             # grep -v "^@" {save}/{name}.sv.sam | cut -f 3 | sort | uniq -c | sort -n >> {save}/{name}_alignment_summary.log
             # echo "===============================" >> {save}/{name}_alignment_summary.log
             # samtools sort -@ {thread} {save}/{name}.sv.sam -o {save}/{name}.sv.sorted.bam
             # samtools index -@ {thread} {save}/{name}.sv.sorted.bam
             # sniffles -t 5 -s 4 -r 2000 -q 20 -d 500 -l 30 -m {save}/{name}.sv.sorted.bam -v {save}/vcf/{name}.sv.vcf         
             # bcftools view -i '(SVTYPE = "DUP" || SVTYPE = "INS" || SVTYPE = "DEL") && ABS(SVLEN) > 49' {save}/vcf/{name}.sv.vcf \
             # > {save}/vcf/{name}.sv.flt.vcf
             
             {call_cmd}
             """.format(save=save,
                        thread=thread,
                        refer=refer,
                        fastq=fastq,
                        name=name,
                        map_mode=map_mode,
                        call_cmd=call_cmd)

    run_logged_command(cmd, save)


if __name__ == '__main__':
    sys.stderr.write('ERROR! This is a module. See commands with `vault -h`')
    sys.exit(1)
