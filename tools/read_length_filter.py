import subprocess


def filter_length(fastq, minlength, maxlength, save_path, file_name):

    if maxlength != 0:
        out_fastq = save_path + file_name + '.' + str(minlength) + '-' + str(maxlength) + '.fastq'

        awk = r'''awk 'BEGIN {FS = "\t" ; OFS = "\n"} 
                    {header = $0 ; getline seq ; getline qheader ; getline qseq 
                    if (length(seq) >= %d && length(seq) <= %d) 
                    {print header, seq, qheader, qseq}}' ''' % (minlength, maxlength)

        cmd = """ cat {fastq} | seqtk seq -l0 | {awk} > {out_fastq}
                """.format(fastq=fastq,
                           awk=awk,
                           out_fastq=out_fastq)
        subprocess.run(cmd, shell=True, check=True, stdout=subprocess.PIPE)

    else:
        out_fastq = save_path + file_name + '.' + str(minlength) + '-x.fastq'

        awk = r'''awk 'BEGIN {FS = "\t" ; OFS = "\n"} 
                            {header = $0 ; getline seq ; getline qheader ; getline qseq 
                            if (length(seq) >= %d) 
                            {print header, seq, qheader, qseq}}' ''' % minlength
        cmd = """ cat {fastq} | seqtk seq -l0 | {awk} > {out_fastq}
                        """.format(fastq=fastq,
                                   awk=awk,
                                   out_fastq=out_fastq)
        subprocess.run(cmd, shell=True, check=True, stdout=subprocess.PIPE)

    return out_fastq
