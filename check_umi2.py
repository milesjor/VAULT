import subprocess
import sys


def run_umi_analysis(name, fastq, thread, end, para, save, error, left, umi, right):
    awk1 = """awk '{print length()}' """
    awk2 = """awk -F "\\t" '{if(NF==11 && $8=="Whole_UMIs_Adapter") {print ">"$1"\\n"$6}}' """
    awk3 = """awk -F "\\t" '{if(NF==11 && $8=="left_flank") {print ">"$1"\\n"$7}}' """
    awk4 = """awk -F ">" '{print $2}' | awk -F "," '{print length($2)","$2","$1}' \
              | awk -F "," '{if($1=='$length') {print $2","$3}}' | awk -F " " '{print $1}' """
    awk5 = """awk -F "," '
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
    
    cmd = """
            if [[ ! -d "{save}/{n}" ]]; then mkdir {save}/{n} 2>/dev/null; fi
            
            mkdir {save}/{n}/{umi_end}end_UMIs
            s="{save}/{n}/{umi_end}end_UMIs"   # Direction to save
            left=$(echo {left} | tr [a-z] [A-Z])
            lrev=$(echo {left} | tr [a-z] [A-Z] | tr [ATGCBDHKMRSVWY] [TACGVHDMKYSBWR] | rev)
            u=$(echo {u} | tr [a-z] [A-Z])
            urev=$(echo {u} | tr [a-z] [A-Z] | tr [ATGCBDHKMRSVWY] [TACGVHDMKYSBWR] | rev)
            right=$(echo {right} | tr [a-z] [A-Z])
            rrev=$(echo {right} | tr [a-z] [A-Z] | tr [ATGCBDHKMRSVWY] [TACGVHDMKYSBWR] | rev)
            umi="$left$u$right"
            umi_rev=$(echo $umi | tr [a-z] [A-Z] | tr [ATGCBDHKMRSVWY] [TACGVHDMKYSBWR] | rev)

            if [[ "{umi_end}" -eq 3 ]];then
                umi=$umi_rev
                left=$rrev
                right=$lrev
            fi

            length="$(echo "{u}" | {awk1})"
            L3="$[$length-1]"
            R3="$[$length+1]"

            # echo "---{umi_end} end UMI start---"
            # echo "{umi_end} end UMIs sequence is $umi"
            # echo "Length of UMIs is $length, group reads based on UMIs length from $L3 to $R3"

            echo "--- {n}.extract.$umi.summary" > $s/{n}.00.cut_adapter.summary

            # Get UMIs adapter sequence
            cutadapt -e {error} --overlap 3 -j {thread} -{para} Whole_UMIs_Adapter=$umi --action=lowercase \
                --discard-untrimmed -o $s/{n}.01.reads_with_adapter.fastq \
                {input} --info-file=$s/{n}.01.extract_reads_with_adapter.detail >> $s/{n}.00.cut_adapter.summary 2>&1

            if [[ -s "$s/{n}.01.extract_reads_with_adapter.detail" ]];then
                :
            else
                echo "====== can't generate extract_reads_with_adapter.detail file, maybe something is wrong with cutadapt? \
                    \n       Please check $s/{n}.00.cut_adapter.summary ======" >&2
                exit 2
            fi

            cat $s/{n}.01.extract_reads_with_adapter.detail | {awk2} > $s/{n}.02.extracted.UMI_adapter.fasta

            # use non-internal adapter '-g adapter' and '-a adapterX' to cut left and cut right
            echo "--- {n}.cut_left_flank.$left.summary" >> $s/{n}.00.cut_adapter.summary
            cutadapt -e {error} -j {thread} -g left_flank=$left --action=lowercase --discard-untrimmed \
                -o $s/{n}.03.adapter_with_left_flank.fasta $s/{n}.02.extracted.UMI_adapter.fasta \
                --info-file=$s/{n}.03.cut_left_flank.detail >> $s/{n}.00.cut_adapter.summary 2>&1

            cat $s/{n}.03.cut_left_flank.detail | {awk3} > $s/{n}.04.adapter_without_left_flank.fasta

            echo "--- {n}.cut_right_flank.$rightX.summary" >> $s/{n}.00.cut_adapter.summary
            cutadapt -e {error} -j {thread} -a right_flank=$right""X --action=lowercase --discard-untrimmed \
                -o $s/{n}.05.adapter_with_right_flank.fasta $s/{n}.04.adapter_without_left_flank.fasta \
                --info-file=$s/{n}.05.cut_right_flank.detail >> /dev/null 2>&1

            cutadapt -e {error} -j {thread} -a right_flank=$right""X --action=trim --discard-untrimmed \
                -o $s/{n}.05.adapter_without_right_flank.fasta $s/{n}.04.adapter_without_left_flank.fasta \
                >> $s/{n}.00.cut_adapter.summary 2>&1

            cat $s/{n}.05.adapter_without_right_flank.fasta | paste -d ","  - - | {awk4} > $s/{n}.06.1umi.2read_name.csv

            cat $s/{n}.06.1umi.2read_name.csv | {awk5} > $s/{n}.07.1read_count.2umi.3read_name.lst

            # echo "---{umi_end} end Done!---"
            # date
            
            rm $s/{n}.01.reads_with_adapter.fastq
            rm $s/{n}.01.extract_reads_with_adapter.detail
             """.format(save=save,
                        n=name,
                        umi_end=end,
                        left=left,
                        u=umi,
                        right=right,
                        error=error,
                        thread=thread,
                        para=para,
                        input=fastq,
                        awk1=awk1,
                        awk2=awk2,
                        awk3=awk3,
                        awk4=awk4,
                        awk5=awk5)

    a = subprocess.run(cmd, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if a.stdout != b'':
        print("stdout from " + save + " :\n" + a.stdout.decode('utf-8'))
    if a.stderr != b'':
        print("stderr from " + save + " :\n" + a.stderr.decode('utf-8'))


if __name__ == '__main__':
    sys.stderr.write('ERROR! This is a module. See commands with `vault -h`')
    sys.exit(1)
