#!/bin/bash
# Chongwei 20181011
# Email: chongwei.bi@kaust.edu.sa


function usage()
{
cat << endaaa
-------------------------------------------------------------------------------------------------------------------------------------
need:   cutadapt, version=2.0

sh group_umi.sh -n {name} -l {left_flank} -u {umi_NN} -r {right_flank} -q {dir_fastq} -s {dir_save} -e {error_rate} -t {thread}
                -5 :<5 end UMI>
                -3 :<3 end UMI>

~/bic/jobs/group_umi.sh -n mt_bc2 -l CATCTTACGATTACGCCAACCACTG -u NNNTGNNN -r CTCCCGAATCAACCCTGACCC -q ./file.fastq -s ./ -e 0.1 -g

                                                                                                By Chongwei
-------------------------------------------------------------------------------------------------------------------------------------
endaaa
}


function umi_analysis() # $l $u $r ${umi_end} ${g/a} ${thread}
{
    mkdir ${s}/${n}/${umi_end}end_UMIs
    s="${s}/${n}/${umi_end}end_UMIs"   # Direction to save

    l=$(echo ${l} | tr [a-z] [A-Z])
    lrev=$(echo ${l} | tr [a-z] [A-Z] | tr [ATGC] [TACG] | rev)
    u=$(echo ${u} | tr [a-z] [A-Z])
    urev=$(echo ${u} | tr [a-z] [A-Z] | tr [ATGC] [TACG] | rev)
    r=$(echo ${r} | tr [a-z] [A-Z])
    rrev=$(echo ${r} | tr [a-z] [A-Z] | tr [ATGC] [TACG] | rev)
    umi="${l}${u}${r}"
    umi_rev=$(echo ${umi} | tr [a-z] [A-Z] | tr [ATGC] [TACG] | rev)

    if [[ "${umi_end}" -eq 3 ]];then
        umi=${umi_rev}
        l=${rrev}
        r=${lrev}
    fi

    length="$(echo "${u}" | awk '{print length()}')"
    L3="$[${length}-1]"
    R3="$[${length}+1]"

    echo "---${umi_end} end UMI start---"
    echo "${umi_end} end UMIs sequence is ${umi}"
    echo "Length of UMIs is ${length}, group reads based on UMIs length from ${L3} to ${R3}"

    #source activate bioconda
    echo "--- ${n}.extract.${umi}.summary" > ${s}/${n}.00.cut_adapter.summary

    # Get UMIs adapter sequence
    cutadapt -e ${error} --overlap 3 -j ${thread} -${para} Whole_UMIs_Adapter=${umi} --action=lowercase --discard-untrimmed -o ${s}/${n}.01.reads_with_adapter.fastq \
        ${input} --info-file=${s}/${n}.01.extract_reads_with_adapter.detail >> ${s}/${n}.00.cut_adapter.summary 2>&1

    if [[ -s "${s}/${n}.01.extract_reads_with_adapter.detail" ]];then
        :
    else
        echo "====== can't generate extract_reads_with_adapter.detail file, maybe something is wrong with cutadapt? \
             Please check ${s}/${n}.00.cut_adapter.summary ======" >&2
        exit 2
    fi

    cat ${s}/${n}.01.extract_reads_with_adapter.detail | awk -F "\t" '{if(NF==11 && $8=="Whole_UMIs_Adapter") {print ">"$1"\n"$6}}' > \
        ${s}/${n}.02.extracted.UMI_adapter.fasta

    # use non-internal adapter {-g adapter} and {-a adapterX} to cut left and cut right
    echo "--- ${n}.cut_left_flank.${l}.summary" >> ${s}/${n}.00.cut_adapter.summary
    cutadapt -e ${error} -j ${thread} -g left_flank=${l} --action=lowercase --discard-untrimmed -o ${s}/${n}.03.adapter_with_left_flank.fasta \
        ${s}/${n}.02.extracted.UMI_adapter.fasta --info-file=${s}/${n}.03.cut_left_flank.detail >> ${s}/${n}.00.cut_adapter.summary 2>&1

    cat ${s}/${n}.03.cut_left_flank.detail | awk -F "\t" '{if(NF==11 && $8=="left_flank") {print ">"$1"\n"$7}}' > ${s}/${n}.04.adapter_without_left_flank.fasta

    echo "--- ${n}.cut_right_flank.${r}X.summary" >> ${s}/${n}.00.cut_adapter.summary
    cutadapt -e ${error} -j ${thread} -a right_flank=${r}X --action=lowercase --discard-untrimmed -o ${s}/${n}.05.adapter_with_right_flank.fasta \
        ${s}/${n}.04.adapter_without_left_flank.fasta --info-file=${s}/${n}.05.cut_right_flank.detail >> /dev/null 2>&1

    cutadapt -e ${error} -j ${thread} -a right_flank=${r}X --action=trim --discard-untrimmed -o ${s}/${n}.05.adapter_without_right_flank.fasta \
        ${s}/${n}.04.adapter_without_left_flank.fasta >> ${s}/${n}.00.cut_adapter.summary 2>&1

    cat ${s}/${n}.05.adapter_without_right_flank.fasta | paste -d ","  - - | awk -F ">" '{print $2}' | awk -F "," '{print length($2)","$2","$1}' \
        | awk -F "," '{if($1=='${length}') {print $2","$3}}' | awk -F " " '{print $1}' > ${s}/${n}.06.1umi.2read_name.csv

    cat ${s}/${n}.06.1umi.2read_name.csv | \
    awk -F "," '
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
    }' | awk '{print NF-1" "$0}' > ${s}/${n}.07.1read_count.2umi.3read_name.lst

    echo "---${umi_end} end Done!---"
    date
}


if [[ $# -lt 1 ]]; then
    usage
    exit
fi

#==========================#
#= get required arguments =#
#==========================#

n=""           # name
s=""           # save path
l=""           # UMIs left flank
u=""           # UMIs: NNNNNNNN
r=""           # UMIs right flank
input=""       # Direction to fastq
error=""
umi_end=""
para=""        # parameter for umi adapter extraction
thread="0"

while getopts ":n:q:l:u:r:s:e:t:35" option;
do
    case ${option} in
        n)
            n=$OPTARG
        ;;
        q)
            input=$OPTARG
        ;;
        l)
            l=$OPTARG
        ;;
        u)
            u=$OPTARG
        ;;
        r)
            r=$OPTARG
        ;;
        s)
            s=$OPTARG
        ;;
        e)
            error=$OPTARG
        ;;
        t)
            thread=$OPTARG
        ;;
        3)
            umi_end="3"
            para="a"
        ;;
        5)
            umi_end="5"
            para="g"
        ;;
        ?)
            echo "Invalid option: -$OPTARG" >&2
            exit 1
        ;;
        :)
            echo "Option -$OPTARG requires an argument." >&2
            exit 1
        ;;
    esac
done

if [[ ! -d "${s}/${n}" ]];then
    mkdir ${s}/${n} 2>/dev/null
fi

################### process UMI analysis ###################
umi_analysis

rm ${s}/${n}.01.reads_with_adapter.fastq
rm ${s}/${n}.01.extract_reads_with_adapter.detail