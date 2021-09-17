# VAULT
 Variant Analysis with UMI for Long-read Technology (VAULT)

VAULT is a tool for analyzing UMI-labeled reads, works for both error-prone long reads and accurate single-end/paired-end short reads.
<p align="center">
<img src="https://github.com/milesjor/VAULT/blob/master/example/pic/compare.png" width="600"/>
</p>

More detail: [Long-read Individual-molecule Sequencing Reveals CRISPR-induced Genetic Heterogeneity in Human ESCs](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-02143-8)

## INSTALLATION
### Prerequisites
Anaconda (https://www.anaconda.com/distribution/) or Miniconda (https://conda.io/miniconda.html).
For example, users may download and install the following Anaconda3 package:
```
wget https://repo.anaconda.com/archive/Anaconda3-2019.10-Linux-x86_64.sh
bash Anaconda3-2019.10-Linux-x86_64.sh
```

### Download the VAULT package
```
git clone https://github.com/milesjor/VAULT.git
cd ./VAULT/
```

### Install all required modules
```
conda env create --name vault --file ./tools/vault_env.yaml
conda activate vault
```

Optional: add VAULT to your PATH
```
printf "\n# VAULT added \nexport PATH="\$PATH:$(pwd)"\n\n" >> ~/.bashrc
source ~/.bashrc
```

## USAGE
<p align="center">
<img src="https://github.com/milesjor/VAULT/blob/master/example/pic/vault_pipeline.png" width="600"/>
</p>

### The parameters of VAULT
Users can use the following command to print out the full usage:
```
./vault -h
```
If VAULT is added in your PATH:
```
vault -h
```

A list of available parameters:
```
usage: vault [-h] [-v] [-u UMI_ADAPTER] [-s SAVE_PATH] [-r REFER] [-q FASTQ]
             [-e ERROR] [-t THREAD] [-T THRESHOLD] [-b BASH_THREAD]
             [-F ALLELE_FREQ] [-f SV_FREQ] [-p PE_FASTQ]
             [-a {sr,map-ont,map-pb}] [--minlength MINLENGTH]
             [--maxlength MAXLENGTH] [--unmapped_reads] [--group_filter]
             {summarize,consensus,position,circos,filter,vaf} ...

This is for analyzing UMI labeled reads in IDMseq and iMiGseq.

optional arguments:
  -h, --help            show this help message and exit
  -v, --version         show the current version

subcommands:
  valid subcommands

  {consensus,position,circos,filter,vaf}
                        additional help
    consensus           get consensus sequence from VAULT result
    position            correct the position and chr_name in vcf file
    circos              prepare data for circos, used in IMTseq
    filter              filter out low-confidence UMI groups
    vaf                 calculate variant allele frequency (VAF) based on UMI
                        group number

Required options:
  -u UMI_ADAPTER, --umi_adapter UMI_ADAPTER
                        UMI sequence, automatically detect NNNATGCNNN as UMI
  -s SAVE_PATH, --save_path SAVE_PATH
                        path/to/save/
  -r REFER, --refer REFER
                        path/to/ref.fa
  -q FASTQ, --fastq FASTQ
                        path/to/reads.fastq, or fastq.gz

Optional options:
  -e ERROR, --error ERROR               # related to raw read error rate
                        error tolerate rate in umi analysis (depend on read
                        error rate) [0.11]
  -t THREAD, --thread THREAD            # thread for processing UMI groups
                        parallel process number [5]
  -T THRESHOLD, --threshold THRESHOLD   # ignore UMI group with read number <= [5]
                        threshold of read number for snp analysis [5]
  -b BASH_THREAD, --bash_thread BASH_THREAD    # thread in variant calling
                        thread for running bash cmd [1]
  -F ALLELE_FREQ, --allele_freq ALLELE_FREQ
                        filter SNVs and SVs by allele frequency [0.67]
  -p PE_FASTQ, --pe_fastq PE_FASTQ
                        read2.fastq for illumina pair-end sequencing
  -a {sr,map-ont,map-pb}, --align_mode {sr,map-ont,map-pb}
                        parameter in alignment, minimap2 -ax [sr|map-ont|map-
                        pb]
  --minlength MINLENGTH
                        filter fastq file to remove reads with length less
                        than [int]
  --maxlength MAXLENGTH
                        filter fastq file to remove reads with length more
                        than [int]
  --unmapped_reads      extract mapped reads before UMI analysis  
  --group_filter        filter out low-confidence UMI groups      
```

## EXAMPLE

### Command
```
./vault -u CATCTTACGATTACGCCAACCACTGCGGNNNNNTGNNNNNGACACATTCTCCCAGGCCCTACTT \
        -q ./example/nanopore_reads.fastq.gz \
        -s ./example/result \
        -r ./example/reference.fa \
        -e 0.11 \
        -a map-ont \
        --minlength 300 \
        --maxlength 20000 \
        --unmapped_reads \
        --group_filter \
        -t 4 \
        -b 1
```

### Results
```
./result/
├── 20210622_08.46.11_vault.log             # VAULT log
├── nanopore_reads.300-20000_alignment_summary.log   # raw reads alignment summary
├── nanopore_reads.300-20000.bam
├── nanopore_reads.300-20000.bam.bai
├── nanopore_reads.300-20000.fastq
├── nanopore_reads.300-20000.mapped.fastq   # length filtered and alignment mapped reads used in the VAULT analysis
├── nanopore_reads.300-20000.mapped.lst
├── nanopore_reads.300-20000.sam
├── grouped_reads                           # fastq reads for every UMI groups
│   └── perfect_umi
├── snp  # folder for variant analysis
│   ├── all_snp_from_perfect_umi.vcf        # raw variant calling result for small variants
│   ├── all_snp_from_perfect_umi.pcent.vcf  # add variant allele frequency (supported reads percentage) in the info field of vcf file
│   ├── all_snp_from_perfect_umi.pcent.rem.vcf              # remove wrong UMI group
│   ├── all_snp_from_perfect_umi.pcent.rem.flt.vcf          # filter by depth, quality, allele frequency
│   ├── all_sv_from_perfect_umi.vcf         # raw variant calling result for large variants (>= 30bp)
│   ├── all_sv_from_perfect_umi.filtered.0.67.vcf           # FINAL SV result (remove wrong UMI group and filter SVs by allele frequency [0.5])
│   ├── all_sv_from_perfect_umi.filtered.0.67.sorted.vcf    # sort by position
│   ├── coverage.3plus.txt                  # The region with coverage >= 3 in each UMI groups, can be used to filter out UMI groups
│   ├── pass.group.lst                      # UMI groups that pass group_filter
│   ├── pass_snp_from_perfect_umi.flt.vcf   # FINAL SNV and InDel result (remove wrong UMI group and filter by depth, quality, allele frequency)
│   ├── umi_group.flt.summary.txt           # intermediate file in [--group_filter]
│   ├── wrong.group.lst                     # UMI groups that fail in [--group_filter]
│   ├── wrong.group.summary.txt             # intermediate file in [--group_filter] for wrong UMI groups
│   └── perfect_umi                         # individual UMI analysis result for 5' and 3' end of reads
└── umi_analysis
    ├── 3end_UMIs
    └── 5end_UMIs
```
In the "snp" folder, the files with "all" prefix mean that they are from all UMI groups in ./snp/perfect_umi folder. The files with "pass" prefix mean they are from "--group_filter" passed UMI groups, which is shown in ./snp/pass.group.lst
### Individual UMI group folder
```
./result/snp/perfect_umi/
├── 14_ATCGATGATTTT_AAAATCATCGAT    # 14 reads in this group, 5' UMI is ATCGATGATTTT, 3' is AAAATCATCGAT
├── 33_GACATTGTCTGG_CCAGACAATGTC
├── 35_5end_AACAGTGCTGCT            # 35 reads in this group, all reads from 5' UMI
├── 5_3end_AAAAACATGGCA             # 5 reads in this group, all reads from 3' UMI
├── 7_5end_ATTCTTGGTGTC
├── 7_CTATGTGAAGAA_TTCTTCACATAG
├── 8_3end_ACAAGCAAAAAA
├── 8_AGTTGTGCCATA_TATGGCACAACT
├── 8_CCGCGTGAGATG_CATCTCACGCGG
├── 8_CGTTGTGTTACT_AGTAACACAACG
├── 8_CTATTTGTCACT_AGTGACAAATAG
├── 8_GGGTTTGGTTTG_CAAACCAAACCC
├── 8_GTGGGTGACGGG_CCCGTCACCCAC
├── 8_GTGTTTGTTAGA_TCTAACAAACAC
├── 8_TCAATTGCAGAA_TTCTGCAATTGA
├── 8_TTACTTGATTTT_AAAATCAAGTAA
├── 8_TTGGATGGAAGT_ACTTCCATCCAA
├── 9_AAAGATGCGCGT_ACGCGCATCTTT
├── 9_AGAAATGATAGC_GCTATCATTTCT
├── 9_ATCGATGGTGCG_CGCACCATCGAT
├── 9_ATGTTTGCCAAT_ATTGGCAAACAT
├── 9_CTAACTGCTTAT_ATAAGCAGTTAG
├── 9_GAGAATGAGTAC_GTACTCATTCTC
└── 9_GTTTATGTACAT_ATGTACATAAAC
```

## TIPS

### 5' end UMI and 3' end UMI
<p align="center">
<img src="https://github.com/milesjor/VAULT/blob/master/example/pic/connect_two_end_umi.png" width="600"/>
</p>

### Tips on increasing the analyzable reads (reads with UMIs)
To potentially improve the usable reads, VAULT offered a parameter [-e] to define the error tolerant threshold. The default value [0.11] is set based on the average sequencing error rate of Nanopore. When increasing this value, we will have more reads with UMIs. But this may compromise the confidence of UMI groups.

The other parameter can be adjusted is [-u], it asks for a DNA sequence with UMI in the middle and the BLAST algorithm (cutadapt) will search for this sequence to identify UMIs. User can shorten the DNA sequence input to potentially obtain more reads with UMIs. The minimum sequence should ensure that the two flank regions (regions next to UMI, ***CATCTTACGATTACGCCAACCACTGCGG*** NNNNNTGNNNNN ***GACACATTCTCCCAGGCCCTACTT***) are unique in the amplicon.

### Customized variant calling analysis pipeline
The variant calling for each UMI group works in parallelism. One can modify the ***variants_calling.py*** file to implement customized data analysis pipeline. By default, VAULT applies bcftools for SNV and InDel calling, and Sniffles for SV calling. To change the variant calling pipeline, user can modify the ***cmd*** variable in ***snp_calling*** and ***sv_calling*** as you want. The ***cmd*** is in bash command format. VAULT provides "***save***(save_path), ***thread***, ***ax***(alignment mode of minimap2), ***refer***(reference sequence), ***fastq***(input fastq of UMI group), ***name***(the prefix of file)" for those bash commands.

### --group_filter
The UMI group filter will eliminate problematic UMI groups by surveying the read consistency within every UMI group. A UMI group is defined as a bin of reads with the same UMI. In theory, reads within the UMI group represent the same original molecule, thus exist the same sequence. However, the sequencing errors in UMI region will lead to a wrong separation of reads, which result in some UMI groups with reads from different molecules. A feature of such groups is, after SNV calling, they existed SNVs with various allelic frequencies. The UMI group filter analyzed the allelic frequencies of SNVs detected in every UMI group, and filter UMI groups based on variant read number, SNV number, and SNV allelic frequency. The filter will remove potentially problematic UMI groups to improve the confidence of final results, while sacrifice the number of usable UMI groups. Besides, there is no guarantee that all problematic UMI groups will be removed. Our in-house test showed that the UMI group filter will remove 10% to 40% UMI groups for different data sets. The repeatability of results was improved in random read sampling experiments with ***--group_filter***.

### Generate a summary of the analysis
Try the below command to generate a summary of the run, it includes several useful number as shown below. Please use the same parameters as the previous vault analysis. You can find your previous vault command in ***./example/result/date_time_vault.log***
> vault summarize -q ./example/nanopore_reads.fastq.gz -s ./example/result -r ./example/reference.fa --unmapped_reads -T 0.01

Output:
```
raw_read_number is: 252
used_read_number is: 246
reads_with_umi is: 246
detected_molecule_number is: 24          # detected UMI groups(molecule)
detected_passed_molecule_number is: 21   # filtered-passed UMI groups(molecule)
refer_seq_length is: 7077
covered_region_of_molecule(avg,median,min,max) is: 5784.67,6858,1582,7077   # length of regions with >=3 depth in each UMI group
p95_coverage_molecule is: 12             # UMI groups(molecule) with more than 95% of regions covered by >=3 depth
molecule_with_snv is: 16
total_snv_number is: 63
unique_snv_number is: 14
normalized_snv_number_per_SNVContainingMolecule(avg,median,min,max) is: 4.00683,4.04804,1.00493,9.00194
total_somatic_snv_number is: 0           # somatic SNV is defined as SNVs with VAF < $threshold (defined by -T)
unique_somatic_snv_number is: 0
somatic_snv_load_per_Mbp is: 0
molecule_with_sv is: 16/76.19%
total_sv_number is: 32
unique_sv_number is: 13
molecule_with_deletion is: 15/71.43%     # Below shows only deletion, insertion, inversion and duplication. For more information of SVs, please check ./example/result/snp/summary/all_sv.1count.2pos.3type.4length.txt
total_deletion is: 30
molecule_with_insertion is: 1/4.76%
total_insertion is: 2
molecule_with_inversion is: 0/0.00%
total_inversion is: 0
molecule_with_duplication is: 0/0.00%
total_duplication is: 0
```


