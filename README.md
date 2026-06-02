# VAULT
 Variant Analysis with UMI for Long-read Technology (VAULT)

VAULT is a tool for analyzing UMI-labeled reads, works for both error-prone long reads and accurate single-end/paired-end short reads.
<p align="center">
<img src="https://github.com/milesjor/VAULT/blob/master/example/pic/compare.png" width="600"/>
</p>

More detail: [Long-read Individual-molecule Sequencing Reveals CRISPR-induced Genetic Heterogeneity in Human ESCs](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-02143-8)

## INSTALLATION
### Prerequisites
Miniconda or Anaconda is required. Miniconda is recommended for a lightweight
installation. For Linux x86_64, install the latest Miniconda release:
```
mkdir -p "$HOME/miniconda3"
cd "$HOME/miniconda3"
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh \
        -b \
        -u \
        -p "$HOME/miniconda3"
rm Miniconda3-latest-Linux-x86_64.sh
"$HOME/miniconda3/bin/conda" init bash
source "$HOME/.bashrc"
```

### Download the VAULT package
```
git clone https://github.com/milesjor/VAULT.git
cd ./VAULT/
```

### Install all required modules
VAULT keeps only portable source code in the repository. Platform-specific or
machine-specific copies of cutadapt and Sniffles are not bundled. External
command-line tools should be installed in a conda environment, available in
`PATH`, or passed explicitly with command-line options.

Generic installation from the provided environment file:
```
conda env create --name vault --file ./tools/vault_env.yaml
conda activate vault
```

### Legacy tool installation
For the original VAULT/iMiGseq workflow, `cutadapt 2.7` and
`Sniffles1 1.0.12b` are the recommended legacy tools. These versions best match
the UMI-grouped-read design of VAULT and the behavior used when the pipeline was
originally tuned.

Newer versions, such as cutadapt 5.x and Sniffles2, can be used through
`--cutadapt_path`, `--sv_caller`, and `--sniffles_path`, but they are not the
preferred defaults for reproducing the original VAULT behavior. cutadapt 5.x may
change adapter-matching or output behavior in edge cases, and Sniffles2 is a
substantially redesigned SV caller optimized for modern long-read sample-level
or population-level SV analysis. VAULT instead performs per-UMI analysis, where
each grouped FASTQ may contain only a small number of reads derived from the
same original molecule. For this reason, Sniffles1 is usually a better match for
VAULT's small UMI-group SV calling strategy.

Linux: install `cutadapt 2.7` and `Sniffles1 1.0.12` directly in a
Python 3.7.8 conda environment:
```
conda create -n vault_legacy_py37 \
        -c conda-forge \
        -c bioconda \
        --override-channels \
        python=3.7.8 \
        cutadapt=2.7 \
        sniffles=1.0.12

conda activate vault_legacy_py37
cutadapt --version
sniffles -h | head

export VAULT_CUTADAPT27="$HOME/miniconda3/envs/vault_legacy_py37/bin/cutadapt"
export VAULT_SNIFFLES1="$HOME/miniconda3/envs/vault_legacy_py37/bin/sniffles"
"$VAULT_CUTADAPT27" --version
"$VAULT_SNIFFLES1" -h | head
```

macOS Apple Silicon: `cutadapt 2.7` is not available as a native `osx-arm64`
Bioconda package. Create an `osx-64` conda environment instead. Rosetta is
required to run the x86_64 executable on Apple Silicon.
```
CONDA_SUBDIR=osx-64 conda create -n vault_legacy_py37 \
        -c conda-forge \
        -c bioconda \
        --override-channels \
        python=3.7.8 \
        cutadapt=2.7

conda activate vault_legacy_py37
conda config --env --set subdir osx-64
cutadapt --version
export VAULT_CUTADAPT27="$HOME/miniconda3/envs/vault_legacy_py37/bin/cutadapt"
```

macOS Apple Silicon: `Sniffles1 1.0.12b` usually needs a local source build.
Install the binary inside the user's `VAULT/legacy_tools` folder. Homebrew
LLVM/OpenMP helps avoid Apple Clang OpenMP build failures.
```
export VAULT_DIR="$HOME/path/to/VAULT"
cd "$VAULT_DIR"

brew install llvm libomp cmake

mkdir -p legacy_tools
git clone https://github.com/fritzsedlazeck/Sniffles.git \
        legacy_tools/Sniffles-v1.0.12b
cd legacy_tools/Sniffles-v1.0.12b
git checkout v1.0.12b

perl -0pi -e \
        's/#if defined\(MACOS\) \|\| defined\(TARGET_OS_MAC\)/#if defined(MACOS)/' \
        lib/zlib-1.2.7/zutil.h

export CC=/opt/homebrew/opt/llvm/bin/clang
export CXX=/opt/homebrew/opt/llvm/bin/clang++
export PATH="/opt/homebrew/opt/llvm/bin:$PATH"

cmake \
        -S . \
        -B build \
        -DCMAKE_BUILD_TYPE=Release \
        -DCMAKE_POLICY_VERSION_MINIMUM=3.5 \
        -DCMAKE_C_COMPILER=/opt/homebrew/opt/llvm/bin/clang \
        -DCMAKE_CXX_COMPILER=/opt/homebrew/opt/llvm/bin/clang++ \
        -DCMAKE_CXX_STANDARD=11 \
        -DCMAKE_CXX_STANDARD_REQUIRED=ON \
        -DOpenMP_C_FLAGS="-Xpreprocessor -fopenmp -I/opt/homebrew/opt/libomp/include" \
        -DOpenMP_CXX_FLAGS="-Xpreprocessor -fopenmp -I/opt/homebrew/opt/libomp/include" \
        -DOpenMP_C_LIB_NAMES=omp \
        -DOpenMP_CXX_LIB_NAMES=omp \
        -DOpenMP_omp_LIBRARY=/opt/homebrew/opt/libomp/lib/libomp.dylib \
        -DCMAKE_EXE_LINKER_FLAGS="-L/opt/homebrew/opt/libomp/lib -lomp"

cmake --build build --parallel
mkdir -p ../sniffles1-v1.0.12b/bin
cp "$(find . -type f -perm -111 -name sniffles | head -n 1)" \
        ../sniffles1-v1.0.12b/bin/sniffles
cd ../..

export VAULT_SNIFFLES1="$VAULT_DIR/legacy_tools/sniffles1-v1.0.12b/bin/sniffles"
"$VAULT_SNIFFLES1" -h | head
```

You can keep these exports in the shell before running VAULT:
```
export VAULT_DIR="$HOME/path/to/VAULT"
export VAULT_CUTADAPT27="$HOME/miniconda3/envs/vault_legacy_py37/bin/cutadapt"
export VAULT_SNIFFLES1="$VAULT_DIR/legacy_tools/sniffles1-v1.0.12b/bin/sniffles"
python ./vault ... --sv_caller sniffles1
```

Or pass the executables explicitly for a single run:
```
python ./vault ... \
        --cutadapt_path "$HOME/miniconda3/envs/vault_legacy_py37/bin/cutadapt" \
        --sv_caller sniffles1 \
        --sniffles_path "$VAULT_DIR/legacy_tools/sniffles1-v1.0.12b/bin/sniffles"
```

For Sniffles2, use the same pattern with `--sv_caller sniffles2`,
`--sniffles_path /path/to/sniffles2`, or `$VAULT_SNIFFLES2`.

Check the core tools:
```
python ./vault --help
cutadapt --version
minimap2 --version
samtools --version
bcftools --version
seqtk 2>&1 | head
```

Without explicit paths, VAULT searches:
```
cutadapt:  --cutadapt_path, then $VAULT_CUTADAPT27, then cutadapt in PATH
Sniffles1: --sniffles_path, then $VAULT_SNIFFLES1, then sniffles in PATH
Sniffles2: --sniffles_path, then $VAULT_SNIFFLES2, then sniffles in PATH
```

On macOS, or whenever multiple Python environments are installed, explicit
Python invocation can be safer:
```
cd /path/to/VAULT
conda run -n vault python ./vault --help
```

Optional: add VAULT to your PATH
```
cat >> "$HOME/.bashrc" <<'EOF'

# VAULT added
# !!! -> modify below path to match your installation!
export VAULT_DIR="/path/to/VAULT"
export PATH="$PATH:$VAULT_DIR"
export VAULT_CUTADAPT27="$HOME/miniconda3/envs/vault_legacy_py37/bin/cutadapt"
export VAULT_SNIFFLES1="$HOME/miniconda3/envs/vault_legacy_py37/bin/sniffles"
# export VAULT_SNIFFLES1="$VAULT_DIR/legacy_tools/sniffles1-v1.0.12b/bin/sniffles"
EOF

source "$HOME/.bashrc"
```

## PACKAGE STRUCTURE
Current local VAULT folder layout:
```
VAULT/
├── vault                         # command-line entry script
├── vault_main.py                 # main VAULT workflow
├── arguments.py                  # argparse definitions and subcommands
├── variants_calling.py           # per-UMI SNV/SV calling commands
├── check_umi.sh                  # shell-based UMI extraction helper
├── tools/
│   ├── call_consensus.py         # consensus subcommand
│   ├── change_VCF_pos.py         # position subcommand
│   ├── draw_circos.py            # circos subcommand
│   ├── summarize_data.py         # summarize subcommand
│   ├── umi_group_filter.py       # filter subcommand
│   ├── vaf_calculator.py         # vaf subcommand
│   ├── per_umi_variant.py        # perumi subcommand
│   ├── vdj_mrd.py                # vdj subcommand
│   ├── cutadapt27.py             # cutadapt resolver/helper
│   ├── sniffles1.py              # Sniffles1/Sniffles2 resolver/helper
│   ├── extract_mapped_reads.py
│   ├── filter_sv.py
│   ├── read_length_filter.py
│   ├── snv_location.py
│   ├── vcf_to_matrix.py
│   ├── pS_pN_count.py
│   ├── vault_env.yaml
│   └── vertebrate_mitochondrial_code.txt
├── legacy_tools/                  # optional user-built legacy tools
│   └── sniffles1-v1.0.12b/bin/sniffles
└── example/
    ├── nanopore_reads.fastq.gz
    ├── reference.fa
    └── pic/
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

On the current macOS workstation, the safest non-interactive form is:
```
conda run -n py310 python ./vault --help
```

A list of available parameters:
```
usage: vault [-h] [-v] [-u UMI_ADAPTER] [-s SAVE_PATH] [-r REFER] [-q FASTQ]
             [-e ERROR] [-t THREAD] [-T THRESHOLD] [-b BASH_THREAD]
             [-F ALLELE_FREQ] [-f SV_FREQ] [-p PE_FASTQ]
             [-a {lr:hq,sr,map-ont,map-pb}] [--minlength MINLENGTH]
             [--maxlength MAXLENGTH] [--unmapped_reads] [--group_filter]
             [--silence] [--only_group_reads]
             [--cutadapt_path CUTADAPT_PATH]
             [--sv_caller {sniffles1,sniffles2}]
             [--sniffles_path SNIFFLES_PATH]
             [--sniffles_min_support SNIFFLES_MIN_SUPPORT]
             {summarize,consensus,position,circos,filter,vaf,perumi,vdj} ...

VAULT: analysis of UMI-labeled reads in IDMseq and iMiGseq.

options:
  -h, --help            show this help message and exit
  -v, --version         show the current version

subcommands:
  valid subcommands

  {perumi,vdj,summarize,consensus,position,circos,filter,vaf}
                        additional help
    perumi              run per-UMI variant calling on existing grouped FASTQ files
    vdj                 detect IGH VDJ clonotypes (MRD) using IgBLAST
    summarize           summarize the result of a VAULT run
    consensus           get consensus sequence from VAULT result
    position            correct position and chromosome name in a VCF file
    circos              prepare data for circos (used in iMiGseq)
    filter              filter out low-confidence UMI groups
    vaf                 calculate variant allele frequency (VAF) based on UMI group counts

Required options:
  -u UMI_ADAPTER, --umi_adapter UMI_ADAPTER
                        UMI adapter sequence; NNNATGCNNN-like region is
                        automatically treated as UMI
  -s SAVE_PATH, --save_path SAVE_PATH
                        path/to/save/ (VAULT result folder)
  -r REFER, --refer REFER
                        path/to/ref.fa (reference sequence; not required when
                        --only_group_reads is set)
  -q FASTQ, --fastq FASTQ
                        path/to/reads.fastq or .fastq.gz

Optional options:
  -e ERROR, --error ERROR
                        error tolerance rate in UMI analysis (depends on
                        sequencing error rate) [0.08]
  -t THREAD, --thread THREAD
                        number of parallel workers (Python/multiprocessing)
                        [5]
  -T THRESHOLD, --threshold THRESHOLD
                        minimal read number per UMI group for variant analysis
                        [5]
  -b BASH_THREAD, --bash_thread BASH_THREAD
                        thread count for bash-based tools (minimap2, bcftools,
                        sniffles, etc.) [1]
  -F ALLELE_FREQ, --allele_freq ALLELE_FREQ
                        allele frequency cutoff for SNV filtering in final
                        aggregation [0.67]
  -f SV_FREQ, --sv_freq SV_FREQ
                        allele frequency cutoff for SV filtering in final
                        aggregation [0.5]
  -p PE_FASTQ, --pe_fastq PE_FASTQ
                        read2.fastq for Illumina paired-end sequencing
  -a {lr:hq,sr,map-ont,map-pb}, --align_mode {lr:hq,sr,map-ont,map-pb}
                        alignment mode for minimap2: -ax
                        [lr:hq|sr|map-ont|map-pb] [map-ont]
  --minlength MINLENGTH
                        filter FASTQ reads shorter than this length [0:
                        disabled]
  --maxlength MAXLENGTH
                        filter FASTQ reads longer than this length [0:
                        disabled]
  --unmapped_reads      extract mapped reads (via minimap2) before UMI
                        analysis
  --group_filter        filter out low-confidence UMI groups in final
                        aggregation
  --silence             hide per-UMI variant-calling progress messages
  --only_group_reads    only perform UMI analysis and per-UMI FASTQ grouping;
                        skip downstream per-UMI variant calling and final
                        summarization
  --cutadapt_path CUTADAPT_PATH
                        path to a specific cutadapt executable; if not
                        provided, VAULT will try $VAULT_CUTADAPT27, then
                        cutadapt in PATH
  --sv_caller {sniffles1,sniffles2}
                        SV caller command style for long-read modes
                        [sniffles1]
  --sniffles_path SNIFFLES_PATH
                        path to a specific sniffles executable
  --sniffles_min_support SNIFFLES_MIN_SUPPORT
                        minimum supporting reads for Sniffles SV calling [3]

See https://github.com/milesjor/VAULT for more details.
```

### Subcommands
Current VAULT subcommands:

| Subcommand | Purpose |
| --- | --- |
| `perumi` | Run per-UMI variant calling on existing grouped FASTQ files without redoing UMI extraction/grouping. |
| `vdj` | Detect IGH VDJ clonotypes for MRD workflows using IgBLAST; optional consensus uses minimap2/racon and medaka. |
| `summarize` | Summarize a finished VAULT run, including read counts, molecule counts, SNV/SV loads, and coverage statistics. |
| `consensus` | Generate per-UMI consensus sequences from VAULT results using canu and medaka. |
| `position` | Correct chromosome name and coordinate offsets in a VCF file; optionally reverse coordinates and reverse-complement alleles. |
| `circos` | Prepare depth and VCF-derived files for circos-style iMiGseq plots. |
| `filter` | Filter low-confidence UMI groups after VAULT analysis. |
| `vaf` | Calculate SNV variant allele frequency from UMI group counts. |

To print help for each subcommand:
```
python ./vault summarize --help
python ./vault consensus --help
python ./vault position --help
python ./vault circos --help
python ./vault filter --help
python ./vault vaf --help
python ./vault perumi --help
python ./vault vdj --help
```

## EXAMPLE

### Command
```
python ./vault -u CATCTTACGATTACGCCAACCACTGCGGNNNNNTGNNNNNGACACATTCTCCCAGGCCCTACTT \
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

After a full main run finishes, VAULT automatically runs:
```
python ./vault summarize \
        -q ./example/nanopore_reads.fastq.gz \
        -s ./example/result \
        -r ./example/reference.fa \
        --unmapped_reads \
        -T 0.01
```
and prints the summary to the screen. Users can also run `vault summarize`
manually later with the same FASTQ, reference, result folder, and
`--unmapped_reads` setting used in the main run. The automatic summary is
skipped when `--only_group_reads` is used.

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
├── grouped_reads                           # FASTQ/read-name files for UMI groups
│   ├── perfect_umi                         # analyzed per-UMI FASTQ files (*.fastq.gz)
│   └── perfect_umi_unanalyzed              # UMI read-name lists below analysis threshold (*.lst)
├── per_umi_process  # folder for variant analysis
│   ├── all_snv_from_perfect_umi.vcf        # raw variant calling result for small variants
│   ├── all_snv_from_perfect_umi.pcent.vcf  # add variant allele frequency (supported reads percentage) in the info field of vcf file
│   ├── all_snv_from_perfect_umi.pcent.rem.vcf              # remove wrong UMI group
│   ├── all_snv_from_perfect_umi.pcent.rem.flt.vcf          # filter by depth, quality, allele frequency
│   ├── all_sv_from_perfect_umi.vcf         # raw variant calling result for large variants (>= 30bp)
│   ├── all_sv_from_perfect_umi.filtered.0.67.vcf           # FINAL SV result (remove wrong UMI group and filter SVs by allele frequency [0.5])
│   ├── all_sv_from_perfect_umi.filtered.0.67.sorted.vcf    # sort by position
│   ├── coverage.3plus.txt                  # The region with coverage >= 3 in each UMI groups, can be used to filter out UMI groups
│   ├── pass.group.lst                      # UMI groups that pass group_filter
│   ├── pass_snv_from_perfect_umi.flt.vcf   # FINAL SNV and InDel result (remove wrong UMI group and filter by depth, quality, allele frequency)
│   ├── umi_group.flt.summary.txt           # intermediate file in [--group_filter]
│   ├── wrong.group.lst                     # UMI groups that fail in [--group_filter]
│   ├── wrong.group.summary.txt             # intermediate file in [--group_filter] for wrong UMI groups
│   └── perfect_umi                         # individual UMI analysis result for 5' and 3' end of reads
└── umi_analysis
    ├── 3end_UMIs
    └── 5end_UMIs
```

### Individual UMI group folder
```
./result/per_umi_process/perfect_umi/
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

### Additional subcommand examples
Summarize a completed VAULT run:
```
python ./vault summarize \
        -q ./example/nanopore_reads.fastq.gz \
        -s ./example/result \
        -r ./example/reference.fa \
        --unmapped_reads \
        -T 0.01
```

Run per-UMI variant calling again from existing grouped FASTQ files. This is
useful when the main run finished UMI grouping but the variant-calling step
failed, or when you want to change SV/SNV caller settings without repeating UMI
extraction. Keep `grouped_reads/` and remove only the current variant-calling
folder (`per_umi_process`) before rerunning:
```
rm -rf ./example/result/per_umi_process
```

```
python ./vault perumi \
        -s ./example/result \
        -r ./example/reference.fa \
        -a map-ont \
        -t 5 \
        -b 1 \
        --group_filter \
        --sv_caller sniffles1 \
        --sniffles_path "$VAULT_DIR/legacy_tools/sniffles1-v1.0.12b/bin/sniffles"
```

Prepare circos input files:
```
python ./vault circos \
        -s ./example/circos \
        -d ./example/result/per_umi_process/all.coverage.3plus.length.txt \
        -v ./example/result/per_umi_process/pass_snv_from_perfect_umi.flt.vcf \
        -g hchrM.F6.UMIs
```

Correct VCF positions for downstream annotation:
```
python ./vault position \
        -v ./example/result/per_umi_process/pass_snv_from_perfect_umi.flt.vcf \
        -c chrM \
        -p +0 \
        -s ./example/position_corrected
```

Run IGH VDJ clonotype detection with IgBLAST. 
```
python ./vault vdj \
        -s /path/to/VAULT_result_folder \
        -u 'CATCTTACGATTACGCCAACCACTGNNNNNTGNNNNNCTTACCTGAGGAGACGGTGACC' \
        --igblastn /path/to/igblastn \
        --germV /path/to/human_IGHV \
        --germD /path/to/human_IGHD \
        --germJ /path/to/human_IGHJ \
        --auxiliary_data /path/to/human_gl.aux \
        --mode auto \
        --top_k 50
```

## TIPS

### 5' end UMI and 3' end UMI
<p align="center">
<img src="https://github.com/milesjor/VAULT/blob/master/example/pic/connect_two_end_umi.png" width="600"/>
</p>

### Tips on increasing the analyzable reads (reads with UMIs)
To potentially improve the usable reads, VAULT offers a parameter [-e] to define the error tolerant threshold. The default value [0.08] is set based on the average sequencing error rate of Nanopore. When increasing this value, we will have more reads with UMIs. But this may compromise the confidence of UMI groups.

The other parameter can be adjusted is [-u], it asks for a DNA sequence with UMI in the middle and the BLAST algorithm (cutadapt) will search for this sequence to identify UMIs. User can shorten the DNA sequence input to potentially obtain more reads with UMIs. The minimum sequence should ensure that the two flank regions (regions next to UMI, ***CATCTTACGATTACGCCAACCACTGCGG*** NNNNNTGNNNNN ***GACACATTCTCCCAGGCCCTACTT***) are unique in the amplicon.

### Customized variant calling analysis pipeline
The variant calling for each UMI group works in parallelism. One can modify the ***variants_calling.py*** file to implement customized data analysis pipeline. By default, VAULT applies bcftools for SNV and InDel calling, and Sniffles for SV calling. To change the variant calling pipeline, user can modify the ***cmd*** variable in ***snv_calling*** and ***sv_calling*** as you want. The ***cmd*** is in bash command format. VAULT provides "***save***(save_path), ***thread***, ***ax***(alignment mode of minimap2), ***refer***(reference sequence), ***fastq***(input fastq of UMI group), ***name***(the prefix of file)" for those bash commands.

### --group_filter
The UMI group filter will eliminate problematic UMI groups by surveying the read consistency within every UMI group. A UMI group is defined as a bin of reads with the same UMI. In theory, reads within the UMI group represent the same original molecule, thus exist the same sequence. However, the sequencing errors in UMI region will lead to a wrong separation of reads, which result in some UMI groups with reads from different molecules. A feature of such groups is, after SNV calling, they existed SNVs with various allelic frequencies. The UMI group filter analyzed the allelic frequencies of SNVs detected in every UMI group, and filter UMI groups based on variant read number, SNV number, and SNV allelic frequency. The filter will remove potentially problematic UMI groups to improve the confidence of final results, while sacrifice the number of usable UMI groups. Besides, there is no guarantee that all problematic UMI groups will be removed. Our in-house test showed that the UMI group filter will remove 10% to 40% UMI groups for different data sets. The repeatability of results was improved in random read sampling experiments with ***--group_filter***.

### Generate a summary of the analysis
After a full main run, VAULT automatically runs `vault summarize` and prints the
summary to the screen. You can also run the summary step manually later. Please
use the same FASTQ, reference, result folder, and `--unmapped_reads` setting as
the previous VAULT analysis. You can find your previous VAULT command in
`./example/result/date_time_vault.log`.

```
python ./vault summarize \
        -q ./example/nanopore_reads.fastq.gz \
        -s ./example/result \
        -r ./example/reference.fa \
        --unmapped_reads \
        -T 0.01
```

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
molecule_with_deletion is: 15/71.43%     # Below shows only deletion, insertion, inversion and duplication. For more information of SVs, please check ./example/result/snv/summary/all_sv.1count.2pos.3type.4length.txt
total_deletion is: 30
molecule_with_insertion is: 1/4.76%
total_insertion is: 2
molecule_with_inversion is: 0/0.00%
total_inversion is: 0
molecule_with_duplication is: 0/0.00%
total_duplication is: 0
```
