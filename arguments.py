import argparse
import os
import sys
from _version import vault_version
from tools import (
    umi_group_filter,
    call_consensus,
    change_VCF_pos,
    draw_circos,
    vaf_calculator,
    summarize_data,
    per_umi_variant,
    vdj_mrd,
)


def validate_file(x):
    if not os.path.exists(x):
        raise argparse.ArgumentTypeError("{0} does not exist".format(x))
    return x


def positive_int(x):
    value = int(x)
    if value < 1:
        raise argparse.ArgumentTypeError("must be an integer >= 1")
    return value


def percentage(x):
    value = float(x)
    if value < 0 or value > 100:
        raise argparse.ArgumentTypeError("must be between 0 and 100")
    return value


def get_argparse():
    parser = argparse.ArgumentParser(
        description='VAULT: analysis of UMI-labeled reads in IDMseq and iMiGseq.',
        epilog='See https://github.com/milesjor/VAULT for more details.'
    )

    subparsers = parser.add_subparsers(
        title='subcommands',
        description='valid subcommands',
        help='additional help',
        dest='subcommands_name'
    )
    parser.add_argument('-v', '--version', action='version',
                        version=vault_version, help='show the current version')

    # ======================
    # main VAULT required options
    # ======================
    req = parser.add_argument_group(title='Required options')
    req.add_argument('-u', '--umi_adapter', type=str,
                     help='UMI adapter sequence; NNNATGCNNN-like region is automatically treated as UMI')
    req.add_argument('-s', '--save_path', type=str,
                     help='path/to/save/ (VAULT result folder)')
    req.add_argument('-r', '--refer', type=validate_file,
                     help='path/to/ref.fa (reference sequence; not required when --only_group_reads is set)')
    req.add_argument('-q', '--fastq', type=validate_file,
                     help='path/to/reads.fastq or .fastq.gz')

    # ======================
    # main VAULT optional options
    # ======================
    opt = parser.add_argument_group(title='Optional options')
    opt.add_argument('-e', '--error', type=str, default='0.08',
                     help='error tolerance rate in UMI analysis (depends on sequencing error rate) [0.08]')
    opt.add_argument('-t', '--thread', type=int, default=5,
                     help='number of parallel workers (Python/multiprocessing) [5]')
    opt.add_argument(
        '-T', '--threshold', type=int, default=None,
        help='minimal read number per UMI group [5 for full analysis; '
             '2 with --only_group_reads]'
    )
    opt.add_argument('-b', '--bash_thread', type=int, default=1,
                     help='thread count for bash-based tools (minimap2, bcftools, sniffles, etc.) [1]')
    opt.add_argument('-F', '--allele_freq', type=float, default=0.67,
                     help='allele frequency cutoff for SNV filtering in final aggregation [0.67]')
    opt.add_argument('-f', '--sv_freq', type=float, default=0.5,
                     help='allele frequency cutoff for SV filtering in final aggregation [0.5]')
    opt.add_argument('-p', '--pe_fastq', type=str,
                     help='read2.fastq for Illumina paired-end sequencing')
    opt.add_argument('-a', '--align_mode', type=str, default='map-ont',
                     choices=['lr:hq', 'sr', 'map-ont', 'map-pb'],
                     help='alignment mode for minimap2: -ax [lr:hq|sr|map-ont|map-pb] [map-ont]')
    opt.add_argument('--minlength', type=int, default=0,
                     help='filter FASTQ reads shorter than this length [0: disabled]')
    opt.add_argument('--maxlength', type=int, default=0,
                     help='filter FASTQ reads longer than this length [0: disabled]')
    opt.add_argument('--unmapped_reads', action='store_true',
                     help='extract mapped reads (via minimap2) before UMI analysis')
    opt.add_argument('--group_filter', action='store_true',
                     help='filter out low-confidence UMI groups in final aggregation')
    opt.add_argument('--silence', action='store_true',
                     help='hide per-UMI variant-calling progress messages')
    opt.add_argument(
        '--only_group_reads',
        dest='only_group_reads',
        action='store_true',
        help='only perform UMI analysis and per-UMI FASTQ grouping; '
             'skip downstream per-UMI variant calling and final summarization'
    )
    opt.add_argument(
        '--only_group_fastq',
        dest='only_group_fastq',
        action='store_true',
        help=argparse.SUPPRESS
    )

    opt.add_argument(
        '--cutadapt_path',
        type=str,
        help='path to a specific cutadapt executable; '
             'if not provided, VAULT will try $VAULT_CUTADAPT27, then cutadapt in PATH'
    )
    opt.add_argument(
        '--sv_caller',
        choices=['sniffles1', 'sniffles2'],
        default='sniffles1',
        help='SV caller command style for long-read modes [sniffles1]'
    )
    opt.add_argument(
        '--sniffles_path',
        type=str,
        help='path to a specific sniffles executable'
    )
    opt.add_argument(
        '--sniffles_min_support',
        type=int,
        default=3,
        help='minimum supporting reads for Sniffles SV calling [3]'
    )

    # ======================
    # subcommands: perumi
    # ======================

    perumi = subparsers.add_parser(
        'perumi',
        help='run per-UMI variant calling on existing grouped FASTQ files',
        description=(
            'Run per-UMI variant calling on pre-generated grouped FASTQ files under '
            './grouped_reads/perfect_umi, and aggregate VCFs, without re-doing UMI extraction or grouping.'
        )
    )
    perumi.set_defaults(func=per_umi_variant.perumi_main)
    perumi.add_argument('-s', '--save_path', type=str, required=True,
                        help='path/to/VAULT_result_folder/ (folder containing grouped_reads/perfect_umi)')
    perumi.add_argument('-r', '--refer', type=validate_file, required=True,
                        help='path/to/ref.fa (same reference as in the main VAULT run)')
    perumi.add_argument('-t', '--thread', type=int, default=5,
                        help='number of parallel workers over UMI groups [5]')
    perumi.add_argument('-b', '--bash_thread', type=int, default=1,
                        help='thread count for minimap2/bcftools/sniffles inside each worker [1]')
    perumi.add_argument('-a', '--align_mode', type=str, default='map-ont',
                        choices=['sr', 'map-ont', 'map-pb'],
                        help='minimap2 -ax mode used for per-UMI variant calling [map-ont]')
    perumi.add_argument('-F', '--allele_freq', type=float, default=0.67,
                        help='allele frequency cutoff for SNV filtering in final aggregation [0.67]')
    perumi.add_argument('-f', '--sv_freq', type=float, default=0.5,
                        help='allele frequency cutoff for SV filtering in final aggregation [0.5]')
    perumi.add_argument('--group_filter', action='store_true',
                        help='apply UMI group filter in final_clean_up (same behavior as --group_filter in main VAULT)')
    perumi.add_argument('--sv_caller', choices=['sniffles1', 'sniffles2'],
                        default='sniffles1',
                        help='SV caller command style for long-read modes [sniffles1]')
    perumi.add_argument('--sniffles_path', type=str,
                        help='path to a specific sniffles executable')
    perumi.add_argument('--sniffles_min_support', type=int, default=3,
                        help='minimum supporting reads for Sniffles SV calling [3]')
    perumi.add_argument('--silence', action='store_true',
                        help='hide per-UMI variant-calling progress messages')

    # ======================
    # subcommand: vdj
    # ======================
    vdj = subparsers.add_parser(
        'vdj',
        help='detect IGH VDJ clonotypes (MRD) using IgBLAST',
        description=(
            'Detect IGH VDJ clonotypes from raw FASTQ or VAULT grouped reads. '
            'Designed for MRD workflows with UMI-tagged JH primer. '
            'Per-UMI clonotypes are called by majority vote over per-read IgBLAST AIRR results. '
            'This subcommand requires -u/--umi_adapter (full primer with N-pattern).'
        )
    )
    vdj.set_defaults(func=vdj_mrd.vdj_main)

    # ---- basic IO / mode ----
    vdj.add_argument('-s', '--save_path', type=str, required=True,
                     help='path/to/VAULT_result_folder/ (used to locate grouped_reads; also default output base)')
    vdj.add_argument('-o', '--out_dir', type=str, default=None,
                     help='output folder [default: <save_path>/vdj]')
    vdj.add_argument('--mode', type=str, default='auto',
                     choices=['auto', 'perumi', 'raw'],
                     help='auto: use perumi if grouped_reads exists; otherwise raw [auto]')
    vdj.add_argument('-q', '--fastq', type=validate_file, default=None,
                     help='raw FASTQ(.gz) for --mode raw (required in raw mode)')
    vdj.add_argument('--grouped_dir', type=str, default=None,
                     help='override grouped FASTQ directory [default: <save_path>/grouped_reads/perfect_umi]')

    # ---- REQUIRED: UMI primer (full sequence with N-pattern) ----
    vdj.add_argument(
        '-u', '--umi_adapter', type=str, required=True,
        help=(
            'Full UMI primer sequence used in VAULT -u, e.g. '
            'CATCTTACGATTACGCCAACCACTGNNNNNTGNNNNNCTTACCTGAGGAGACGGTGACC . '
            'vdj_mrd will derive 5\' adapter and 3\' (revcomp) adapter from this.'
        )
    )

    # ---- optional adapter overrides (now implemented in vdj_mrd.py) ----
    vdj.add_argument('--adapter_5', type=str, default=None,
                     help='Override 5\' adapter used for trimming (default: derived from -u)')
    vdj.add_argument('--adapter_3', type=str, default=None,
                     help='Override 3\' adapter used for trimming (default: revcomp-derived from -u)')
    # vdj.add_argument('--umi_len', type=int, default=None,
    #                  help='Deprecated!! Override UMI length (for logging/sanity; trimming uses full primer so not required)')

    # ---- cutadapt options ----
    vdj.add_argument('--cutadapt_path', type=str,
                     help='path to a specific cutadapt executable/script; '
                          'if not provided, VAULT will try $VAULT_CUTADAPT27, then cutadapt in PATH')
    vdj.add_argument('--cutadapt_error_rate', type=float, default=0.10,
                     help='cutadapt -e error rate for primer trimming [0.10]')
    vdj.add_argument('--cutadapt_revcomp', dest='cutadapt_revcomp',
                     action='store_true', default=True,
                     help='also consider reverse-complement orientation when trimming (ONT recommended) [enabled]')
    vdj.add_argument('--no_cutadapt_revcomp', dest='cutadapt_revcomp',
                     action='store_false', help=argparse.SUPPRESS)

    # ---- IgBLAST core ----
    vdj.add_argument('--igblastn', type=str, default='igblastn',
                     help='igblastn executable name or full path [igblastn]')
    vdj.add_argument('--germV', type=str, required=True,
                     help='IgBLAST V germline BLAST DB prefix (e.g. /path/human_IGHV)')
    vdj.add_argument('--germD', type=str, required=True,
                     help='IgBLAST D germline BLAST DB prefix (e.g. /path/human_IGHD)')
    vdj.add_argument('--germJ', type=str, required=True,
                     help='IgBLAST J germline BLAST DB prefix (e.g. /path/human_IGHJ)')
    vdj.add_argument('--auxiliary_data', type=validate_file, required=True,
                     help='IgBLAST auxiliary data file (e.g. human_gl.aux)')
    vdj.add_argument('--organism', type=str, default='human',
                     help='IgBLAST organism [human]')
    vdj.add_argument('--threads', type=int, default=8,
                     help='IgBLAST threads [-num_threads] [8]')

    # ---- per-group calling filters / knobs ----
    vdj.add_argument('--require_productive', action='store_true', default=False,
                     help='If set, clonotype calling and consensus extraction start from productive rearrangements '
                          'only (stop_codon=F, vj_in_frame=T, productive=T). Default: [False]')
    vdj.add_argument('--min_reads_per_umi', type=int, default=1,
                     help='skip UMI groups with reads < this threshold [1]')
    vdj.add_argument('--min_support', type=int, default=1,
                     help='minimum supporting reads for a clonotype within a UMI group [1]')
    vdj.add_argument('--top_k', type=int, default=50,
                     help='report top K clonotypes [50]')
    vdj.add_argument('--max_reads_per_group', type=int, default=0,
                     help='cap reads per UMI group when building merged query FASTQ [0: no cap]')

    # ---- identity thresholds ----
    vdj.add_argument('--min_v_identity', type=float, default=0.0,
                     help='minimum V alignment identity; accept 0~1 (fraction) or 0~100 (percent). 0 disables [0]')
    vdj.add_argument('--min_j_identity', type=float, default=0.0,
                     help='minimum J alignment identity; accept 0~1 (fraction) or 0~100 (percent). 0 disables [0]')

    # ---- clonotype key controls ----
    vdj.add_argument('--include_d', action='store_true',
                     help='include D gene in clonotype key (default: V+J+CDR3)')
    vdj.add_argument('--use_allele', action='store_true',
                     help='use allele-level calls (e.g. IGHV1-69*13) instead of gene-level (IGHV1-69)')
    vdj.add_argument('--cdr3_mode', type=str, default='aa',
                     choices=['aa', 'nt'],
                     help='use CDR3 amino-acid (aa) or nucleotide (nt) as clonotype tag [aa]')

    # ---- consensus calling pipeline controls ----
    vdj.add_argument('--do_consensus', action='store_true',
                     help='build consensus for selected clonotypes (amplicon + VDJ-only) (default by racon; optionally medaka)')
    vdj.add_argument(
        '--consensus_top_k',
        type=positive_int,
        default=None,
        help='limit consensus to the top K clonotypes by UMI count; '
             'if neither consensus selector is provided, K defaults to 3'
    )
    vdj.add_argument(
        '--consensus_top_p',
        type=percentage,
        default=None,
        help='call consensus only for clonotypes with UMI fraction greater '
             'than P percent; if neither consensus selector is provided, '
             'P defaults to 10'
    )
    vdj.add_argument('--racon_rounds', type=int, default=2,
                     help='number of racon polishing rounds [2]')
    vdj.add_argument('--consensus_threads', type=int, default=8,
                     help='threads for minimap2/racon/medaka in consensus stage [8]')
    vdj.add_argument('--consensus_fallback_min_reads', type=int, default=10,
                     help='When --do_consensus and --require_productive are set: '
                          'if a top clonotype has < N reads; will ignore --require_productive Default: [10]')
    vdj.add_argument('--enable_medaka', action='store_true',
                     help='run medaka polishing after racon (recommended for ONT)')
    vdj.add_argument('--medaka_model', type=str, default='auto',
                     help='medaka model name; "auto" lets vdj_mrd choose/skip model-specific handling [auto]')

    # NEW: allow specifying medaka_consensus path
    vdj.add_argument('--medaka_consensus', type=str, default=None,
                     help='path to medaka_consensus executable; if not provided, use medaka_consensus in PATH')

    # ======================
    # subcommands: summarize
    # ======================

    summarize = subparsers.add_parser(
        'summarize',
        help='summarize the result of a VAULT run',
        description=(
            'Summarize key metrics from a finished VAULT analysis, including read counts, '
            'molecule counts, SNV/SV loads, and coverage statistics.'
        )
    )
    summarize.set_defaults(func=summarize_data.summarize_main)
    summarize.add_argument('-s', '--save_path', type=str, required=True,
                          help='path/to/VAULT_result_folder/ (the same folder given to vault -s)')
    summarize.add_argument('-r', '--refer', type=validate_file, required=True,
                          help='path/to/ref.fa (reference sequence used in VAULT analysis)')
    summarize.add_argument('-q', '--fastq', type=validate_file, required=True,
                          help='path/to/raw_reads.fastq(.gz) analyzed by VAULT')
    summarize.add_argument('-T', '--somatic_VAF', type=float, default=0.01,
                          help='VAF threshold to define a somatic SNV [0.01]')
    summarize.add_argument('--unmapped_reads', action='store_true',
                          help='use this if you also used --unmapped_reads in the main VAULT analysis')

    # ======================
    # subcommands: consensus
    # ======================

    consensus = subparsers.add_parser(
        'consensus',
        help='get consensus sequence from VAULT result',
        description=(
            'Generate consensus sequences from VAULT per-UMI results using canu + medaka. '
            'You only need to provide the VAULT result folder (-s), which contains '
            './per_umi_process, ./grouped_reads, and ./umi_analysis.'
        )
    )
    consensus.set_defaults(func=call_consensus.consensus_main)
    consensus.add_argument('-s', '--save_path', type=str, required=True,
                           help='path/to/VAULT_result_folder/')
    consensus.add_argument('-t', '--thread', type=int, default=6,
                           help='number of parallel workers [6]')
    consensus.add_argument('-T', '--sub_thread', type=int, default=4,
                           help='thread count for each consensus worker [4]')
    consensus.add_argument('--threshold', type=int, default=20,
                           help='process UMI groups with read count >= this threshold [20]')

    # ======================
    # subcommands: position
    # ======================

    position = subparsers.add_parser(
        'position',
        help='liftover VAULT VCF coordinates to a standard reference',
        description=(
            'Correct VAULT VCF coordinates by aligning the VAULT analysis '
            'reference to a standard reference genome.'
        )
    )
    position.set_defaults(func=change_VCF_pos.position_main)
    position.add_argument('-v', '--vcf', dest='vcf_file', type=validate_file,
                          help='path/to/file.vcf or file.vcf.gz')
    position.add_argument('--vault_result', type=str,
                          help='path/to/VAULT_result_folder; convert VCFs directly under per_umi_process')
    position.add_argument('-c', '--chr_name', type=str,
                          help='force all converted CHROM/CHR2 names to this value')
    position.add_argument('-s', '--save_path', type=str,
                          help='optional output directory')
    position.add_argument('--analysis_ref', type=validate_file, required=True,
                          help='reference FASTA used by VAULT analysis')
    position.add_argument('--standard_ref', type=validate_file, required=True,
                          help='standard reference FASTA for downstream annotation')
    position.add_argument('--minimap2_path', type=str, default='minimap2',
                          help='minimap2 executable path [minimap2]')

    # ======================
    # subcommands: circos
    # ======================

    circos = subparsers.add_parser(
        'circos',
        help='prepare data for circos (used in iMiGseq)',
        description=(
            'Prepare data for circos plots. For VCF files, it will correct positions and remove indels. '
            'For depth files, it will bin coverage by user-defined bin size.'
        )
    )
    circos.set_defaults(func=draw_circos.circos_main)
    circos.add_argument('-n', '--name', type=str, default="circos",
                        help='prefix of output files [circos]')
    circos.add_argument('-d', '--depth_file', type=validate_file,
                        help='depth file generated by `samtools depth`')
    circos.add_argument('-s', '--save_path', type=str, required=True,
                        help='path/to/save/')
    circos.add_argument('-v', '--vcf_file', type=validate_file,
                        help='path/to/file.vcf')
    circos.add_argument('-b', '--bin_size', type=int, default=30,
                        help='bin size in base pairs for depth binning [30]')
    circos.add_argument('-g', '--genome', type=str,
                        choices=[
                            'mmt_nod_F6_10N_to_C57',
                            'hchrM.F9.UMIs',
                            'hchrM.F6.UMIs',
                            'hchrM.bamh1.UMIs',
                            'mmt_c57_F6_10N'
                        ],
                        help='predefined genome used in VAULT (for automatic chr_name handling)')
    circos.add_argument('-c', '--chr_name', type=str,
                        help='chromosome name shown in VCF file. If -g is provided, '
                             '-c will be set automatically (chrM for human, chrMT for mouse).')
    circos.add_argument('-A', '--keep_1alt', action='store_true',
                        help='keep only 1 ALT allele per site in VCF (for easier SNV annotation)')

    group_filter = subparsers.add_parser(
        'filter',
        help='filter out low-confidence UMI groups',
        description=(
            'Filter out low-confidence UMI groups after VAULT analysis. '
            'This is equivalent to `vault --group_filter` applied post-hoc.'
        )
    )
    group_filter.set_defaults(func=umi_group_filter.filter_main)
    group_filter.add_argument('-s', '--save_path', type=str, required=True,
                              help='path/to/per_umi_process (folder generated by VAULT)')
    group_filter.add_argument('-F', '--allele_freq', type=float, default=0.67,
                              help='allele frequency cutoff for SNV filtering [0.67]')

    vaf = subparsers.add_parser(
        'vaf',
        help='calculate variant allele frequency (VAF) based on UMI group counts',
        description=(
            'Calculate VAF for individual SNVs: '
            'VAF = (UMI group count with the SNV) / (UMI groups covering the SNV position).'
        )
    )
    vaf.set_defaults(func=vaf_calculator.get_vaf)
    vaf.add_argument('-s', '--save_path', type=str, required=True,
                     help='path/to/per_umi_process (folder generated by VAULT)')
    vaf.add_argument('-n', '--name', type=str, default='VAFcalculator',
                     help='output file name prefix [VAFcalculator]')

    # ======================
    # parse & dispatch
    # ======================
    args = parser.parse_args()

    if getattr(args, "only_group_fastq", False):
        args.only_group_reads = True
    args.only_group_fastq = bool(getattr(args, "only_group_reads", False))

    if args.subcommands_name is None and args.threshold is None:
        args.threshold = 2 if args.only_group_reads else 5

    if args.subcommands_name in ["summarize", "consensus", "position", "circos",
                                 "filter", "vaf", "perumi", "vdj"]:
        args.func(args)
        sys.exit(0)
    else:
        # main VAULT run
        if args.pe_fastq is not None:
            args.align_mode = 'sr'

        require_refer = not bool(getattr(args, "only_group_reads", False))

        missing = []
        if args.umi_adapter is None:
            missing.append("-u/--umi_adapter")
        if args.save_path is None:
            missing.append("-s/--save_path")
        if require_refer and args.refer is None:
            missing.append("-r/--refer")
        if args.fastq is None:
            missing.append("-q/--fastq")

        if missing:
            sys.stderr.write(
                "usage: vault [-h] [-v] -u UMI_ADAPTER -s SAVE_PATH "
                + ("-r REFER " if require_refer else "")
                + "-q FASTQ\n"
                "             [-e ERROR] [-t THREAD] [-T THRESHOLD] [-b BASH_THREAD]\n"
                "             [-F ALLELE_FREQ] [-f SV_FREQ] [-p PE_FASTQ]\n"
                "             [-a {sr,map-ont,map-pb}]\n"
                "             [--minlength MINLENGTH] [--maxlength MAXLENGTH]\n"
                "             [--cutadapt_path CUTADAPT_PATH]\n"
                "             [--sv_caller {sniffles1,sniffles2}]\n"
                "             [--sniffles_path SNIFFLES_PATH]\n"
                "             [--unmapped_reads] [--group_filter]\n"
                "             [--only_group_reads]\n"
                "subcommands: \n"
                "       vault {summarize,consensus,position,circos,filter,vaf,perumi,vdj} ...\n\n"
                "ERROR: the following arguments are required: "
                + ", ".join(missing) + "\n"
            )
            sys.exit(1)

    return args
