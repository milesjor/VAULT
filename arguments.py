import argparse
import os
import sys
from _version import vault_version
from tools import umi_group_filter, call_consensus, change_VCF_pos, draw_circos, vaf_calculator, summarize_data


def get_argparse():
    args = argparse.ArgumentParser(description='This is for analyzing UMI labeled reads in IDMseq and iMiGseq.',
                                   epilog='Refer github https://github.com/milesjor/VAULT for more detail')
    subparsers = args.add_subparsers(title='subcommands',
                                     description='valid subcommands',
                                     help='additional help',
                                     dest='subcommands_name')
    args.add_argument('-v', '--version', action='version', version=vault_version, help='show the current version')

    req = args.add_argument_group(title='Required options')
    req.add_argument('-u', '--umi_adapter', type=str, help='UMI sequence, automatically detect NNNATGCNNN as UMI')
    req.add_argument('-s', '--save_path', type=str, help='path/to/save/')
    req.add_argument('-r', '--refer', type=validate_file, help='path/to/ref.fa')
    req.add_argument('-q', '--fastq', type=validate_file, help='path/to/reads.fastq, or fastq.gz')

    opt = args.add_argument_group(title='Optional options')
    opt.add_argument('-e', '--error', type=str, default='0.11', help='error tolerate rate in umi analysis '
                                                                     '(depend on read error rate) [0.11]')
    opt.add_argument('-t', '--thread', type=int, default='5', help='parallel process number [5]')
    opt.add_argument('-T', '--threshold', type=int, default='5', help='threshold of read number for snp analysis [5]')
    opt.add_argument('-b', '--bash_thread', type=int, default='1', help='thread for running bash cmd [1]')
    opt.add_argument('-F', '--allele_freq', type=float, default='0.67', help='filter SNVs '
                                                                             'by allelic frequency [0.67]')
    opt.add_argument('-f', '--sv_freq', type=float, default='0.5', help='filter SVs by allelic frequency [0.5]')
    opt.add_argument('-p', '--pe_fastq', type=str, help='read2.fastq for illumina pair-end sequencing')
    opt.add_argument('-a', '--align_mode', type=str, default='map-ont',
                     choices=['sr', 'map-ont', 'map-pb'], help='parameter in alignment, '
                                                               'minimap2 -ax [sr|map-ont|map-pb]')
    opt.add_argument('--minlength', type=int, default='0',
                     help='filter fastq file to remove reads with length less than [int]')
    opt.add_argument('--maxlength', type=int, default='0',
                     help='filter fastq file to remove reads with length more than [int]')
    opt.add_argument('--unmapped_reads', action='store_true', help='extract mapped reads before UMI analysis')
    opt.add_argument('--group_filter', action='store_true', help='filter out low-confidence UMI groups')
    opt.add_argument('--py_only', action='store_true', help='use only python script in the umi analysis')

    # subcommands
    # summarize the result
    summarize = subparsers.add_parser('summarize',
                                      help='summarize the result of VAULT',
                                      description='This will generate several useful number to give an overview of the sequencing data')
    summarize.set_defaults(func=summarize_data.summarize_main)
    summarize.add_argument('-s', '--save_path', type=str, required=True, help='path/to/VAULT_result_folder/')
    summarize.add_argument('-r', '--refer', type=validate_file, required=True,
                           help='path/to/ref.fa, reference sequence in VAULT analysis')
    summarize.add_argument('-q', '--fastq', type=validate_file, required=True,
                           help='path/to/reads.fastq, raw sequencing reads analyzed by VAULT')
    summarize.add_argument('-T', '--somatic_VAF', type=float, default='0.1', help='VAF threshold to define a somatic SNV')
    summarize.add_argument('--unmapped_reads', action='store_true', help='use --unmapped_reads here if you used it in vault analysis')


    # get consensus sequence
    consensus = subparsers.add_parser('consensus',
                                      help='get consensus sequence from VAULT result',
                                      description='This is for generating consensus sequence from VAULT result '
                                                  'using canu + medaka. '
                                                  'Just input the path of VAULT result folder '
                                                  '--> the folder of vault -s (containing ./snp ./grouped_reads ./umi_analysis). ')
    consensus.set_defaults(func=call_consensus.consensus_main)
    consensus.add_argument('-s', '--save_path', type=str, required=True, help='path/to/VAULT_result_folder/')
    consensus.add_argument('-t', '--thread', type=int, default='6', help='parallel worker [6]')
    consensus.add_argument('-T', '--sub_thread', type=int, default='4', help='thread for every parallel worker [4]')
    consensus.add_argument('--threshold', type=int, default='20', help='process UMI group with reads more than [20]')

    # change vcf position
    position = subparsers.add_parser('position',
                                     help='correct the position and chr_name in vcf file',
                                     description='It will correct the position and chr_name in vcf file '
                                                 'to enable variant annotation. It can also reverse the coordinate'
                                                 ' and DNA base if the used referene in VAULT '
                                                 'is reverse complimentary of the referene genome.')
    position.set_defaults(func=change_VCF_pos.position_main)
    position.add_argument('-v', '--vcf_file', type=validate_file, required=True, help='path/to/file.vcf')
    position.add_argument('-c', '--chr_name', type=str, required=True, help='chromosome name')
    position.add_argument('-p', '--pos_change', type=str, required=True, help='position change, e.g. +12300 or -55789')
    position.add_argument('-b', '--total_base', type=int,
                          help='the length of reference genome used. When provided, it will reverse the '
                               'coordinate(position) and also do reverse complimentary of the DNA bases in vcf file')
    position.add_argument('-s', '--save_path', type=str, default="./", help='path/to/save/')

    # prepare data for circos
    circos = subparsers.add_parser('circos',
                                   help='prepare data for circos, used in iMiGseq',
                                   description='This is for preparing data for circos. It is used in iMiGseq. '
                                               'For vcf file, it will correct position and remove indel. '
                                               'For depth file, it will bin depth based on '
                                               'user defined bin size.')
    circos.set_defaults(func=draw_circos.circos_main)
    circos.add_argument('-n', '--name', type=str, default="circos", help='prefix of output file [circos]')
    circos.add_argument('-d', '--depth_file', type=validate_file, help='depth file from samtools depth')
    circos.add_argument('-s', '--save_path', type=str, required=True, help='path/to/save/')
    circos.add_argument('-v', '--vcf_file', type=validate_file, help='path/to/file.vcf')
    circos.add_argument('-b', '--bin_size', type=int, default="30", help='how many bases per bin [30]')
    circos.add_argument('-g', '--genome', type=str,
                        choices=['mmt_nod_F6_10N_to_C57', 'hchrM.F9.UMIs', 'hchrM.F6.UMIs', 'hchrM.bamh1.UMIs',
                                 'mmt_c57_F6_10N'],
                        help='the genome used in VAULT')
    circos.add_argument('-c', '--chr_name', type=str,
                        help='chromosome name showed in vcf file, '
                             'when providing -g, -c will be set automatically as [chrM] for human and [chrMT] for mouse')
    circos.add_argument('-A', '--keep_1alt', action='store_true',
                        help='leave only 1 Alt in vcf file for helping SNV annotation')

    # filter out low-confidence UMI groups
    group_filter = subparsers.add_parser('filter',
                                         help='filter out low-confidence UMI groups',
                                         description='This is for filtering out low-confidence UMI groups after '
                                                     'VAULT analysis. It is the same as <vault --group_filter>')
    group_filter.set_defaults(func=umi_group_filter.filter_main)
    group_filter.add_argument('-s', '--save_path', type=str, required=True,
                              help='path/to/snp (the path to snp folder generated by VAULT)')
    group_filter.add_argument('-F', '--allele_freq', type=float, default='0.67', help='filter SNVs '
                                                                                      'by allelic frequency [0.67]')

    # calculate VAF
    vaf = subparsers.add_parser('vaf',
                                help='calculate variant allele frequency (VAF) based on UMI group number',
                                description='calculate VAF for individual SNV, '
                                            'VAF = (UMI group number with the same SNV) / '
                                            '(UMI groups number covering the SNV position)')
    vaf.set_defaults(func=vaf_calculator.get_vaf)
    vaf.add_argument('-s', '--save_path', type=str, required=True,
                     help='path/to/snp (the path to snp folder generated by VAULT)')
    vaf.add_argument('-n', '--name', type=str, default='VAFcalculator', help='File name prefix [VAFcalculator]')

    args = args.parse_args()
    if args.subcommands_name in ["summarize", "consensus", "position", "circos", "filter", "vaf"]:
        args.func(args)

    else:
        if args.pe_fastq is not None:
            args.align_mode = 'sr'
        if args.umi_adapter is None or args.save_path is None or args.refer is None or args.fastq is None:
            sys.stderr.write("usage: vault [-h] [-v] -u UMI_ADAPTER -s SAVE_PATH -r REFER -q FASTQ\n"
                             "             [-e ERROR] [-t THREAD] [-T THRESHOLD] [-b BASH_THREAD]\n"
                             "             [-F ALLELE_FREQ] [-f SV_FREQ] [-p PE_FASTQ]\n"
                             "             [-a {sr,map-ont,map-pb}] [--minlength MINLENGTH]\n"
                             "             [--maxlength MAXLENGTH] [--unmapped_reads] [--group_filter]\n"
                             "             {summarize,consensus,position,circos,filter,vaf} ...\n"
                             "vault: error: the following arguments are required: "
                             "-u/--umi_adapter, -s/--save_path, -r/--refer, -q/--fastq\n")
            sys.exit(1)

    return args


def validate_file(x):
    if not os.path.exists(x):
        raise argparse.ArgumentTypeError("{0} does not exist".format(x))
    return x
