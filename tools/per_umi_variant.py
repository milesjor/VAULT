#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Run per-UMI variant calling on existing grouped FASTQ files.
#
# Usage example:
#   ./vault perumi -s ./example/result -r ./example/reference.fa \
#       -t 4 -b 1 -a map-ont --group_filter -F 0.67 -f 0.5
#

import sys
import os
import logging
from datetime import datetime


def perumi_main(args):
    """
    Run per-UMI variant calling on existing grouped FASTQ files only.

    This subcommand assumes that the main VAULT pipeline has already been run
    (at least up to per-UMI FASTQ grouping), and that the following structure exists:

        <save_path>/grouped_reads/perfect_umi/*.fastq.gz

    Steps:
      1) Call ProcessUmi.run_variant_on_grouped_fastq(args) to perform
         per-UMI variant calling in parallel for each UMI group.
      2) Call final_clean_up(args) to aggregate all per-UMI VCFs and apply
         UMI-group-based filtering (optional) and allele frequency filters.

    Output VCFs are written under:

        <save_path>/per_umi_process/perfect_umi/<group>/vcf/*.vcf
        <save_path>/per_umi_process/*.vcf  (aggregated)
    """
    if not os.path.isdir(args.save_path):
        sys.stderr.write("ERROR: save_path does not exist: %s\n" % args.save_path)
        sys.exit(1)

    ctime = datetime.now().strftime("%Y%m%d_%H.%M.%S")
    log_file = os.path.join(args.save_path, f"{ctime}_vault_perumi.log")

    # log to file and also to stdout
    logging.basicConfig(
        level=logging.DEBUG,
        format="%(message)s",
        filename=log_file,
        filemode="w"
    )
    logging.getLogger().addHandler(logging.StreamHandler(sys.stdout))

    if args.align_mode != "sr":
        from tools import sniffles1
        try:
            args.sniffles_bin = sniffles1.resolve_sniffles(
                getattr(args, "sniffles_path", None),
                getattr(args, "sv_caller", "sniffles1"),
            )
        except FileNotFoundError as e:
            sys.exit(str(e))

    logging.info(" ".join(sys.argv))
    logging.info("\n-------------- per-UMI variant calling (grouped FASTQ only) -------------")
    logging.info("-vault path:         %s", args.save_path)
    logging.info("-reference sequence: %s", args.refer)
    logging.info("-thread (UMI group): %s", args.thread)
    logging.info("-bash thread:        %s", args.bash_thread)
    logging.info("-align mode:         %s", args.align_mode)
    logging.info("-allele_freq:        %s", args.allele_freq)
    logging.info("-sv_freq:            %s", args.sv_freq)
    logging.info("-group_filter?:      %s", getattr(args, "group_filter", False))
    logging.info("-sv caller:          %s", getattr(args, "sv_caller", "sniffles1"))
    logging.info("-sniffles bin:       %s", getattr(args, "sniffles_bin", ""))
    logging.info("-sniffles min support: %s",
                 getattr(args, "sniffles_min_support", 3))
    logging.info("---------------------------------------------------------------------\n")

    # delayed import to avoid circular imports at module load time
    from vault_main import ProcessUmi, final_clean_up

    # Step1: per-UMI variant calling on grouped FASTQ
    logging.info("%s    Start per-UMI variant calling ..."
                 % datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
    ProcessUmi.run_variant_on_grouped_fastq(args)
    logging.info("%s    Per-UMI variant calling finished"
                 % datetime.now().strftime("%Y-%m-%d %H:%M:%S"))

    # Step2: aggregate per-UMI VCFs and apply final filters
    logging.info("%s    Start aggregating per-UMI VCFs and final filtering ..."
                 % datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
    final_clean_up(args)
    logging.info("%s    Final VCF aggregation finished"
                 % datetime.now().strftime("%Y-%m-%d %H:%M:%S"))

    logging.info("\n------------ Jobs done! ------------\n")
    sys.exit(0)
