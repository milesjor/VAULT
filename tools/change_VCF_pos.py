#!/usr/bin/env python3
# Chongwei 20191027
# bicwei@gmail.com

import argparse
import gzip
import os
import random
import re
import subprocess
import sys


MIN_MAPQ = 20
MIN_QUERY_COVERAGE = 0.50
QC_RANDOM_SEED = 20260603


def get_argparse():
    parser = argparse.ArgumentParser(
        description=(
            "Correct VAULT VCF coordinates by aligning the VAULT analysis "
            "reference to a standard reference genome."
        )
    )
    parser.add_argument(
        "-v",
        "--vcf",
        dest="vcf_file",
        type=validate_file,
        help="path/to/file.vcf or file.vcf.gz",
    )
    parser.add_argument(
        "--vault_result",
        type=str,
        help=(
            "path/to/VAULT_result_folder; convert VCFs directly under "
            "per_umi_process"
        ),
    )
    parser.add_argument(
        "--analysis_ref",
        type=validate_file,
        required=True,
        help="reference FASTA used by VAULT analysis",
    )
    parser.add_argument(
        "--standard_ref",
        type=validate_file,
        required=True,
        help="standard reference FASTA for downstream annotation",
    )
    parser.add_argument(
        "--chr_name",
        type=str,
        help=(
            "force all converted VCF CHROM/CHR2/contig names to this value; "
            "default uses names from --standard_ref"
        ),
    )
    parser.add_argument(
        "--minimap2_path",
        type=str,
        default="minimap2",
        help="minimap2 executable path [minimap2]",
    )
    parser.add_argument(
        "-s",
        "--save_path",
        type=str,
        help=(
            "optional output directory; default is per_umi_process/"
            "pos_correct_vcf for --vault_result, or the input VCF directory "
            "for --vcf"
        ),
    )

    return parser.parse_args()


def validate_file(path):
    if not os.path.exists(path):
        raise argparse.ArgumentTypeError("%s does not exist" % path)
    return path


def validate_args(args):
    has_vcf = args.vcf_file is not None
    has_vault_result = args.vault_result is not None
    if has_vcf == has_vault_result:
        sys.stderr.write("ERROR: provide exactly one of --vcf or --vault_result\n")
        sys.exit(1)

    if has_vault_result:
        per_umi = os.path.join(args.vault_result, "per_umi_process")
        if not os.path.isdir(per_umi):
            sys.stderr.write(
                "ERROR: --vault_result does not contain per_umi_process: %s\n"
                % args.vault_result
            )
            sys.exit(1)


def open_text(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path)


def read_fasta_headers(path):
    lengths = {}
    name = None
    length = 0
    with open(path) as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    lengths[name] = length
                name = line[1:].split()[0]
                length = 0
            else:
                length += len(line)
    if name is not None:
        lengths[name] = length
    return lengths


def read_fasta_selected(path, selected_names=None):
    selected = set(selected_names or [])
    keep_all = selected_names is None
    seqs = {}
    name = None
    chunks = []
    keep = False
    with open(path) as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None and keep:
                    seqs[name] = "".join(chunks).upper()
                name = line[1:].split()[0]
                chunks = []
                keep = keep_all or name in selected
            elif keep:
                chunks.append(line)
    if name is not None and keep:
        seqs[name] = "".join(chunks).upper()
    return seqs


def revcomp(seq):
    table = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return seq.translate(table)[::-1]


def revcomp_allele(allele):
    if allele == "." or allele.startswith("<"):
        return allele
    if any(base not in "ACGTNacgtn" for base in allele):
        return allele
    return revcomp(allele).upper()


def revcomp_alt(alt):
    return ",".join(revcomp_allele(item) for item in alt.split(","))


def parse_cigar(cigar):
    return [
        (int(length), op)
        for length, op in re.findall(r"(\d+)([MIDNSHP=X])", cigar)
    ]


def run_command(cmd, error_label):
    try:
        return subprocess.run(
            cmd,
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
    except FileNotFoundError:
        sys.stderr.write(
            "ERROR: minimap2 not found. Use --minimap2_path or add minimap2 "
            "to PATH.\n"
        )
        sys.exit(1)
    except subprocess.CalledProcessError as exc:
        sys.stderr.write("ERROR: %s failed:\n%s\n" % (error_label, exc.stderr))
        sys.exit(exc.returncode)


def print_minimap2_version(args):
    proc = run_command([args.minimap2_path, "--version"], "minimap2 --version")
    version = (proc.stdout or proc.stderr).strip()
    print("minimap2 version: %s" % (version or "unknown"))


def run_minimap2(args):
    cmd = [
        args.minimap2_path,
        "-cx",
        "asm5",
        "--cs",
        args.standard_ref,
        args.analysis_ref,
    ]
    proc = run_command(cmd, "minimap2 reference alignment")
    paf_lines = [line for line in proc.stdout.splitlines() if line.strip()]
    if not paf_lines:
        sys.stderr.write(
            "ERROR: no alignment found between --analysis_ref and "
            "--standard_ref\n"
        )
        sys.exit(1)
    return paf_lines, proc.stderr


def parse_paf(paf_lines, standard_lengths):
    blocks = []
    for line in paf_lines:
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 12:
            continue

        tags = {}
        for tag in fields[12:]:
            parts = tag.split(":", 2)
            if len(parts) == 3:
                tags[parts[0]] = parts[2]

        block = {
            "analysis": fields[0],
            "analysis_len": int(fields[1]),
            "analysis_start": int(fields[2]) + 1,
            "analysis_end": int(fields[3]),
            "strand": fields[4],
            "standard": fields[5],
            "standard_len": int(fields[6]),
            "standard_start": int(fields[7]) + 1,
            "standard_end": int(fields[8]),
            "match": int(fields[9]),
            "block": int(fields[10]),
            "mapq": int(fields[11]),
            "cigar": tags.get("cg", "%sM" % (int(fields[3]) - int(fields[2]))),
        }
        standard_lengths.setdefault(block["standard"], block["standard_len"])
        blocks.append(block)
    if not blocks:
        sys.stderr.write("ERROR: minimap2 returned no parseable PAF alignments\n")
        sys.exit(1)
    return blocks


def validate_alignment_blocks(blocks, analysis_lengths):
    analysis_names = {block["analysis"] for block in blocks}
    standard_names = {block["standard"] for block in blocks}
    strands = {block["strand"] for block in blocks}

    if len(analysis_names) != 1:
        sys.stderr.write(
            "ERROR: multiple analysis reference sequences were aligned: %s\n"
            % ", ".join(sorted(analysis_names))
        )
        sys.exit(1)
    if len(standard_names) != 1:
        sys.stderr.write(
            "ERROR: alignment hits multiple standard chromosomes: %s. "
            "Please use a more specific --standard_ref or inspect the "
            "reference alignment.\n" % ", ".join(sorted(standard_names))
        )
        sys.exit(1)
    if len(strands) != 1:
        sys.stderr.write(
            "ERROR: alignment contains both forward and reverse-strand hits; "
            "this coordinate conversion is ambiguous.\n"
        )
        sys.exit(1)
    if any(block["mapq"] < MIN_MAPQ for block in blocks):
        sys.stderr.write(
            "ERROR: low-confidence reference alignment detected. Minimum MAPQ "
            "is %s; inspect minimap2 output before coordinate conversion.\n"
            % MIN_MAPQ
        )
        sys.exit(1)

    blocks.sort(key=lambda block: (block["analysis_start"], block["analysis_end"]))
    for prev, curr in zip(blocks, blocks[1:]):
        if curr["analysis_start"] <= prev["analysis_end"]:
            sys.stderr.write(
                "ERROR: overlapping analysis-reference alignment blocks were "
                "detected. This is ambiguous and was not converted.\n"
            )
            sys.exit(1)

    analysis_name = blocks[0]["analysis"]
    analysis_len = analysis_lengths.get(analysis_name, blocks[0]["analysis_len"])
    mapped_query = sum(
        block["analysis_end"] - block["analysis_start"] + 1
        for block in blocks
    )
    query_coverage = mapped_query / analysis_len if analysis_len else 0
    if query_coverage < MIN_QUERY_COVERAGE:
        sys.stderr.write(
            "ERROR: only %.2f%% of --analysis_ref is aligned to --standard_ref; "
            "coordinate conversion is likely unsafe.\n"
            % (query_coverage * 100)
        )
        sys.exit(1)

    if len(blocks) > 2:
        sys.stderr.write(
            "ERROR: more than two reference alignment blocks were detected. "
            "The current VAULT position converter supports one linear block "
            "or two clean split blocks, such as a circular mtDNA reference.\n"
        )
        sys.exit(1)

    if len(blocks) == 2:
        first, second = blocks
        target_overlap = not (
            second["standard_end"] < first["standard_start"]
            or second["standard_start"] > first["standard_end"]
        )
        if target_overlap:
            sys.stderr.write(
                "ERROR: split alignment blocks overlap on the standard "
                "reference. This is ambiguous and was not converted.\n"
            )
            sys.exit(1)

    return blocks, query_coverage


def choose_mapping(existing, candidate):
    if existing is None:
        return candidate
    if candidate["mapq"] > existing["mapq"]:
        return candidate
    if candidate["mapq"] == existing["mapq"] and candidate["block"] > existing["block"]:
        return candidate
    return existing


def build_position_map(blocks):
    mapping = {}
    target_order = []
    for block in blocks:
        if block["standard"] not in target_order:
            target_order.append(block["standard"])

        qpos = block["analysis_start"]
        if block["strand"] == "+":
            tpos = block["standard_start"]
        else:
            tpos = block["standard_end"]

        for length, op in parse_cigar(block["cigar"]):
            if op in ("M", "=", "X"):
                for i in range(length):
                    q_coordinate = qpos + i
                    if block["strand"] == "+":
                        t_coordinate = tpos + i
                    else:
                        t_coordinate = tpos - i
                    candidate = {
                        "chrom": block["standard"],
                        "pos": t_coordinate,
                        "strand": block["strand"],
                        "mapq": block["mapq"],
                        "block": block["block"],
                    }
                    mapping[q_coordinate] = choose_mapping(
                        mapping.get(q_coordinate),
                        candidate,
                    )
                qpos += length
                if block["strand"] == "+":
                    tpos += length
                else:
                    tpos -= length
            elif op in ("I", "S"):
                qpos += length
            elif op in ("D", "N"):
                if block["strand"] == "+":
                    tpos += length
                else:
                    tpos -= length
            elif op in ("H", "P"):
                continue
    return mapping, target_order


def build_liftover_map(args):
    paf_lines, minimap2_stderr = run_minimap2(args)
    analysis_lengths = read_fasta_headers(args.analysis_ref)
    standard_lengths = read_fasta_headers(args.standard_ref)
    blocks = parse_paf(paf_lines, standard_lengths)
    blocks, query_coverage = validate_alignment_blocks(blocks, analysis_lengths)
    mapping, target_order = build_position_map(blocks)
    return {
        "mapping": mapping,
        "standard_lengths": standard_lengths,
        "target_order": target_order,
        "blocks": blocks,
        "query_coverage": query_coverage,
        "minimap2_stderr": minimap2_stderr,
        "forced_chrom": args.chr_name,
    }


def output_chrom(liftover, chrom):
    return liftover["forced_chrom"] or chrom


def output_standard_length(liftover, chrom):
    return liftover["standard_lengths"].get(chrom)


def parse_info(info):
    if info in ("", "."):
        return []
    parsed = []
    for item in info.split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            parsed.append([key, value])
        else:
            parsed.append([item, None])
    return parsed


def format_info(items):
    if not items:
        return "."
    out = []
    for key, value in items:
        if value is None:
            out.append(key)
        else:
            out.append("%s=%s" % (key, value))
    return ";".join(out)


def set_info_value(items, key, value):
    for item in items:
        if item[0] == key:
            item[1] = value
            return
    items.append([key, value])


def append_flag(items, key):
    if not any(item[0] == key for item in items):
        items.append([key, None])


def get_info_value(items, key):
    for item in items:
        if item[0] == key:
            return item[1]
    return None


def split_vcf_name(vcf):
    base = os.path.basename(vcf)
    if base.endswith(".vcf.gz"):
        return base[:-7]
    if base.endswith(".vcf"):
        return base[:-4]
    return base


def write_header_line(outfile, line, liftover, wrote_contigs):
    if line.startswith("##contig="):
        return wrote_contigs

    if line.startswith("#CHROM"):
        if not wrote_contigs:
            seen_names = set()
            for chrom in liftover["target_order"]:
                out_chrom = output_chrom(liftover, chrom)
                if out_chrom in seen_names:
                    continue
                seen_names.add(out_chrom)
                length = output_standard_length(liftover, chrom)
                if length:
                    outfile.write("##contig=<ID=%s,length=%s>\n" % (out_chrom, length))
                else:
                    outfile.write("##contig=<ID=%s>\n" % out_chrom)
            outfile.write(
                "##INFO=<ID=ANALYSIS_CHROM,Number=1,Type=String,"
                "Description=\"Original VAULT analysis reference name\">\n"
            )
            outfile.write(
                "##INFO=<ID=ANALYSIS_POS,Number=1,Type=Integer,"
                "Description=\"Original VAULT analysis reference position\">\n"
            )
            outfile.write(
                "##INFO=<ID=ANALYSIS_END,Number=1,Type=Integer,"
                "Description=\"Original VAULT analysis reference END position\">\n"
            )
            outfile.write(
                "##INFO=<ID=STANDARD_SPAN,Number=1,Type=String,"
                "Description=\"Standard-reference span; circular spans are represented as start-length|1-end\">\n"
            )
            outfile.write(
                "##INFO=<ID=CIRCULAR_SPAN,Number=0,Type=Flag,"
                "Description=\"Variant span crosses the origin of a circular standard reference\">\n"
            )
            wrote_contigs = True
        outfile.write(line)
        return wrote_contigs

    outfile.write(line)
    return wrote_contigs


def convert_record(fields, liftover):
    old_chrom = fields[0]
    old_pos = int(fields[1])
    pos_map = liftover["mapping"].get(old_pos)
    if pos_map is None:
        return None

    info_items = parse_info(fields[7]) if len(fields) > 7 else []
    old_end_value = get_info_value(info_items, "END")
    end_map = None
    if old_end_value and old_end_value.isdigit():
        end_map = liftover["mapping"].get(int(old_end_value))

    new_fields = list(fields)
    new_fields[0] = output_chrom(liftover, pos_map["chrom"])
    new_fields[1] = str(pos_map["pos"])

    if pos_map["strand"] == "-":
        new_fields[3] = revcomp_allele(new_fields[3])
        new_fields[4] = revcomp_alt(new_fields[4])

    set_info_value(info_items, "ANALYSIS_CHROM", old_chrom)
    set_info_value(info_items, "ANALYSIS_POS", str(old_pos))

    if end_map is not None:
        set_info_value(info_items, "ANALYSIS_END", old_end_value)
        set_info_value(info_items, "END", str(end_map["pos"]))
        set_info_value(info_items, "CHR2", output_chrom(liftover, end_map["chrom"]))

        if end_map["chrom"] == pos_map["chrom"]:
            if end_map["pos"] < pos_map["pos"]:
                append_flag(info_items, "CIRCULAR_SPAN")
                span = "%s-%s|1-%s" % (
                    pos_map["pos"],
                    output_standard_length(liftover, pos_map["chrom"]) or "?",
                    end_map["pos"],
                )
            else:
                span = "%s-%s" % (pos_map["pos"], end_map["pos"])
            set_info_value(info_items, "STANDARD_SPAN", span)
        else:
            set_info_value(
                info_items,
                "STANDARD_SPAN",
                "%s:%s-%s:%s" % (
                    output_chrom(liftover, pos_map["chrom"]),
                    pos_map["pos"],
                    output_chrom(liftover, end_map["chrom"]),
                    end_map["pos"],
                ),
            )
    elif old_end_value:
        set_info_value(info_items, "ANALYSIS_END", old_end_value)

    if len(new_fields) > 7:
        new_fields[7] = format_info(info_items)

    sort_key = (
        liftover["target_order"].index(pos_map["chrom"])
        if pos_map["chrom"] in liftover["target_order"] else 999999,
        int(new_fields[1]),
        new_fields[2],
    )
    return sort_key, new_fields


def output_dir_for_vcf(vcf, args):
    if args.save_path:
        return args.save_path
    return os.path.dirname(os.path.abspath(vcf)) or "."


def output_dir_for_vault(args):
    if args.save_path:
        return args.save_path
    return os.path.join(args.vault_result, "per_umi_process", "pos_correct_vcf")


def correct_vcf_by_reference(vcf, save_path, liftover):
    if not os.path.isdir(save_path):
        os.makedirs(save_path, exist_ok=True)

    file_name = split_vcf_name(vcf)
    outvcf = os.path.join(save_path, file_name + "_standard_ref.vcf")
    unmapped_vcf = os.path.join(save_path, file_name + "_standard_ref.unmapped.vcf")

    records = []
    unmapped_records = []
    header_lines = []

    with open_text(vcf) as infile:
        for line in infile:
            if line.startswith("#"):
                header_lines.append(line)
                continue
            fields = line.rstrip("\n").split("\t")
            converted = convert_record(fields, liftover)
            if converted is None:
                unmapped_records.append(line)
            else:
                records.append(converted)

    records.sort(key=lambda item: item[0])

    with open(outvcf, "w") as outfile:
        wrote_contigs = False
        for line in header_lines:
            wrote_contigs = write_header_line(
                outfile,
                line,
                liftover,
                wrote_contigs,
            )
        for _sort_key, fields in records:
            outfile.write("\t".join(fields) + "\n")

    unmapped_output = ""
    if unmapped_records:
        with open(unmapped_vcf, "w") as outfile:
            for line in header_lines:
                outfile.write(line)
            for line in unmapped_records:
                outfile.write(line)
        unmapped_output = unmapped_vcf

    return {
        "input": vcf,
        "output": outvcf,
        "unmapped_output": unmapped_output,
        "records": len(records),
        "unmapped": len(unmapped_records),
    }


def find_vault_result_vcfs(vault_result):
    per_umi = os.path.join(vault_result, "per_umi_process")
    out = []
    for name in sorted(os.listdir(per_umi)):
        path = os.path.join(per_umi, name)
        if not os.path.isfile(path):
            continue
        if name.endswith(".vcf") or name.endswith(".vcf.gz"):
            out.append(path)
    return out


def print_reference_mapping(liftover):
    print("\n-------------- reference alignment schematic --------------")
    print("query coverage: %.2f%%" % (liftover["query_coverage"] * 100))
    for i, block in enumerate(liftover["blocks"], 1):
        direction = "-->" if block["strand"] == "+" else "<--"
        print(
            "block %s: %s:%s-%s %s %s:%s-%s "
            "(strand=%s, mapq=%s, cigar=%s)"
            % (
                i,
                block["analysis"],
                block["analysis_start"],
                block["analysis_end"],
                direction,
                block["standard"],
                block["standard_start"],
                block["standard_end"],
                block["strand"],
                block["mapq"],
                block["cigar"],
            )
        )
    if len(liftover["blocks"]) == 2:
        print(
            "note: two clean split blocks were detected. This is expected for "
            "a circular reference cut at a different origin, such as mtDNA."
        )
    elif len(liftover["blocks"]) == 1:
        print("note: one linear reference alignment block was detected.")


def reference_name_sets(liftover):
    analysis_names = {block["analysis"] for block in liftover["blocks"]}
    standard_names = {block["standard"] for block in liftover["blocks"]}
    return analysis_names, standard_names


def fetch_seq(seqs, chrom, pos, length):
    seq = seqs.get(chrom)
    if seq is None:
        return None
    start = pos - 1
    end = start + length
    if start < 0 or end > len(seq):
        return None
    return seq[start:end]


def parse_info_dict(info):
    out = {}
    for key, value in parse_info(info):
        out[key] = True if value is None else value
    return out


def qc_one_record(fields, liftover, analysis_seqs, standard_seqs):
    chrom, pos, vid, ref, alt = fields[:5]
    pos = int(pos)
    info = parse_info_dict(fields[7] if len(fields) > 7 else ".")
    analysis_pos = info.get("ANALYSIS_POS")
    analysis_chrom = info.get("ANALYSIS_CHROM")
    if not analysis_pos or not analysis_chrom:
        return "id=%s qc=SKIP reason=no_ANALYSIS_POS" % vid

    analysis_pos = int(analysis_pos)
    pos_map = liftover["mapping"].get(analysis_pos)
    if pos_map is None:
        return "id=%s qc=FAIL reason=analysis_pos_not_in_liftover" % vid

    expected_chrom = output_chrom(liftover, pos_map["chrom"])
    pos_ok = expected_chrom == chrom and pos_map["pos"] == pos

    ref_check = "SKIP"
    std_ref = "-"
    expected_ref = "-"
    if re.fullmatch(r"[ACGTNacgtn]+", ref):
        std_ref = fetch_seq(standard_seqs, pos_map["chrom"], pos, len(ref))
        analysis_ref = fetch_seq(
            analysis_seqs,
            str(analysis_chrom),
            analysis_pos,
            len(ref),
        )
        if analysis_ref is not None and pos_map["strand"] == "-":
            analysis_ref = revcomp(analysis_ref)
        expected_ref = analysis_ref or "-"
        if std_ref is None:
            ref_check = "FAIL:no_standard_sequence"
        elif std_ref.upper() == ref.upper():
            ref_check = "PASS"
        else:
            ref_check = "FAIL:standard_REF=%s" % std_ref

    end_check = "NA"
    if "ANALYSIS_END" in info and "END" in info:
        end_map = liftover["mapping"].get(int(info["ANALYSIS_END"]))
        if end_map is None:
            end_check = "FAIL:analysis_end_not_in_liftover"
        elif str(end_map["pos"]) == str(info["END"]):
            end_check = "PASS"
        else:
            end_check = "FAIL:expected_END=%s" % end_map["pos"]

    return (
        "id=%s %s:%s->%s:%s strand=%s pos=%s ref=%s "
        "analysis_ref=%s standard_ref=%s end=%s"
        % (
            vid,
            analysis_chrom,
            analysis_pos,
            chrom,
            pos,
            pos_map["strand"],
            "PASS" if pos_ok else "FAIL",
            ref_check,
            expected_ref,
            std_ref,
            end_check,
        )
    )


def run_qc_for_output(result, liftover, args, max_records=20):
    converted_records = []
    with open(result["output"]) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            converted_records.append(line.rstrip("\n").split("\t"))

    if not converted_records:
        print("QC for %s: no converted records to check" % result["output"])
        return

    sample_size = min(max_records, len(converted_records))
    rng = random.Random(QC_RANDOM_SEED)
    sampled = rng.sample(converted_records, sample_size)

    analysis_names, standard_names = reference_name_sets(liftover)
    analysis_seqs = read_fasta_selected(args.analysis_ref, analysis_names)
    standard_seqs = read_fasta_selected(args.standard_ref, standard_names)

    print("\n-------------- coordinate/DNA QC --------------")
    print("file: %s" % result["output"])
    print("random seed: %s" % QC_RANDOM_SEED)
    print("checked records: %s/%s" % (sample_size, len(converted_records)))
    for fields in sampled:
        print(qc_one_record(fields, liftover, analysis_seqs, standard_seqs))


def position_main(args):
    validate_args(args)
    print("\n-------------- received arguments -------------")
    print("analysis ref:  %s" % args.analysis_ref)
    print("standard ref:  %s" % args.standard_ref)
    print("minimap2 bin:  %s" % args.minimap2_path)
    if args.chr_name:
        print("force chr name: %s" % args.chr_name)
    print_minimap2_version(args)

    liftover = build_liftover_map(args)
    print_reference_mapping(liftover)

    if args.vault_result is not None:
        vcfs = find_vault_result_vcfs(args.vault_result)
        save_path = output_dir_for_vault(args)
    else:
        vcfs = [args.vcf_file]
        save_path = output_dir_for_vcf(args.vcf_file, args)

    if not vcfs:
        sys.stderr.write("ERROR: no VCF files found for conversion\n")
        sys.exit(1)

    print("\n-------------- converting VCF files -----------")
    results = []
    for vcf in vcfs:
        result = correct_vcf_by_reference(vcf, save_path, liftover)
        results.append(result)
        print("input:    %s" % result["input"])
        print("output:   %s" % result["output"])
        print("records:  %s" % result["records"])
        print("unmapped: %s" % result["unmapped"])
        if result["unmapped_output"]:
            print("unmapped output: %s" % result["unmapped_output"])
        print("")

    for result in results:
        run_qc_for_output(result, liftover, args, max_records=20)

    print("\nchange VCF file pos and chr done!\n")
    sys.exit(0)


if __name__ == "__main__":
    position_main(get_argparse())
