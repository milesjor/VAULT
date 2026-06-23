#!/usr/bin/env python3
"""
VDJ detection for MRD (Minimal Residual Disease) using IgBLAST, with optional ONT-style consensus.

Perumi mode workflow (recommended for VAULT):
1) Collect per-UMI grouped FASTQ reads (grouped_reads/perfect_umi/*.fastq[.gz]).
2) Trim the full UMI primer/adaptor sequence using cutadapt:
   - front (5'):  -g ^<ADAPTER_5>
   - back  (3'):  -a <ADAPTER_3>$
   Default:
     ADAPTER_5 = full -u sequence
     ADAPTER_3 = revcomp(full -u sequence)

3) Run igblastn on trimmed reads and parse AIRR TSV (outfmt 19).
4) Call per-UMI clonotype by majority vote across reads in the same group.
5) Summarize dominant clonotypes by UMI count.

Optional: --do_consensus
- Select clonotypes with --consensus_top_k and/or --consensus_top_p, then
  extract reads and generate two DNA consensus:
  (A) amplicon consensus:  reads.amplicon.fastq.gz -> minimap2 -> racon -> (optional medaka)
      reads.amplicon.fastq.gz is extracted from query.raw.fastq.gz (NO adapter trimming),
      so amplicon.consensus.fa contains primer/adaptor sequences (real amplicon).
  (B) VDJ-only consensus: reads.vdj.fastq.gz      -> minimap2 -> racon -> (optional medaka)

Direction policy:
- During consensus extraction, reads are re-oriented to the VDJ direction using IgBLAST AIRR field 'rev_comp'.
  If rev_comp == 'T', we reverse-complement the sequence AND reverse quality (so both amplicon & VDJ are aligned
  to the same VDJ direction).
- After consensus is produced, we run QC:
  (a) check substring containment (vdj in amplicon)
  (b) run igblast on each consensus and confirm V/J/CDR3 match the original clonotype key

NEW:
- Generate vdj_from_amplicon.consensus.fa:
  Run IgBLAST on amplicon.consensus.fa, take V start / J end coords and extract VDJ directly from
  amplicon consensus (ground-truth anchor). This helps debug cases where vdj.consensus.fa is not contained.

Outputs:
- query.raw.fastq.gz
- query.trimmed.fastq.gz
- query.untrimmed.fastq.gz
- query.trimmed.fa
- cutadapt.report.txt
- igblast.airr.tsv
- umi_group_calls.tsv
- clonotype.summary.tsv
- qc.summary.txt
- consensus/<clonotype_id>/reads.amplicon.fastq.gz
- consensus/<clonotype_id>/reads.vdj.fastq.gz
- consensus/<clonotype_id>/amplicon.consensus.fa
- consensus/<clonotype_id>/vdj.consensus.fa
- consensus/<clonotype_id>/vdj_from_amplicon.consensus.fa   (NEW)
- consensus/consensus.summary.tsv
- vdj.log
"""
from __future__ import annotations

import csv
import gzip
import os
import re
import shutil
import subprocess
import sys
import time
from collections import Counter, defaultdict
from dataclasses import dataclass
from typing import Dict, Iterator, List, Tuple, Iterable, Optional

# Prefer VAULT's own cutadapt resolver if available.
try:
    from tools import cutadapt27  # type: ignore
except Exception:  # pragma: no cover
    cutadapt27 = None  # type: ignore


# ---------------------------
# Small utilities
# ---------------------------

def revcomp(seq: str) -> str:
    comp = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return seq.translate(comp)[::-1]


def revcomp_seq_qual(seq: str, qual: str) -> Tuple[str, str]:
    """Reverse-complement a DNA sequence and reverse its FASTQ quality string."""
    return revcomp(seq), qual[::-1]


def open_maybe_gzip(path: str, mode: str = "rt"):
    return gzip.open(path, mode) if path.endswith(".gz") else open(path, mode)


def iter_fastq_records(path: str) -> Iterator[Tuple[str, str, str]]:
    """Yield (read_id, seq, qual) from a FASTQ/FASTQ.GZ file."""
    with open_maybe_gzip(path, "rt") as f:
        while True:
            h = f.readline()
            if not h:
                break
            seq = f.readline().strip()
            _ = f.readline()
            qual = f.readline().strip()
            if not qual:
                break
            rid = h.strip()[1:].split()[0]
            yield rid, seq, qual


def write_fastq_gz(path: str, records: Iterable[Tuple[str, str, str]]) -> int:
    os.makedirs(os.path.dirname(path), exist_ok=True)
    n = 0
    with gzip.open(path, "wt") as out:
        for rid, seq, qual in records:
            out.write(f"@{rid}\n{seq}\n+\n{qual}\n")
            n += 1
    return n


def fastq_to_fasta_gz(in_fastq_gz: str, out_fa: str, group_counter: Optional[Counter] = None) -> int:
    """
    Convert gz-fastq to fasta. If group_counter is provided, count group_id occurrences:
      group_id = read_id.split("|", 1)[0]
    """
    os.makedirs(os.path.dirname(out_fa), exist_ok=True)
    n = 0
    with gzip.open(in_fastq_gz, "rt") as f, open(out_fa, "w") as out:
        while True:
            h = f.readline()
            if not h:
                break
            seq = f.readline().strip()
            _ = f.readline()
            qual = f.readline()
            if not qual:
                break
            rid = h.strip()[1:].split()[0]
            out.write(f">{rid}\n{seq}\n")
            n += 1
            if group_counter is not None:
                gid = rid.split("|", 1)[0]
                group_counter[gid] += 1
    return n


def parse_group_read_count_from_name(group_name: str) -> Optional[int]:
    """VAULT grouped FASTQ base name often starts with read_count, e.g. '32_5end_UMI...'."""
    m = re.match(r"^(\d+)_", group_name)
    if m:
        try:
            return int(m.group(1))
        except ValueError:
            return None
    return None


def read_fasta_one(path: str) -> Tuple[str, str]:
    """Read the first FASTA record (header, seq). Assumes single-record FASTA for consensus."""
    header = ""
    parts: List[str] = []
    with open(path, "r") as f:
        for ln in f:
            ln = ln.rstrip("\n")
            if not ln:
                continue
            if ln.startswith(">"):
                if header:
                    break
                header = ln[1:].strip()
                continue
            parts.append(ln.strip())
    return header, "".join(parts).upper()


def write_fasta_one(path: str, header: str, seq: str, width: int = 80) -> None:
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as f:
        f.write(f">{header}\n")
        for i in range(0, len(seq), max(1, width)):
            f.write(seq[i:i + width] + "\n")


def revcomp_fasta_inplace(path: str, tag: str = "revcomp") -> None:
    h, s = read_fasta_one(path)
    h2 = f"{h}|{tag}" if h else tag
    write_fasta_one(path, h2, revcomp(s))


def collapse_gene(call_field: str, use_allele: bool = False) -> str:
    """
    Convert IgBLAST v_call/j_call/d_call like:
      'IGHV1-69*13,IGHV1-69*15'
    into:
      - allele-level (use_allele=True): 'IGHV1-69*13'
      - gene-level   (use_allele=False): 'IGHV1-69'
    """
    if not call_field:
        return ""
    first = call_field.split(",")[0].strip()
    if use_allele:
        return first
    return first.split("*")[0].strip()


def infer_umi_region(adapter: str) -> Tuple[str, int, int]:
    """
    Infer UMI region from the full adapter sequence:
    - find first 'N' and last 'N'
    - umi_region = adapter[firstN:lastN+1]
    Returns (umi_region, umi_region_len, N_count).
    """
    a = adapter.upper()
    if "N" not in a:
        return ("", 0, 0)
    i = a.find("N")
    j = a.rfind("N")
    umi_region = a[i:j + 1]
    return umi_region, len(umi_region), umi_region.count("N")


def _safe_int(x: str) -> Optional[int]:
    try:
        return int(str(x).strip())
    except Exception:
        return None


# ---------------------------
# Logging
# ---------------------------

class Logger:
    def __init__(self, log_path: str):
        self.log_path = log_path
        os.makedirs(os.path.dirname(log_path), exist_ok=True)

    def _ts(self) -> str:
        return time.strftime("%Y-%m-%d %H:%M:%S")

    def info(self, msg: str) -> None:
        line = f"{self._ts()}\t{msg}"
        print(line, flush=True)
        with open(self.log_path, "a") as f:
            f.write(line + "\n")

    def warn(self, msg: str) -> None:
        self.info(f"[WARN] {msg}")

    def error(self, msg: str) -> None:
        self.info(f"[ERROR] {msg}")


def _is_exe(path: str) -> bool:
    return os.path.isfile(path) and os.access(path, os.X_OK)


def resolve_exe(name_or_path: str) -> Optional[str]:
    if not name_or_path:
        return None
    if os.path.sep in name_or_path or name_or_path.startswith("."):
        return os.path.abspath(name_or_path) if _is_exe(name_or_path) else None
    return shutil.which(name_or_path)


def check_blast_db_prefix(prefix: str) -> bool:
    exts = [".nhr", ".nin", ".nsq"]
    return any(os.path.exists(prefix + e) for e in exts)


# ---------------------------
# IgBLAST
# ---------------------------

def run_igblast_airr(
    igblastn: str,
    query_fa: str,
    out_tsv: str,
    germ_v: str,
    germ_d: str,
    germ_j: str,
    auxiliary_data: str,
    threads: int = 1,
    organism: str = "human",
) -> None:
    cmd = [
        igblastn,
        "-query", query_fa,
        "-out", out_tsv,
        "-outfmt", "19",
        "-ig_seqtype", "Ig",
        "-organism", organism,
        "-domain_system", "imgt",
        "-germline_db_V", germ_v,
        "-germline_db_D", germ_d,
        "-germline_db_J", germ_j,
        "-auxiliary_data", auxiliary_data,
        "-num_threads", str(threads),
        "-strand", "both",
    ]
    subprocess.run(cmd, check=True)


def parse_airr_tsv(path: str) -> Iterator[Dict[str, str]]:
    with open(path, "r", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            yield row


def _as_float(x: str) -> float:
    try:
        return float(x)
    except Exception:
        return 0.0


def _pass_identity(val: str, min_identity: float) -> bool:
    if min_identity <= 0:
        return True
    thr = float(min_identity)
    if 0 < thr <= 1.0:
        thr *= 100.0
    v = _as_float(val)
    if v <= 1.0:
        v *= 100.0
    return v >= thr


def _row_pass_productive(row: Dict[str, str]) -> bool:
    prod = (row.get("productive") or "").strip().upper()
    in_frame = (row.get("vj_in_frame") or "").strip().upper()
    stop = (row.get("stop_codon") or "").strip().upper()
    return prod == "T" and in_frame == "T" and stop == "F"


def _cdr3_tag_from_row(row: Dict[str, str], cdr3_mode: str) -> str:
    cdr3aa = (row.get("cdr3_aa") or "").strip()
    cdr3nt = (row.get("cdr3") or "").strip()
    if cdr3_mode == "aa":
        return cdr3aa or cdr3nt
    return cdr3nt


def _row_pass_basic_validity(
    row: Dict[str, str],
    include_d: bool,
    cdr3_mode: str,
    min_v_identity: float,
    min_j_identity: float,
) -> bool:
    v_call = (row.get("v_call") or "").strip()
    j_call = (row.get("j_call") or "").strip()
    d_call = (row.get("d_call") or "").strip()

    if not (v_call and j_call):
        return False

    cdr3_tag = _cdr3_tag_from_row(row, cdr3_mode=cdr3_mode)
    if not cdr3_tag:
        return False

    if include_d and not d_call:
        return False

    if not _pass_identity((row.get("v_identity") or "").strip(), min_v_identity):
        return False
    if not _pass_identity((row.get("j_identity") or "").strip(), min_j_identity):
        return False

    return True


def build_read_hit_maps(
    airr_tsv: str,
    require_productive: bool,
    include_d: bool,
    use_allele: bool,
    cdr3_mode: str,
    min_v_identity: float,
    min_j_identity: float,
) -> Tuple[Dict[str, tuple], Dict[str, Tuple[int, int]], Dict[str, bool]]:
    """
    Build read-level maps:
      - read_to_key[sequence_id] = clonotype key
      - read_to_coords[sequence_id] = (v_sequence_start, j_sequence_end) for VDJ-only cropping
      - read_to_revcomp[sequence_id] = True if AIRR 'rev_comp' == 'T'
    """
    read_to_key: Dict[str, tuple] = {}
    read_to_coords: Dict[str, Tuple[int, int]] = {}
    read_to_revcomp: Dict[str, bool] = {}

    for row in parse_airr_tsv(airr_tsv):
        rid = (row.get("sequence_id") or "").strip()
        if not rid:
            continue

        if not _row_pass_basic_validity(
            row=row,
            include_d=include_d,
            cdr3_mode=cdr3_mode,
            min_v_identity=min_v_identity,
            min_j_identity=min_j_identity,
        ):
            continue

        if require_productive and not _row_pass_productive(row):
            continue

        v_gene = collapse_gene((row.get("v_call") or "").strip(), use_allele=use_allele)
        j_gene = collapse_gene((row.get("j_call") or "").strip(), use_allele=use_allele)
        d_gene = collapse_gene((row.get("d_call") or "").strip(), use_allele=use_allele) if include_d else ""
        cdr3_tag = _cdr3_tag_from_row(row, cdr3_mode=cdr3_mode)

        key = (v_gene, d_gene, j_gene, cdr3_tag) if include_d else (v_gene, j_gene, cdr3_tag)
        read_to_key[rid] = key

        v_start = _safe_int(row.get("v_sequence_start") or "")
        j_end = _safe_int(row.get("j_sequence_end") or "")
        if v_start is not None and j_end is not None and v_start > 0 and j_end > 0 and j_end >= v_start:
            read_to_coords[rid] = (v_start, j_end)

        rc = (row.get("rev_comp") or "").strip().upper()
        read_to_revcomp[rid] = (rc == "T")

    return read_to_key, read_to_coords, read_to_revcomp


def igblast_first_row(
    igblastn_exe: str,
    query_fa: str,
    out_tsv: str,
    germ_v: str,
    germ_d: str,
    germ_j: str,
    auxiliary_data: str,
    organism: str,
    threads: int = 1,
) -> Optional[Dict[str, str]]:
    """Run igblast on query_fa and return the first AIRR row (or None)."""
    run_igblast_airr(
        igblastn=igblastn_exe,
        query_fa=query_fa,
        out_tsv=out_tsv,
        germ_v=germ_v,
        germ_d=germ_d,
        germ_j=germ_j,
        auxiliary_data=auxiliary_data,
        threads=threads,
        organism=organism,
    )
    for row in parse_airr_tsv(out_tsv):
        return row
    return None


def normalize_consensus_orientation_by_igblast(
    cons_fa: str,
    label: str,
    outdir: str,
    igblastn_exe: str,
    germ_v: str,
    germ_d: str,
    germ_j: str,
    auxiliary_data: str,
    organism: str,
    logger: Logger,
) -> Tuple[Optional[Dict[str, str]], bool]:
    """
    Ensure consensus fasta is in VDJ direction.
    If igblast says rev_comp == 'T', reverse-complement the fasta in-place and re-run igblast.
    Returns (row, flipped).
    """
    tsv = os.path.join(outdir, f"qc.igblast.{label}.airr.tsv")
    try:
        row = igblast_first_row(
            igblastn_exe=igblastn_exe,
            query_fa=cons_fa,
            out_tsv=tsv,
            germ_v=germ_v,
            germ_d=germ_d,
            germ_j=germ_j,
            auxiliary_data=auxiliary_data,
            organism=organism,
            threads=1,
        )
    except Exception as e:
        logger.warn(f"[QC] igblast failed on {label} consensus ({cons_fa}): {e}")
        return None, False

    if not row:
        logger.warn(f"[QC] igblast returned no rows for {label} consensus ({cons_fa})")
        return None, False

    rc = (row.get("rev_comp") or "").strip().upper()
    if rc == "T":
        revcomp_fasta_inplace(cons_fa, tag="revcomp_by_igblast")
        logger.warn(f"[QC] {label} consensus rev_comp=T => reverse-complemented in-place to enforce VDJ direction.")
        try:
            row2 = igblast_first_row(
                igblastn_exe=igblastn_exe,
                query_fa=cons_fa,
                out_tsv=tsv,
                germ_v=germ_v,
                germ_d=germ_d,
                germ_j=germ_j,
                auxiliary_data=auxiliary_data,
                organism=organism,
                threads=1,
            )
            return row2, True
        except Exception as e:
            logger.warn(f"[QC] igblast re-run failed after flipping {label} consensus: {e}")
            return row, True

    return row, False


# ---------------------------
# Group call
# ---------------------------

@dataclass
class GroupCall:
    group_id: str
    read_count: int
    productive_reads: int
    called: bool
    v_gene: str = ""
    d_gene: str = ""
    j_gene: str = ""
    cdr3_tag: str = ""
    support: int = 0
    support_frac: float = 0.0
    note: str = ""


def choose_group_call(
    group_id: str,
    read_count: int,
    per_read_rows: List[Dict[str, str]],
    require_productive: bool = True,
    min_support: int = 1,
    include_d: bool = False,
    use_allele: bool = False,
    cdr3_mode: str = "aa",
    min_v_identity: float = 0.0,
    min_j_identity: float = 0.0,
) -> GroupCall:
    valid_rows: List[Dict[str, str]] = []
    productive_rows: List[Dict[str, str]] = []

    for row in per_read_rows:
        if not _row_pass_basic_validity(
            row=row,
            include_d=include_d,
            cdr3_mode=cdr3_mode,
            min_v_identity=min_v_identity,
            min_j_identity=min_j_identity,
        ):
            continue

        valid_rows.append(row)
        if _row_pass_productive(row):
            productive_rows.append(row)

    productive_n = len(productive_rows)
    considered = productive_rows if require_productive else valid_rows

    if not considered:
        return GroupCall(group_id, read_count, productive_n, False, note="NO_VALID_HITS")

    counter = Counter()
    score_sum: Dict[tuple, float] = defaultdict(float)

    for row in considered:
        v_gene = collapse_gene((row.get("v_call") or "").strip(), use_allele=use_allele)
        d_gene = collapse_gene((row.get("d_call") or "").strip(), use_allele=use_allele) if include_d else ""
        j_gene = collapse_gene((row.get("j_call") or "").strip(), use_allele=use_allele)
        cdr3_tag = _cdr3_tag_from_row(row, cdr3_mode=cdr3_mode)

        key = (v_gene, d_gene, j_gene, cdr3_tag) if include_d else (v_gene, j_gene, cdr3_tag)
        counter[key] += 1
        score_sum[key] += _as_float(row.get("v_score") or "") + _as_float(row.get("j_score") or "")

    best_key = None
    best_support = -1
    best_score = -1.0
    for key, n in counter.items():
        sc = score_sum.get(key, 0.0)
        if n > best_support or (n == best_support and sc > best_score):
            best_key, best_support, best_score = key, n, sc

    if best_support < min_support:
        return GroupCall(group_id, read_count, productive_n, False, note=f"LOW_SUPPORT<{min_support}")

    if include_d:
        v_gene, d_gene, j_gene, cdr3_tag = best_key  # type: ignore[misc]
    else:
        v_gene, j_gene, cdr3_tag = best_key  # type: ignore[misc]
        d_gene = ""

    support_frac = best_support / max(1, len(considered))
    note = "PASS" if require_productive else "PASS_NONPRODUCTIVE_ALLOWED"
    return GroupCall(
        group_id=group_id,
        read_count=read_count,
        productive_reads=productive_n,
        called=True,
        v_gene=v_gene,
        d_gene=d_gene,
        j_gene=j_gene,
        cdr3_tag=cdr3_tag,
        support=best_support,
        support_frac=support_frac,
        note=note,
    )


def write_tsv(path: str, header: List[str], rows: Iterable[List[str]]) -> None:
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", newline="") as f:
        f.write("\t".join(header) + "\n")
        for r in rows:
            f.write("\t".join(r) + "\n")


# ---------------------------
# cutadapt trimming
# ---------------------------

CutadaptCmd = List[str]


def resolve_cutadapt_cmd(args, logger: Logger) -> CutadaptCmd:
    def as_cmd(path: str) -> CutadaptCmd:
        if _is_exe(path):
            return [os.path.abspath(path)]
        if os.path.isfile(path):
            return [sys.executable, os.path.abspath(path)]
        raise RuntimeError(f"cutadapt not runnable: {path}")

    if hasattr(args, "cutadapt_path") and args.cutadapt_path:
        return as_cmd(args.cutadapt_path)

    if cutadapt27 is not None:
        try:
            p = cutadapt27.resolve_cutadapt(getattr(args, "cutadapt_path", None))
            if p:
                try:
                    return as_cmd(p)
                except Exception:
                    q = resolve_exe(p)
                    if q:
                        return [q]
        except Exception as e:
            logger.warn(f"VAULT cutadapt resolver failed, falling back to PATH. Detail: {e}")

    p = resolve_exe("cutadapt")
    if not p:
        raise RuntimeError("cutadapt not found. Provide --cutadapt_path or install cutadapt.")
    return [p]


def cutadapt_supports_revcomp(cutadapt_cmd: CutadaptCmd, logger: Logger) -> bool:
    """
    Detect whether this cutadapt supports --revcomp.
    """
    try:
        out = subprocess.check_output(cutadapt_cmd + ["--help"], stderr=subprocess.STDOUT, text=True)
        return "--revcomp" in out
    except Exception as e:
        logger.warn(f"cutadapt --help detection failed ({e}). Assume --revcomp NOT supported.")
        return False


def run_cutadapt_full_primer(
    cutadapt_cmd: CutadaptCmd,
    in_fastq_gz: str,
    out_trimmed_fastq_gz: str,
    out_untrimmed_fastq_gz: str,
    report_path: str,
    adapter_5: str,
    adapter_3: str,
    threads: int = 4,
    error_rate: float = 0.10,
    enable_revcomp: bool = True,
) -> Dict[str, int]:
    fwd = adapter_5.upper()
    a3 = adapter_3.upper()

    base_opts = [
        "-j", str(max(1, threads)),
        "-e", str(error_rate),
        "-g", f"^{fwd}",
        "-a", f"{a3}$",
        "-o", out_trimmed_fastq_gz,
        "--untrimmed-output", out_untrimmed_fastq_gz,
        in_fastq_gz,
    ]
    cmd = cutadapt_cmd + (["--revcomp"] if enable_revcomp else []) + base_opts

    os.makedirs(os.path.dirname(report_path), exist_ok=True)
    with open(report_path, "w") as rep:
        subprocess.run(cmd, check=True, stdout=rep, stderr=subprocess.STDOUT)

    metrics = {"cutadapt_kept": 0, "cutadapt_untrimmed": 0}
    try:
        txt = open(report_path, "r").read()
        m = re.search(r"Reads written \(passing filters\):\s+([\d,]+)", txt)
        if m:
            metrics["cutadapt_kept"] = int(m.group(1).replace(",", ""))
        m2 = re.search(r"Reads discarded as untrimmed:\s+([\d,]+)", txt)
        if m2:
            metrics["cutadapt_untrimmed"] = int(m2.group(1).replace(",", ""))
    except Exception:
        pass

    return metrics


# ---------------------------
# Consensus utilities
# ---------------------------

def build_longest_read_draft(reads_fastq_gz: str, draft_fa: str) -> int:
    longest = ""
    longest_id = ""
    for rid, seq, _qual in iter_fastq_records(reads_fastq_gz):
        if len(seq) > len(longest):
            longest = seq
            longest_id = rid
    os.makedirs(os.path.dirname(draft_fa), exist_ok=True)
    with open(draft_fa, "w") as f:
        f.write(f">draft_from_{longest_id}\n{longest}\n")
    return len(longest)


def consensus_minimap2_racon(
    reads_fastq_gz: str,
    out_dir: str,
    threads: int,
    racon_rounds: int,
    minimap2_exe: str,
    racon_exe: str,
    logger: Logger,
) -> str:
    os.makedirs(out_dir, exist_ok=True)
    draft = os.path.join(out_dir, "draft.0.fa")
    build_longest_read_draft(reads_fastq_gz, draft)

    cur = draft
    for i in range(1, racon_rounds + 1):
        paf = os.path.join(out_dir, f"aln.{i}.paf")
        new = os.path.join(out_dir, f"draft.racon{i}.fa")

        with open(paf, "w") as f:
            subprocess.run(
                [minimap2_exe, "-x", "map-ont", "-t", str(max(1, threads)), cur, reads_fastq_gz],
                check=True,
                stdout=f,
                stderr=subprocess.DEVNULL,
            )

        with open(new, "w") as f:
            subprocess.run(
                [racon_exe, "-t", str(max(1, threads)), reads_fastq_gz, paf, cur],
                check=True,
                stdout=f,
                stderr=subprocess.DEVNULL,
            )

        cur = new
        logger.info(f"racon round {i} done: {cur}")

    return cur


def consensus_medaka(
    reads_fastq_gz: str,
    draft_fa: str,
    out_dir: str,
    threads: int,
    medaka_consensus_exe: str,
    medaka_model: str,
    logger: Logger,
) -> str:
    os.makedirs(out_dir, exist_ok=True)
    cmd = [
        medaka_consensus_exe,
        "-i", reads_fastq_gz,
        "-d", draft_fa,
        "-o", out_dir,
        "-t", str(max(1, threads)),
        "-m", medaka_model,
    ]
    logger.info("RUN: " + " ".join(cmd))
    subprocess.run(cmd, check=True)
    out_fa = os.path.join(out_dir, "consensus.fasta")
    if not os.path.exists(out_fa):
        raise RuntimeError("medaka_consensus did not produce consensus.fasta")
    return out_fa


def crop_vdj(seq: str, qual: str, v_start: int, j_end: int) -> Optional[Tuple[str, str]]:
    """
    Crop read to VDJ-only region using 1-based inclusive coords:
      [v_sequence_start .. j_sequence_end]
    """
    if v_start <= 0 or j_end <= 0 or j_end < v_start:
        return None
    if len(seq) != len(qual):
        return None
    s = max(0, v_start - 1)
    e = min(len(seq), j_end)
    if e <= s:
        return None
    return seq[s:e], qual[s:e]


def _run_consensus_pair(
    reads_fq: str,
    outdir: str,
    label: str,
    minimap2_exe: str,
    racon_exe: str,
    threads: int,
    racon_rounds: int,
    logger: Logger,
    medaka_enabled: bool,
    medaka_exe: Optional[str],
    medaka_model: str,
) -> str:
    racon_out = consensus_minimap2_racon(
        reads_fastq_gz=reads_fq,
        out_dir=os.path.join(outdir, f"racon.{label}"),
        threads=threads,
        racon_rounds=racon_rounds,
        minimap2_exe=minimap2_exe,
        racon_exe=racon_exe,
        logger=logger,
    )
    final = racon_out
    if medaka_enabled:
        assert medaka_exe is not None
        final = consensus_medaka(
            reads_fastq_gz=reads_fq,
            draft_fa=racon_out,
            out_dir=os.path.join(outdir, f"medaka.{label}"),
            threads=threads,
            medaka_consensus_exe=medaka_exe,
            medaka_model=medaka_model,
            logger=logger,
        )
    out_fa = os.path.join(outdir, f"{label}.consensus.fa")
    shutil.copyfile(final, out_fa)
    return out_fa


def minimap2_best_strand(target_fa: str, query_fa: str, minimap2_exe: str) -> Optional[str]:
    """Return '+' or '-' from the first PAF line of minimap2 mapping (or None if no alignment)."""
    try:
        out = subprocess.check_output(
            [minimap2_exe, "-x", "asm5", target_fa, query_fa],
            stderr=subprocess.DEVNULL,
            text=True,
        )
    except Exception:
        return None
    out = out.strip()
    if not out:
        return None
    first = out.splitlines()[0]
    fields = first.split("\t")
    if len(fields) < 5:
        return None
    return fields[4].strip()


def extract_vdj_from_amplicon_consensus(
    amplicon_cons_fa: str,
    amp_row: Dict[str, str],
    out_fa: str,
    logger: Logger,
) -> Optional[str]:
    """
    NEW:
    Using IgBLAST AIRR row on amplicon consensus, extract [v_sequence_start..j_sequence_end]
    from the (already oriented) amplicon consensus sequence to get a ground-truth VDJ.
    """
    v_start = _safe_int(amp_row.get("v_sequence_start") or "")
    j_end = _safe_int(amp_row.get("j_sequence_end") or "")
    if v_start is None or j_end is None or v_start <= 0 or j_end <= 0 or j_end < v_start:
        logger.warn(f"[QC] vdj_from_amplicon: invalid coords v_start={v_start} j_end={j_end}")
        return None

    h, seq = read_fasta_one(amplicon_cons_fa)
    if not seq:
        logger.warn("[QC] vdj_from_amplicon: empty amplicon consensus sequence")
        return None

    s = max(0, v_start - 1)
    e = min(len(seq), j_end)
    if e <= s:
        logger.warn(f"[QC] vdj_from_amplicon: slice empty (len={len(seq)} s={s} e={e})")
        return None

    vdj_seq = seq[s:e]
    hdr = f"vdj_from_amplicon|{h}|v_sequence_start={v_start}|j_sequence_end={j_end}"
    write_fasta_one(out_fa, hdr, vdj_seq)
    return out_fa


# ---------------------------
# Pipeline modules
# ---------------------------

def _setup_outdir_and_logger(args) -> Tuple[str, Logger]:
    if getattr(args, "out_dir", None) is None:
        args.out_dir = os.path.join(args.save_path, "vdj")
    os.makedirs(args.out_dir, exist_ok=True)
    logger = Logger(os.path.join(args.out_dir, "vdj.log"))
    return args.out_dir, logger


def _resolve_adapter_params(args, logger: Logger) -> Tuple[str, str, str]:
    umi_adapter = getattr(args, "umi_adapter", None)
    if not umi_adapter:
        raise ValueError("vault vdj requires -u/--umi_adapter (full primer with N-pattern).")
    umi_adapter = umi_adapter.strip().upper()

    adapter_5 = (getattr(args, "adapter_5", None) or umi_adapter).strip().upper()
    adapter_3 = (getattr(args, "adapter_3", None) or revcomp(umi_adapter)).strip().upper()

    umi_region, umi_region_len, umi_n_count = infer_umi_region(umi_adapter)
    logger.info(f"UMI adapter (-u): {umi_adapter}")
    logger.info(f"Inferred UMI region: {umi_region} (len={umi_region_len}, N_count={umi_n_count})")
    logger.info(f"Adapter_5 used for cutadapt -g: {adapter_5}")
    logger.info(f"Adapter_3 used for cutadapt -a: {adapter_3}")
    return umi_adapter, adapter_5, adapter_3


def _resolve_tools_and_validate(args, logger: Logger) -> Tuple[str, CutadaptCmd]:
    igblastn_exe = resolve_exe(args.igblastn) if getattr(args, "igblastn", None) else None
    if not igblastn_exe:
        raise RuntimeError(f"igblastn not found or not executable: {args.igblastn}")

    if not check_blast_db_prefix(args.germV):
        logger.warn(f"V germline DB prefix seems missing BLAST files: {args.germV}(.nhr/.nin/.nsq)")
    if not check_blast_db_prefix(args.germD):
        logger.warn(f"D germline DB prefix seems missing BLAST files: {args.germD}(.nhr/.nin/.nsq)")
    if not check_blast_db_prefix(args.germJ):
        logger.warn(f"J germline DB prefix seems missing BLAST files: {args.germJ}(.nhr/.nin/.nsq)")
    if not os.path.exists(args.auxiliary_data):
        raise RuntimeError(f"IgBLAST auxiliary_data not found: {args.auxiliary_data}")

    cutadapt_cmd = resolve_cutadapt_cmd(args, logger)
    logger.info(f"Using cutadapt command: {' '.join(cutadapt_cmd)}")
    return igblastn_exe, cutadapt_cmd


def _determine_mode(args, logger: Logger) -> str:
    mode = args.mode
    if mode == "auto":
        default_group_dir = os.path.join(args.save_path, "grouped_reads", "perfect_umi")
        mode = "perumi" if os.path.isdir(default_group_dir) else "raw"
    if mode not in ("raw", "perumi"):
        raise ValueError(f"Invalid mode: {mode}")
    logger.info(f"Mode: {mode}")
    return mode


def _build_query_raw_fastq(
    args,
    mode: str,
    out_dir: str,
    logger: Logger,
) -> Tuple[str, Dict[str, int], int]:
    query_raw_fq = os.path.join(out_dir, "query.raw.fastq.gz")
    group_read_counts_raw: Dict[str, int] = {}
    total_reads_written = 0

    def raw_records() -> Iterator[Tuple[str, str, str]]:
        nonlocal total_reads_written
        if mode == "raw":
            if not getattr(args, "fastq", None):
                raise ValueError("In raw mode, --fastq is required.")
            for rid, seq, qual in iter_fastq_records(args.fastq):
                total_reads_written += 1
                yield rid, seq, qual
        else:
            group_dir = args.grouped_dir or os.path.join(args.save_path, "grouped_reads", "perfect_umi")
            if not os.path.isdir(group_dir):
                raise ValueError(f"Grouped FASTQ directory not found: {group_dir}")

            for fn in sorted(os.listdir(group_dir)):
                if not (fn.endswith(".fastq") or fn.endswith(".fastq.gz") or fn.endswith(".fq") or fn.endswith(".fq.gz")):
                    continue
                if fn.endswith(".read2.fastq.gz") or fn.endswith(".read2.fastq"):
                    continue

                group_id = fn
                for suffix in (".fastq.gz", ".fastq", ".fq.gz", ".fq"):
                    if group_id.endswith(suffix):
                        group_id = group_id[: -len(suffix)]
                        break

                fq_path = os.path.join(group_dir, fn)
                read_n = 0
                for i, (_rid0, seq, qual) in enumerate(iter_fastq_records(fq_path), start=1):
                    if args.max_reads_per_group > 0 and i > args.max_reads_per_group:
                        break
                    seq_id = f"{group_id}|r{i:06d}"
                    total_reads_written += 1
                    read_n += 1
                    yield seq_id, seq, qual

                if read_n == 0:
                    prefix_n = parse_group_read_count_from_name(group_id)
                    read_n = prefix_n or 0
                group_read_counts_raw[group_id] = read_n

    total_reads_written = write_fastq_gz(query_raw_fq, raw_records())
    logger.info(f"Wrote merged query FASTQ: {query_raw_fq} (reads={total_reads_written})")
    return query_raw_fq, group_read_counts_raw, total_reads_written


def _run_cutadapt_and_make_query_fa(
    args,
    mode: str,
    query_raw_fq: str,
    adapter_5: str,
    adapter_3: str,
    cutadapt_cmd: CutadaptCmd,
    out_dir: str,
    logger: Logger,
) -> Tuple[str, str, str, Optional[Counter], int, int, bool, float]:
    query_trim_fq = os.path.join(out_dir, "query.trimmed.fastq.gz")
    query_untrim_fq = os.path.join(out_dir, "query.untrimmed.fastq.gz")
    query_fa = os.path.join(out_dir, "query.trimmed.fa")
    cutadapt_report = os.path.join(out_dir, "cutadapt.report.txt")

    error_rate = float(getattr(args, "cutadapt_error_rate", 0.10))
    enable_revcomp = bool(getattr(args, "cutadapt_revcomp", True))
    if enable_revcomp:
        if cutadapt_supports_revcomp(cutadapt_cmd, logger):
            logger.info("cutadapt supports --revcomp: enabled.")
        else:
            logger.warn("cutadapt does NOT support --revcomp => auto-disable revcomp trimming.")
            enable_revcomp = False

    logger.info(
        f"Trim adapters: -g ^{adapter_5} ; -a {adapter_3}$ ; "
        f"revcomp_mode={int(enable_revcomp)} ; -e {error_rate}"
    )

    run_cutadapt_full_primer(
        cutadapt_cmd=cutadapt_cmd,
        in_fastq_gz=query_raw_fq,
        out_trimmed_fastq_gz=query_trim_fq,
        out_untrimmed_fastq_gz=query_untrim_fq,
        report_path=cutadapt_report,
        adapter_5=adapter_5,
        adapter_3=adapter_3,
        threads=int(getattr(args, "threads", 4)),
        error_rate=error_rate,
        enable_revcomp=enable_revcomp,
    )

    group_trimmed_counts: Optional[Counter] = Counter() if mode == "perumi" else None
    trimmed_n = fastq_to_fasta_gz(query_trim_fq, query_fa, group_counter=group_trimmed_counts)
    untrimmed_n = sum(1 for _ in iter_fastq_records(query_untrim_fq))
    logger.info(f"cutadapt done. trimmed_reads={trimmed_n}, untrimmed_reads={untrimmed_n}, report={cutadapt_report}")
    logger.info(f"Built IgBLAST query FASTA: {query_fa} (records={trimmed_n})")

    return query_trim_fq, query_untrim_fq, query_fa, group_trimmed_counts, trimmed_n, untrimmed_n, enable_revcomp, error_rate


def _run_igblast_main(args, igblastn_exe: str, query_fa: str, out_dir: str, logger: Logger) -> str:
    igblast_tsv = os.path.join(out_dir, "igblast.airr.tsv")
    logger.info("Running IgBLAST (AIRR outfmt 19) ...")
    run_igblast_airr(
        igblastn=igblastn_exe,
        query_fa=query_fa,
        out_tsv=igblast_tsv,
        germ_v=args.germV,
        germ_d=args.germD,
        germ_j=args.germJ,
        auxiliary_data=args.auxiliary_data,
        threads=int(getattr(args, "threads", 8)),
        organism=args.organism,
    )
    logger.info(f"IgBLAST finished: {igblast_tsv}")
    return igblast_tsv


def _raw_mode_report(
    args,
    igblast_tsv: str,
    out_dir: str,
    logger: Logger,
    include_d: bool,
    use_allele: bool,
    cdr3_mode: str,
    min_v_identity: float,
    min_j_identity: float,
    require_productive: bool,
    trimmed_n: int,
    untrimmed_n: int,
    error_rate: float,
    enable_revcomp: bool,
) -> None:
    clonotype_tsv = os.path.join(out_dir, "clonotype.summary.tsv")
    qc_txt = os.path.join(out_dir, "qc.summary.txt")

    top = Counter()
    total = 0
    kept = 0
    for row in parse_airr_tsv(igblast_tsv):
        total += 1
        if not _row_pass_basic_validity(row, include_d, cdr3_mode, min_v_identity, min_j_identity):
            continue
        if require_productive and not _row_pass_productive(row):
            continue

        v_gene = collapse_gene((row.get("v_call") or "").strip(), use_allele=use_allele)
        j_gene = collapse_gene((row.get("j_call") or "").strip(), use_allele=use_allele)
        d_gene = collapse_gene((row.get("d_call") or "").strip(), use_allele=use_allele) if include_d else ""
        cdr3_tag = _cdr3_tag_from_row(row, cdr3_mode=cdr3_mode)

        kept += 1
        if include_d:
            top[(v_gene, d_gene, j_gene, cdr3_tag)] += 1
        else:
            top[(v_gene, j_gene, cdr3_tag)] += 1

    denom = max(1, sum(top.values()))
    rows = []
    for rank, (key, n) in enumerate(top.most_common(args.top_k), start=1):
        if include_d:
            v, d, j, c = key
            rows.append([str(rank), str(n), f"{n/denom:.6f}", v, d, j, c])
        else:
            v, j, c = key
            rows.append([str(rank), str(n), f"{n/denom:.6f}", v, j, c])

    header = ["rank", "read_count", "read_fraction", "V_gene"]
    if include_d:
        header += ["D_gene"]
    header += ["J_gene", "CDR3_tag"]
    write_tsv(clonotype_tsv, header, rows)

    with open(qc_txt, "w") as f:
        f.write("mode\traw\n")
        f.write(f"total_records\t{total}\n")
        f.write(f"kept_records\t{kept}\n")
        f.write(f"top_k\t{args.top_k}\n")
        f.write(f"cutadapt_trimmed_reads\t{trimmed_n}\n")
        f.write(f"cutadapt_untrimmed_reads\t{untrimmed_n}\n")
        f.write(f"cutadapt_error_rate\t{error_rate}\n")
        f.write(f"cutadapt_revcomp\t{int(enable_revcomp)}\n")
        f.write(f"require_productive\t{int(require_productive)}\n")

    logger.info("RAW mode done.")


def _perumi_calls_and_summary(
    args,
    igblast_tsv: str,
    group_read_counts_raw: Dict[str, int],
    group_trimmed_counts: Counter,
    include_d: bool,
    use_allele: bool,
    cdr3_mode: str,
    min_v_identity: float,
    min_j_identity: float,
    require_productive: bool,
    logger: Logger,
) -> Tuple[List[GroupCall], List[Tuple[int, int, tuple]], Dict[str, tuple], int]:
    logger.info("Aggregating IgBLAST AIRR rows by UMI group ...")
    group_to_rows: Dict[str, List[Dict[str, str]]] = defaultdict(list)
    for row in parse_airr_tsv(igblast_tsv):
        sid = (row.get("sequence_id") or "").strip()
        if not sid:
            continue
        group_id = sid.split("|", 1)[0]
        group_to_rows[group_id].append(row)

    group_read_counts_eff: Dict[str, int] = {}
    for gid, n in group_trimmed_counts.items():
        group_read_counts_eff[gid] = int(n)
    for gid in group_read_counts_raw.keys():
        group_read_counts_eff.setdefault(gid, 0)

    group_calls: List[GroupCall] = []
    for group_id, eff_count in group_read_counts_eff.items():
        if eff_count < args.min_reads_per_umi:
            group_calls.append(GroupCall(group_id, eff_count, 0, False, note=f"READS<{args.min_reads_per_umi}"))
            continue
        rows = group_to_rows.get(group_id, [])
        call = choose_group_call(
            group_id=group_id,
            read_count=eff_count,
            per_read_rows=rows,
            require_productive=require_productive,
            min_support=int(args.min_support),
            include_d=include_d,
            use_allele=use_allele,
            cdr3_mode=cdr3_mode,
            min_v_identity=min_v_identity,
            min_j_identity=min_j_identity,
        )
        group_calls.append(call)

    clonotypes = defaultdict(lambda: {"umi": 0, "reads": 0})
    group_to_key: Dict[str, tuple] = {}
    for c in group_calls:
        if not c.called:
            continue
        key = (c.v_gene, c.d_gene, c.j_gene, c.cdr3_tag) if include_d else (c.v_gene, c.j_gene, c.cdr3_tag)
        clonotypes[key]["umi"] += 1
        clonotypes[key]["reads"] += c.read_count
        group_to_key[c.group_id] = key

    items = sorted([(v["umi"], v["reads"], k) for k, v in clonotypes.items()], reverse=True)
    return group_calls, items, group_to_key, len(group_read_counts_eff)


def _write_perumi_outputs(
    args,
    out_dir: str,
    group_calls: List[GroupCall],
    items: List[Tuple[int, int, tuple]],
    total_groups: int,
    total_reads_written: int,
    trimmed_n: int,
    untrimmed_n: int,
    error_rate: float,
    enable_revcomp: bool,
    include_d: bool,
    use_allele: bool,
    cdr3_mode: str,
    min_v_identity: float,
    min_j_identity: float,
    require_productive: bool,
    adapter_5: str,
    adapter_3: str,
    logger: Logger,
) -> None:
    group_calls_tsv = os.path.join(out_dir, "umi_group_calls.tsv")
    clonotype_tsv = os.path.join(out_dir, "clonotype.summary.tsv")
    qc_txt = os.path.join(out_dir, "qc.summary.txt")

    group_calls_sorted = sorted(group_calls, key=lambda x: (-x.called, -x.read_count, x.group_id))
    header = ["group_id", "read_count", "productive_reads", "called", "V_gene"]
    if include_d:
        header += ["D_gene"]
    header += ["J_gene", "CDR3_tag", "support_reads", "support_frac", "note"]

    rows_out: List[List[str]] = []
    for c in group_calls_sorted:
        row = [
            c.group_id, str(c.read_count), str(c.productive_reads),
            "1" if c.called else "0",
            c.v_gene,
        ]
        if include_d:
            row.append(c.d_gene)
        row += [c.j_gene, c.cdr3_tag, str(c.support), f"{c.support_frac:.6f}", c.note]
        rows_out.append(row)

    write_tsv(group_calls_tsv, header, rows_out)
    logger.info(f"Wrote per-UMI calls: {group_calls_tsv}")

    total_umi_called = sum(umi for umi, _reads, _k in items) or 1
    header2 = ["rank", "UMI_count", "UMI_fraction", "read_count_sum", "V_gene"]
    if include_d:
        header2 += ["D_gene"]
    header2 += ["J_gene", "CDR3_tag"]

    rows2 = []
    for i, (umi, reads, key) in enumerate(items[:args.top_k], start=1):
        if include_d:
            v, d, j, cdr3 = key
            rows2.append([str(i), str(umi), f"{umi/total_umi_called:.6f}", str(reads), v, d, j, cdr3])
        else:
            v, j, cdr3 = key
            rows2.append([str(i), str(umi), f"{umi/total_umi_called:.6f}", str(reads), v, j, cdr3])

    write_tsv(clonotype_tsv, header2, rows2)
    logger.info(f"Wrote clonotype summary: {clonotype_tsv}")

    called_groups = sum(1 for c in group_calls if c.called)
    with open(qc_txt, "w") as f:
        f.write("mode\tperumi\n")
        f.write(f"total_groups\t{total_groups}\n")
        f.write(f"called_groups\t{called_groups}\n")
        f.write(f"min_reads_per_umi\t{args.min_reads_per_umi}\n")
        f.write(f"min_support\t{args.min_support}\n")
        f.write(f"require_productive\t{int(require_productive)}\n")
        f.write(f"total_reads_merged_raw\t{total_reads_written}\n")
        f.write(f"cutadapt_trimmed_reads\t{trimmed_n}\n")
        f.write(f"cutadapt_untrimmed_reads\t{untrimmed_n}\n")
        f.write(f"cutadapt_error_rate\t{error_rate}\n")
        f.write(f"cutadapt_revcomp\t{int(enable_revcomp)}\n")
        f.write(f"min_v_identity\t{min_v_identity}\n")
        f.write(f"min_j_identity\t{min_j_identity}\n")
        f.write(f"include_d\t{int(include_d)}\n")
        f.write(f"use_allele\t{int(use_allele)}\n")
        f.write(f"cdr3_mode\t{cdr3_mode}\n")
        f.write(f"adapter_5\t{adapter_5}\n")
        f.write(f"adapter_3\t{adapter_3}\n")

    logger.info(f"QC summary written: {qc_txt}")


def select_consensus_items(
    items: List[Tuple[int, int, tuple]],
    top_k: Optional[int],
    top_p: Optional[float],
) -> Tuple[
    List[Tuple[int, int, tuple]],
    Optional[int],
    Optional[float],
    List[List[str]],
]:
    """Select clonotypes for consensus using rank and/or UMI percentage.

    If neither selector is supplied, the combined defaults are top 3 and
    UMI_fraction > 10%. If only one selector is supplied, only that criterion
    is applied. When both are supplied, a clonotype must pass both.
    """
    if top_k is None and top_p is None:
        effective_top_k: Optional[int] = 3
        effective_top_p: Optional[float] = 10.0
    else:
        effective_top_k = int(top_k) if top_k is not None else None
        effective_top_p = float(top_p) if top_p is not None else None

    if effective_top_k is not None and effective_top_k < 1:
        raise ValueError("consensus_top_k must be >= 1")
    if effective_top_p is not None and not 0 <= effective_top_p <= 100:
        raise ValueError("consensus_top_p must be between 0 and 100")

    total_umi = sum(umi for umi, _reads, _key in items)
    selected: List[Tuple[int, int, tuple]] = []
    selection_rows: List[List[str]] = []

    for rank, item in enumerate(items, start=1):
        umi_count, read_count, key = item
        fraction = (umi_count / total_umi) if total_umi else 0.0
        percent = fraction * 100.0
        rank_pass = effective_top_k is None or rank <= effective_top_k
        percent_pass = effective_top_p is None or percent > effective_top_p
        is_selected = rank_pass and percent_pass
        if is_selected:
            selected.append(item)
        selection_rows.append([
            str(rank),
            str(umi_count),
            f"{fraction:.6f}",
            str(read_count),
            "|".join(str(x) for x in key),
            "1" if rank_pass else "0",
            "1" if percent_pass else "0",
            "1" if is_selected else "0",
        ])

    return selected, effective_top_k, effective_top_p, selection_rows


def _run_consensus_pipeline(
    args,
    igblastn_exe: str,
    igblast_tsv: str,
    query_raw_fq: str,
    query_trim_fq: str,
    out_dir: str,
    items: List[Tuple[int, int, tuple]],
    group_to_key: Dict[str, tuple],
    include_d: bool,
    use_allele: bool,
    cdr3_mode: str,
    min_v_identity: float,
    min_j_identity: float,
    require_productive: bool,
    logger: Logger,
) -> None:
    if not bool(getattr(args, "do_consensus", False)):
        logger.info("----- vault vdj finished (no consensus) -----")
        return

    logger.info("Consensus mode enabled: building selected clonotype consensus sequences ...")

    selected_items, effective_top_k, effective_top_p, selection_rows = select_consensus_items(
        items=items,
        top_k=getattr(args, "consensus_top_k", None),
        top_p=getattr(args, "consensus_top_p", None),
    )
    base_cons_dir = os.path.join(out_dir, "consensus")
    os.makedirs(base_cons_dir, exist_ok=True)
    write_tsv(
        os.path.join(base_cons_dir, "consensus.selection.tsv"),
        [
            "rank",
            "UMI_count",
            "UMI_fraction",
            "read_count_sum",
            "clonotype_key",
            "rank_pass",
            "percent_pass",
            "selected",
        ],
        selection_rows,
    )

    criteria = []
    if effective_top_k is not None:
        criteria.append(f"rank <= {effective_top_k}")
    if effective_top_p is not None:
        criteria.append(f"UMI_fraction > {effective_top_p:g}%")
    logger.info(
        f"Consensus selection criteria: {' AND '.join(criteria)}; "
        f"selected={len(selected_items)}/{len(items)}"
    )

    if not selected_items:
        logger.warn(
            "No clonotype satisfies the consensus selection criteria; "
            "skip consensus generation."
        )
        write_tsv(
            os.path.join(base_cons_dir, "consensus.summary.tsv"),
            [
                "clonotype_id",
                "policy",
                "strict_amplicon_reads",
                "strict_vdj_reads",
                "amplicon_reads",
                "vdj_reads",
                "amplicon_consensus_fasta",
                "vdj_consensus_fasta",
                "vdj_from_amplicon_consensus_fasta",
            ],
            [],
        )
        with open(
            os.path.join(base_cons_dir, "consensus.summary.T.tsv"),
            "w",
        ) as f:
            f.write("field\n")
        logger.info("----- vault vdj finished (no clonotype selected for consensus) -----")
        return

    minimap2_exe = resolve_exe("minimap2")
    racon_exe = resolve_exe("racon")
    if not minimap2_exe or not racon_exe:
        raise RuntimeError("Consensus requires minimap2 and racon in PATH. Install them or disable --do_consensus.")

    medaka_enabled = bool(getattr(args, "enable_medaka", False))
    medaka_exe: Optional[str] = None
    if medaka_enabled:
        user_medaka = getattr(args, "medaka_consensus", None)
        medaka_exe = resolve_exe(user_medaka) if user_medaka else resolve_exe("medaka_consensus")
        if not medaka_exe:
            raise RuntimeError("enable_medaka is set but medaka_consensus not found.")
        logger.info(f"Using medaka_consensus: {medaka_exe}")

    consensus_threads = int(getattr(args, "consensus_threads", 8))
    racon_rounds = int(getattr(args, "racon_rounds", 2))
    medaka_model = str(getattr(args, "medaka_model", "auto")).strip()
    consensus_fallback_min_reads = int(getattr(args, "consensus_fallback_min_reads", 10))

    top_keys = [key for _umi, _reads, key in selected_items]

    strict_amp = Counter()
    strict_vdj = Counter()
    if require_productive:
        read_to_key_strict, read_to_coords_strict, read_to_rc_strict = build_read_hit_maps(
            airr_tsv=igblast_tsv,
            require_productive=True,
            include_d=include_d,
            use_allele=use_allele,
            cdr3_mode=cdr3_mode,
            min_v_identity=min_v_identity,
            min_j_identity=min_j_identity,
        )
        read_to_key_relaxed, read_to_coords_relaxed, read_to_rc_relaxed = build_read_hit_maps(
            airr_tsv=igblast_tsv,
            require_productive=False,
            include_d=include_d,
            use_allele=use_allele,
            cdr3_mode=cdr3_mode,
            min_v_identity=min_v_identity,
            min_j_identity=min_j_identity,
        )
        logger.info(
            f"Consensus policy base=STRICT_PRODUCTIVE (require_productive=1). "
            f"Fallback threshold={consensus_fallback_min_reads} (applies to BOTH amplicon+VDJ)."
        )
    else:
        read_to_key_strict, read_to_coords_strict, read_to_rc_strict = {}, {}, {}
        read_to_key_relaxed, read_to_coords_relaxed, read_to_rc_relaxed = build_read_hit_maps(
            airr_tsv=igblast_tsv,
            require_productive=False,
            include_d=include_d,
            use_allele=use_allele,
            cdr3_mode=cdr3_mode,
            min_v_identity=min_v_identity,
            min_j_identity=min_j_identity,
        )
        logger.info("Consensus policy=NONPRODUCTIVE_ALLOWED (require_productive=0). No fallback needed.")

    key_to_outdir: Dict[tuple, str] = {}
    key_to_amp_fq: Dict[tuple, str] = {}
    key_to_vdj_fq: Dict[tuple, str] = {}

    for k in top_keys:
        if include_d:
            v, d, j, cdr3 = k
            cid = f"{v}__{d}__{j}__{cdr3}"
        else:
            v, j, cdr3 = k
            cid = f"{v}__{j}__{cdr3}"
        cid = re.sub(r"[^A-Za-z0-9_.-]+", "_", cid)[:180]
        outd = os.path.join(base_cons_dir, cid)
        os.makedirs(outd, exist_ok=True)
        key_to_outdir[k] = outd
        key_to_amp_fq[k] = os.path.join(outd, "reads.amplicon.fastq.gz")
        key_to_vdj_fq[k] = os.path.join(outd, "reads.vdj.fastq.gz")

    if require_productive:
        top_set = set(top_keys)
        for rid, seq, qual in iter_fastq_records(query_trim_fq):
            group_id = rid.split("|", 1)[0]
            k = group_to_key.get(group_id)
            if k not in top_set:
                continue
            if read_to_key_strict.get(rid) == k:
                strict_amp[k] += 1
                coords = read_to_coords_strict.get(rid)
                if coords is not None:
                    if read_to_rc_strict.get(rid, False):
                        seq, qual = revcomp_seq_qual(seq, qual)
                    if crop_vdj(seq, qual, coords[0], coords[1]) is not None:
                        strict_vdj[k] += 1

    key_policy: Dict[tuple, str] = {}
    if require_productive:
        for k in top_keys:
            sa = strict_amp.get(k, 0)
            sv = strict_vdj.get(k, 0)
            if sa < consensus_fallback_min_reads or sv < consensus_fallback_min_reads:
                key_policy[k] = "FALLBACK_NONPRODUCTIVE_ALLOWED"
                logger.warn(
                    f"Consensus fallback for {os.path.basename(key_to_outdir[k])}: "
                    f"strict_amplicon_reads={sa}, strict_vdj_reads={sv} < {consensus_fallback_min_reads}"
                )
            else:
                key_policy[k] = "STRICT_PRODUCTIVE"
                logger.info(
                    f"Consensus strict for {os.path.basename(key_to_outdir[k])}: "
                    f"strict_amplicon_reads={sa}, strict_vdj_reads={sv}"
                )
    else:
        key_policy = {k: "NONPRODUCTIVE_ALLOWED" for k in top_keys}

    amp_handles: Dict[tuple, gzip.GzipFile] = {k: gzip.open(key_to_amp_fq[k], "wt") for k in top_keys}
    vdj_handles: Dict[tuple, gzip.GzipFile] = {k: gzip.open(key_to_vdj_fq[k], "wt") for k in top_keys}

    written_amp = Counter()
    written_vdj = Counter()
    top_set = set(top_keys)

    for rid, seq, qual in iter_fastq_records(query_raw_fq):
        group_id = rid.split("|", 1)[0]
        k = group_to_key.get(group_id)
        if k not in top_set:
            continue

        policy = key_policy.get(k, "NONPRODUCTIVE_ALLOWED")
        if require_productive and policy == "STRICT_PRODUCTIVE":
            if read_to_key_strict.get(rid) != k:
                continue
            rcflag = read_to_rc_strict.get(rid, False)
        else:
            if read_to_key_relaxed.get(rid) != k:
                continue
            rcflag = read_to_rc_relaxed.get(rid, False)

        if rcflag:
            seq, qual = revcomp_seq_qual(seq, qual)

        amp_handles[k].write(f"@{rid}\n{seq}\n+\n{qual}\n")
        written_amp[k] += 1

    for rid, seq, qual in iter_fastq_records(query_trim_fq):
        group_id = rid.split("|", 1)[0]
        k = group_to_key.get(group_id)
        if k not in top_set:
            continue

        policy = key_policy.get(k, "NONPRODUCTIVE_ALLOWED")
        if require_productive and policy == "STRICT_PRODUCTIVE":
            if read_to_key_strict.get(rid) != k:
                continue
            coords = read_to_coords_strict.get(rid)
            rcflag = read_to_rc_strict.get(rid, False)
        else:
            if read_to_key_relaxed.get(rid) != k:
                continue
            coords = read_to_coords_relaxed.get(rid)
            rcflag = read_to_rc_relaxed.get(rid, False)

        if coords is None:
            continue
        if rcflag:
            seq, qual = revcomp_seq_qual(seq, qual)

        cropped = crop_vdj(seq, qual, coords[0], coords[1])
        if cropped is None:
            continue
        s2, q2 = cropped
        vdj_handles[k].write(f"@{rid}\n{s2}\n+\n{q2}\n")
        written_vdj[k] += 1

    for h in amp_handles.values():
        h.close()
    for h in vdj_handles.values():
        h.close()

    logger.info(
        f"Extracted reads into {len(top_keys)} selected clonotype bins under: "
        f"{base_cons_dir}"
    )
    for k in top_keys:
        logger.info(
            f"  {os.path.basename(key_to_outdir[k])}\tpolicy={key_policy[k]}"
            f"\tamplicon_reads_written(raw)={written_amp.get(k, 0)}"
            f"\tvdj_reads_written(trimmed+cropped)={written_vdj.get(k, 0)}"
        )

    summary_rows = []
    for k in top_keys:
        outd = key_to_outdir[k]
        amp_fq = key_to_amp_fq[k]
        vdj_fq = key_to_vdj_fq[k]
        cid = os.path.basename(outd)

        amp_n = written_amp.get(k, 0)
        vdj_n = written_vdj.get(k, 0)

        amp_cons = ""
        vdj_cons = ""
        vdj_from_amp_cons = ""

        logger.info(f"Consensus for clonotype: {cid} (policy={key_policy[k]})")

        if amp_n >= 1:
            amp_cons = _run_consensus_pair(
                reads_fq=amp_fq,
                outdir=outd,
                label="amplicon",
                minimap2_exe=minimap2_exe,
                racon_exe=racon_exe,
                threads=consensus_threads,
                racon_rounds=racon_rounds,
                logger=logger,
                medaka_enabled=medaka_enabled,
                medaka_exe=medaka_exe,
                medaka_model=medaka_model,
            )
        else:
            logger.warn("  amplicon reads=0 => skip amplicon consensus")

        if vdj_n >= 1:
            vdj_cons = _run_consensus_pair(
                reads_fq=vdj_fq,
                outdir=outd,
                label="vdj",
                minimap2_exe=minimap2_exe,
                racon_exe=racon_exe,
                threads=consensus_threads,
                racon_rounds=racon_rounds,
                logger=logger,
                medaka_enabled=medaka_enabled,
                medaka_exe=medaka_exe,
                medaka_model=medaka_model,
            )
        else:
            logger.warn("  vdj reads=0 => skip vdj consensus")

        expected_v = k[0]
        expected_d = k[1] if include_d else ""
        expected_j = k[2] if include_d else k[1]
        expected_cdr3 = k[3] if include_d else k[2]

        amp_row = None
        vdj_row = None
        if amp_cons and os.path.exists(amp_cons):
            amp_row, _ = normalize_consensus_orientation_by_igblast(
                cons_fa=amp_cons,
                label="amplicon",
                outdir=outd,
                igblastn_exe=igblastn_exe,
                germ_v=args.germV,
                germ_d=args.germD,
                germ_j=args.germJ,
                auxiliary_data=args.auxiliary_data,
                organism=args.organism,
                logger=logger,
            )
        if vdj_cons and os.path.exists(vdj_cons):
            vdj_row, _ = normalize_consensus_orientation_by_igblast(
                cons_fa=vdj_cons,
                label="vdj",
                outdir=outd,
                igblastn_exe=igblastn_exe,
                germ_v=args.germV,
                germ_d=args.germD,
                germ_j=args.germJ,
                auxiliary_data=args.auxiliary_data,
                organism=args.organism,
                logger=logger,
            )

        # NEW: vdj_from_amplicon.consensus.fa
        if amp_cons and amp_row is not None and os.path.exists(amp_cons):
            out_vdj_from_amp = os.path.join(outd, "vdj_from_amplicon.consensus.fa")
            got = extract_vdj_from_amplicon_consensus(
                amplicon_cons_fa=amp_cons,
                amp_row=amp_row,
                out_fa=out_vdj_from_amp,
                logger=logger,
            )
            if got:
                vdj_from_amp_cons = got

        if amp_cons and vdj_cons and os.path.exists(amp_cons) and os.path.exists(vdj_cons):
            _hA, amp_seq = read_fasta_one(amp_cons)
            _hV, vdj_seq = read_fasta_one(vdj_cons)

            vdj_in_amp = (vdj_seq in amp_seq)
            vdjRC_in_amp = (revcomp(vdj_seq) in amp_seq)

            if (not vdj_in_amp) and vdjRC_in_amp:
                revcomp_fasta_inplace(vdj_cons, tag="revcomp_to_match_amplicon")
                _hV, vdj_seq = read_fasta_one(vdj_cons)
                vdj_in_amp = (vdj_seq in amp_seq)
                vdjRC_in_amp = (revcomp(vdj_seq) in amp_seq)
                logger.warn(f"[QC] {cid} vdj was reverse of amplicon by substring => flipped vdj.consensus.fa")

            if not vdj_in_amp:
                strand = minimap2_best_strand(target_fa=amp_cons, query_fa=vdj_cons, minimap2_exe=minimap2_exe)
                if strand == "-":
                    revcomp_fasta_inplace(vdj_cons, tag="revcomp_to_match_amplicon_minimap2")
                    logger.warn(f"[QC] {cid} minimap2 says vdj vs amplicon is '-' => flipped vdj.consensus.fa")
                    _hV, vdj_seq = read_fasta_one(vdj_cons)
                    vdj_in_amp = (vdj_seq in amp_seq)
                    vdjRC_in_amp = (revcomp(vdj_seq) in amp_seq)

            logger.info(
                f"[QC] {cid}\tcontainment"
                f"\tlen_amp={len(amp_seq)}\tlen_vdj={len(vdj_seq)}"
                f"\tvdj_in_amp={int(vdj_in_amp)}\tvdjRC_in_amp={int(vdjRC_in_amp)}"
            )

        def _format_match(row: Optional[Dict[str, str]], label: str) -> str:
            if not row:
                return f"[QC] {cid}\tigblast_{label}\tNO_CALL"
            v_call = collapse_gene((row.get("v_call") or "").strip(), use_allele=use_allele)
            j_call = collapse_gene((row.get("j_call") or "").strip(), use_allele=use_allele)
            d_call = collapse_gene((row.get("d_call") or "").strip(), use_allele=use_allele) if include_d else ""
            ctag = _cdr3_tag_from_row(row, cdr3_mode=cdr3_mode)
            rc = (row.get("rev_comp") or "").strip()

            ok = (v_call == expected_v) and (j_call == expected_j) and (ctag == expected_cdr3)
            if include_d:
                ok = ok and (d_call == expected_d)

            # return (
            #     f"[QC] {cid}\tigblast_{label}"
            #     f"\trev_comp={rc}"
            #     f"\tV={v_call}\tD={d_call}\tJ={j_call}\tCDR3={ctag}"
            #     f"\tmatch_expected={int(ok)}"
            # )
            fields = [
                f"[QC] {cid}",
                f"igblast_{label}",
                f"rev_comp={rc}",
                f"V={v_call}",
            ]
            if include_d:
                fields.append(f"D={d_call}")
            fields += [
                f"J={j_call}",
                f"CDR3={ctag}",
                f"match_expected={int(ok)}",
            ]
            return "\t".join(fields)

        if amp_cons:
            logger.info(_format_match(amp_row, "amplicon"))
        if vdj_cons:
            logger.info(_format_match(vdj_row, "vdj"))

        # optional: igblast check on vdj_from_amplicon (helps sanity check / “健全性检查”)
        if vdj_from_amp_cons and os.path.exists(vdj_from_amp_cons):
            row_from_amp = None
            try:
                row_from_amp = igblast_first_row(
                    igblastn_exe=igblastn_exe,
                    query_fa=vdj_from_amp_cons,
                    out_tsv=os.path.join(outd, "qc.igblast.vdj_from_amplicon.airr.tsv"),
                    germ_v=args.germV,
                    germ_d=args.germD,
                    germ_j=args.germJ,
                    auxiliary_data=args.auxiliary_data,
                    organism=args.organism,
                    threads=1,
                )
            except Exception as e:
                logger.warn(f"[QC] igblast failed on vdj_from_amplicon consensus: {e}")
            logger.info(_format_match(row_from_amp, "vdj_from_amplicon"))

        summary_rows.append([
            cid,
            key_policy[k],
            str(strict_amp.get(k, 0) if require_productive else 0),
            str(strict_vdj.get(k, 0) if require_productive else 0),
            str(amp_n),
            str(vdj_n),
            amp_cons or "",
            vdj_cons or "",
            vdj_from_amp_cons or "",
        ])

    # write_tsv(
    #     os.path.join(base_cons_dir, "consensus.summary.tsv"),
    #     [
    #         "clonotype_id",
    #         "policy",
    #         "strict_amplicon_reads",
    #         "strict_vdj_reads",
    #         "amplicon_reads",
    #         "vdj_reads",
    #         "amplicon_consensus_fasta",
    #         "vdj_consensus_fasta",
    #         "vdj_from_amplicon_consensus_fasta",
    #     ],
    #     summary_rows,
    # )

    def write_tsv_transposed(path: str, header: List[str], rows: List[List[str]]) -> None:
        """
        Transpose a TSV table (rows -> columns) for easier shell viewing.
        First column becomes 'field', columns are clonotype_id(s).
        """
        os.makedirs(os.path.dirname(path), exist_ok=True)
        if not rows:
            with open(path, "w") as f:
                f.write("field\n")
            return

        # clonotype_id column is the first column in rows
        clonotype_ids = [r[0] if len(r) > 0 else "" for r in rows]

        # ensure unique column names (in case of duplicates)
        seen: Dict[str, int] = {}
        cols: List[str] = []
        for cid in clonotype_ids:
            seen[cid] = seen.get(cid, 0) + 1
            cols.append(f"{cid}#{seen[cid]}" if seen[cid] > 1 else cid)

        with open(path, "w") as f:
            f.write("field\t" + "\t".join(cols) + "\n")

            # each original header field becomes one row
            for j, field in enumerate(header[1:], start=1):
                vals = []
                for r in rows:
                    vals.append(r[j] if j < len(r) else "")
                f.write(field + "\t" + "\t".join(vals) + "\n")

    summary_header = [
        "clonotype_id",
        "policy",
        "strict_amplicon_reads",
        "strict_vdj_reads",
        "amplicon_reads",
        "vdj_reads",
        "amplicon_consensus_fasta",
        "vdj_consensus_fasta",
        "vdj_from_amplicon_consensus_fasta",
    ]
    summary_path = os.path.join(base_cons_dir, "consensus.summary.tsv")
    write_tsv(summary_path, summary_header, summary_rows)

    # NEW: transpose version for shell
    write_tsv_transposed(os.path.join(base_cons_dir, "consensus.summary.T.tsv"), summary_header, summary_rows)

    logger.info("----- vault vdj finished (with consensus) -----")


# ---------------------------
# Main entry (short)
# ---------------------------

def vdj_main(args) -> None:
    out_dir, logger = _setup_outdir_and_logger(args)
    logger.info("----- vault vdj started -----")

    _umi_adapter, adapter_5, adapter_3 = _resolve_adapter_params(args, logger)
    igblastn_exe, cutadapt_cmd = _resolve_tools_and_validate(args, logger)
    mode = _determine_mode(args, logger)

    query_raw_fq, group_read_counts_raw, total_reads_written = _build_query_raw_fastq(args, mode, out_dir, logger)

    (
        query_trim_fq,
        _query_untrim_fq,
        query_fa,
        group_trimmed_counts,
        trimmed_n,
        untrimmed_n,
        enable_revcomp,
        error_rate,
    ) = _run_cutadapt_and_make_query_fa(
        args=args,
        mode=mode,
        query_raw_fq=query_raw_fq,
        adapter_5=adapter_5,
        adapter_3=adapter_3,
        cutadapt_cmd=cutadapt_cmd,
        out_dir=out_dir,
        logger=logger,
    )

    igblast_tsv = _run_igblast_main(args, igblastn_exe, query_fa, out_dir, logger)

    include_d = bool(getattr(args, "include_d", False))
    use_allele = bool(getattr(args, "use_allele", False))
    cdr3_mode = str(getattr(args, "cdr3_mode", "aa"))
    min_v_identity = float(getattr(args, "min_v_identity", 0.0))
    min_j_identity = float(getattr(args, "min_j_identity", 0.0))
    require_productive = bool(getattr(args, "require_productive", False))

    if mode == "raw":
        _raw_mode_report(
            args=args,
            igblast_tsv=igblast_tsv,
            out_dir=out_dir,
            logger=logger,
            include_d=include_d,
            use_allele=use_allele,
            cdr3_mode=cdr3_mode,
            min_v_identity=min_v_identity,
            min_j_identity=min_j_identity,
            require_productive=require_productive,
            trimmed_n=trimmed_n,
            untrimmed_n=untrimmed_n,
            error_rate=error_rate,
            enable_revcomp=enable_revcomp,
        )
        logger.info("----- vault vdj finished -----")
        return

    if group_trimmed_counts is None:
        group_trimmed_counts = Counter()

    group_calls, items, group_to_key, total_groups = _perumi_calls_and_summary(
        args=args,
        igblast_tsv=igblast_tsv,
        group_read_counts_raw=group_read_counts_raw,
        group_trimmed_counts=group_trimmed_counts,
        include_d=include_d,
        use_allele=use_allele,
        cdr3_mode=cdr3_mode,
        min_v_identity=min_v_identity,
        min_j_identity=min_j_identity,
        require_productive=require_productive,
        logger=logger,
    )

    _write_perumi_outputs(
        args=args,
        out_dir=out_dir,
        group_calls=group_calls,
        items=items,
        total_groups=total_groups,
        total_reads_written=total_reads_written,
        trimmed_n=trimmed_n,
        untrimmed_n=untrimmed_n,
        error_rate=error_rate,
        enable_revcomp=enable_revcomp,
        include_d=include_d,
        use_allele=use_allele,
        cdr3_mode=cdr3_mode,
        min_v_identity=min_v_identity,
        min_j_identity=min_j_identity,
        require_productive=require_productive,
        adapter_5=adapter_5,
        adapter_3=adapter_3,
        logger=logger,
    )

    _run_consensus_pipeline(
        args=args,
        igblastn_exe=igblastn_exe,
        igblast_tsv=igblast_tsv,
        query_raw_fq=query_raw_fq,
        query_trim_fq=query_trim_fq,
        out_dir=out_dir,
        items=items,
        group_to_key=group_to_key,
        include_d=include_d,
        use_allele=use_allele,
        cdr3_mode=cdr3_mode,
        min_v_identity=min_v_identity,
        min_j_identity=min_j_identity,
        require_productive=require_productive,
        logger=logger,
    )


__all__ = ["vdj_main"]
