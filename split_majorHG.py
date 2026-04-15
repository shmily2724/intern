#!/usr/bin/env python3
import argparse
import re
from pathlib import Path
from collections import defaultdict, Counter


def read_fasta(path):
    """
    Read FASTA safely.
    - keeps all sequence characters exactly as they are (? - N etc.)
    - only strips line breaks / leading-trailing whitespace
    - validates basic FASTA format
    """
    records = []
    header = None
    seq_lines = []
    seen_headers = set()

    with open(path, "r", encoding="utf-8") as fh:
        for lineno, raw in enumerate(fh, start=1):
            line = raw.rstrip("\n\r")

            # skip completely empty lines
            if not line.strip():
                continue

            if line.startswith(">"):
                if header is not None:
                    if not seq_lines:
                        raise ValueError(f"{path}: header '{header}' has no sequence.")
                    records.append((header, "".join(seq_lines)))

                header = line[1:].strip()
                if not header:
                    raise ValueError(f"{path}: empty FASTA header at line {lineno}.")
                if header in seen_headers:
                    raise ValueError(f"{path}: duplicated FASTA header '{header}'.")
                seen_headers.add(header)
                seq_lines = []
            else:
                if header is None:
                    raise ValueError(f"{path}: sequence found before first header at line {lineno}.")
                seq = line.strip()
                if any(ch.isspace() for ch in seq):
                    raise ValueError(f"{path}: whitespace inside sequence line at line {lineno}.")
                seq_lines.append(seq)

    if header is not None:
        if not seq_lines:
            raise ValueError(f"{path}: header '{header}' has no sequence.")
        records.append((header, "".join(seq_lines)))

    if not records:
        raise ValueError(f"{path}: no FASTA records found.")

    return records


def write_fasta(records, out_path, width=80):
    """
    Write FASTA with standard wrapping.
    Sequence characters are preserved exactly; only line wrapping is applied.
    """
    with open(out_path, "w", encoding="utf-8") as out:
        for header, seq in records:
            out.write(f">{header}\n")
            for i in range(0, len(seq), width):
                out.write(seq[i:i + width] + "\n")


def extract_sample_id_from_original(header):
    """
    Example:
    MH448947.1 Homo sapiens isolate CoLao566 haplogroup M10 mitochondrion, complete genome
    -> CoLao566
    """
    m = re.search(r"\bisolate\s+([A-Za-z][A-Za-z0-9_-]*)\b", header, flags=re.IGNORECASE)
    return m.group(1) if m else None


def extract_haplogroup_from_original(header):
    """
    Example:
    ... haplogroup M10 mitochondrion ...
    -> M10
    """
    m = re.search(r"\bhaplogroup\s+([A-Za-z][A-Za-z0-9+._-]*)\b", header, flags=re.IGNORECASE)
    return m.group(1) if m else None


def major_haplogroup(full_hg):
    """
    M10 -> M
    F1a1a -> F
    A+152+16362 -> A
    """
    if not full_hg:
        return None
    m = re.match(r"([A-Za-z])", full_hg)
    return m.group(1).upper() if m else None


def split_fasta_by_major_hg(original_fasta, reorder_fasta, outdir, line_width=80):
    original_fasta = Path(original_fasta)
    reorder_fasta = Path(reorder_fasta)
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    original_records = read_fasta(original_fasta)
    reorder_records = read_fasta(reorder_fasta)

    # build sample_id -> full haplogroup map from original FASTA
    sample_to_hg = {}
    missing_hg_in_original = []

    for header, _seq in original_records:
        sample_id = extract_sample_id_from_original(header)
        if sample_id is None:
            raise ValueError(f"Cannot parse sample ID from original header:\n{header}")

        haplogroup = extract_haplogroup_from_original(header)

        if sample_id in sample_to_hg:
            raise ValueError(f"Duplicated sample ID in original FASTA: {sample_id}")

        sample_to_hg[sample_id] = haplogroup

        if haplogroup is None:
            missing_hg_in_original.append(sample_id)

    # validate reorder headers are unique
    reorder_ids = [header.strip() for header, _seq in reorder_records]
    dup_reorder = [x for x, c in Counter(reorder_ids).items() if c > 1]
    if dup_reorder:
        raise ValueError(
            "Duplicated headers in reorder FASTA:\n" + "\n".join(sorted(dup_reorder))
        )

    # validate alignment lengths in reorder FASTA
    aligned_lengths = {len(seq) for _header, seq in reorder_records}
    if len(aligned_lengths) != 1:
        raise ValueError(
            "Sequences in reorder FASTA do not all have the same length:\n"
            + ", ".join(map(str, sorted(aligned_lengths)))
        )
    aligned_length = next(iter(aligned_lengths))

    grouped = defaultdict(list)
    unresolved = []

    for header, seq in reorder_records:
        short_header = header.strip()   # keep short header exactly from reorder
        full_hg = sample_to_hg.get(short_header)
        major = major_haplogroup(full_hg)

        if major is None:
            unresolved.append((short_header, seq))
        else:
            grouped[major].append((short_header, seq))

    # write one FASTA per major haplogroup
    for major, records in sorted(grouped.items()):
        out_path = outdir / f"major_{major}.fasta"
        write_fasta(records, out_path, width=line_width)

    # always create unresolved.fasta
    unresolved_path = outdir / "unresolved.fasta"
    if unresolved:
        write_fasta(unresolved, unresolved_path, width=line_width)
    else:
        unresolved_path.write_text("", encoding="utf-8")

    # summary to screen
    print("Done.")
    print(f"Original FASTA : {original_fasta}")
    print(f"Reorder FASTA  : {reorder_fasta}")
    print(f"Output dir     : {outdir}")
    print(f"Reorder records: {len(reorder_records)}")
    print(f"Aligned length : {aligned_length}")
    print()

    for major in sorted(grouped):
        print(f"major_{major}.fasta\t{len(grouped[major])}")

    print(f"unresolved.fasta\t{len(unresolved)}")

    if unresolved:
        print("\nUnresolved sample IDs:")
        for header, _seq in unresolved:
            print(f"  - {header}")

    if missing_hg_in_original:
        print("\nSamples present in original FASTA but lacking haplogroup in original header:")
        for sample_id in sorted(missing_hg_in_original):
            print(f"  - {sample_id}")


def main():
    parser = argparse.ArgumentParser(
        description="Split MAFFT-aligned reorder FASTA into one FASTA per major haplogroup."
    )
    parser.add_argument(
        "--original",
        required=True,
        help="Original FASTA with full header containing sample ID and haplogroup."
    )
    parser.add_argument(
        "--reorder",
        required=True,
        help="Reordered/aligned FASTA with short headers (e.g. >Kinh01)."
    )
    parser.add_argument(
        "--outdir",
        required=True,
        help="Output directory for major_*.fasta files."
    )
    parser.add_argument(
        "--width",
        type=int,
        default=80,
        help="Line width for FASTA output (default: 80)."
    )

    args = parser.parse_args()

    if args.width <= 0:
        raise ValueError("--width must be > 0")

    split_fasta_by_major_hg(
        original_fasta=args.original,
        reorder_fasta=args.reorder,
        outdir=args.outdir,
        line_width=args.width
    )


if __name__ == "__main__":
    main()
