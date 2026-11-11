#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
from pathlib import Path
import csv
import math
import subprocess
import tempfile

parser = argparse.ArgumentParser(description="Prepare Duroc positions with MAF >= 0.05.")
parser.add_argument("--du-frequency", type=Path, required=True)
parser.add_argument("--outdir", type=Path, required=True)
args = parser.parse_args()
INPUT = args.du_frequency
OUT_DIR = args.outdir
OUTPUT = OUT_DIR / "DU_MAF005_only_pos.txt"

def main():
    if not INPUT.is_file():
        raise FileNotFoundError(INPUT)

    OUT_DIR.mkdir(parents=True, exist_ok=True)

    scanned = 0
    qualified = 0
    invalid = 0
    repeated_headers = 0

    # Sort on disk to deduplicate positions without loading the full list into memory.
    with tempfile.TemporaryDirectory(
        prefix="DU_MAF005_", dir=OUT_DIR
    ) as temporary:
        temporary = Path(temporary)
        raw = temporary / "qualified_positions.txt"
        sorted_file = temporary / "unique_positions.txt"

        with INPUT.open() as source, raw.open("w") as dest:
            reader = csv.DictReader(source, delimiter="\t")
            required = {"CHROM", "POS", "MAF"}
            if not required.issubset(reader.fieldnames or []):
                raise ValueError(
                    f"Required columns: {required}. Available columns: {reader.fieldnames}"
                )

            for row in reader:
                # Skip repeated header rows.
                if row.get("CHROM", "").strip() == "CHROM":
                    repeated_headers += 1
                    continue

                scanned += 1
                try:
                    maf = float(row["MAF"])
                    pos = int(row["POS"])
                    chrom = row["CHROM"].strip()
                    if chrom.lower().startswith("chr"):
                        chrom = chrom[3:]
                    if chrom.endswith(".0"):
                        chrom = chrom[:-2]
                    if (
                        not math.isfinite(maf)
                        or not 0 <= maf <= 0.5
                        or pos < 1
                        or not chrom
                    ):
                        raise ValueError
                except (ValueError, TypeError, KeyError, AttributeError):
                    invalid += 1
                    continue

                if maf >= 0.05:
                    dest.write(f"{chrom}:{pos}\n")
                    qualified += 1

                if scanned % 5_000_000 == 0:
                    print(
                        f"Scanned {scanned:,} rows. "
                        f"Allele records passing the MAF threshold: {qualified:,}",
                        flush=True
                    )

        print("Sorting and deduplicating positions...", flush=True)
        import os
        env = os.environ.copy()
        env["LC_ALL"] = "C"
        subprocess.run(
            [
                "sort", "-u", "-T", str(temporary),
                "-o", str(sorted_file), str(raw)
            ],
            check=True, env=env
        )

        with sorted_file.open() as handle:
            unique_count = sum(1 for _ in handle)

        # Replace the final output only after sorting succeeds.
        sorted_file.replace(OUTPUT)

    summary = (
        f"Input\t{INPUT}\n"
        "Rule\tAny allele record at a position has MAF >= 0.05\n"
        "Additional_NS_filter\tNone\n"
        "Additional_QC_flag_filter\tNone\n"
        f"Scanned_records\t{scanned}\n"
        f"Qualified_allele_records\t{qualified}\n"
        f"Unique_retained_positions\t{unique_count}\n"
        f"Invalid_records_skipped\t{invalid}\n"
        f"Repeated_headers_skipped\t{repeated_headers}\n"
        f"Output\t{OUTPUT}\n"
    )
    (OUT_DIR / "DU_MAF005_only_summary.tsv").write_text(
        summary, encoding="utf-8"
    )

    print("\nCompleted:", flush=True)
    print(summary, flush=True)

if __name__ == "__main__":
    main()
