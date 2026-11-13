#!/usr/bin/env python3
"""Join per-sample ASE summaries to the exonic variant list."""
import argparse
import gzip
from pathlib import Path

import pandas as pd


def read_variants(path):
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as handle:
        for line in handle:
            if line.startswith("#CHROM\t") or line.startswith("CHROM\t"):
                columns = line.lstrip("#").rstrip("\r\n").split("\t")
                return pd.read_csv(handle, sep="\t", names=columns, header=None,
                                   usecols=["CHROM", "POS", "REF", "ALT"], dtype=str)
    raise ValueError(f"Missing VCF column header: {path}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--variants", type=Path, required=True)
    parser.add_argument("--sample-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    variants = read_variants(args.variants)
    variants.insert(0, "SNP_id", variants.CHROM + ":" + variants.POS + "_" +
                    variants.REF + "_" + variants.ALT)
    files = sorted(args.sample_dir.glob("*_ASEReadCounter_binom_test.txt"))
    if not files:
        raise ValueError("No per-sample ASE result tables found.")
    ids = [path.name.split("_", 1)[0] for path in files]
    if len(ids) != len(set(ids)):
        raise ValueError("Duplicate sample IDs extracted from filenames.")
    samples = []
    for path, sample in zip(files, ids):
        data = pd.read_csv(path, sep="\t", dtype=str)
        keys = data.contig + ":" + data.position + "_" + data.refAllele + "_" + data.altAllele
        result = pd.DataFrame({
            f"{sample}_info": (data.refCount + "_" + data.altCount + "_" + data.FDR).to_numpy(),
            f"{sample}_FDR": data.FDR.to_numpy(),
        }, index=pd.Index(keys, name="SNP_id"))
        if not result.index.is_unique:
            raise ValueError(f"Duplicate variant keys in {path}")
        samples.append(result.reindex(variants.SNP_id).reset_index(drop=True))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    pd.concat([variants, *samples], axis=1).to_csv(
        args.output, sep="\t", index=False, na_rep="NA"
    )


if __name__ == "__main__":
    main()
