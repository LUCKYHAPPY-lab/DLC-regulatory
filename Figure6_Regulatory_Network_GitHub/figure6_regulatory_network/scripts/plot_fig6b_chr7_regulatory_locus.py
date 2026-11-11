#!/usr/bin/env python3
"""Plot the five-layer locus view centered on chr7:23,797,335 C>T.

Show regional GWAS, gene models, cis-eQTLs and LD, MPRA allelic-skew
significance, and five chromatin tracks around three highlighted sites.
"""

import argparse
import re
from pathlib import Path

import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import ConnectionPatch



parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--gwas", type=Path, required=True)
parser.add_argument("--eqtl", type=Path, required=True)
parser.add_argument("--ld", type=Path, required=True)
parser.add_argument("--gtf", type=Path, required=True)
parser.add_argument("--mpra", type=Path, required=True)
parser.add_argument("--atac", type=Path, required=True)
parser.add_argument("--ctcf", type=Path, required=True)
parser.add_argument("--h3k27ac", type=Path, required=True)
parser.add_argument("--h3k4me1", type=Path, required=True)
parser.add_argument("--h3k4me3", type=Path, required=True)
parser.add_argument("--outdir", type=Path, required=True)
args = parser.parse_args()
import pyBigWig

WORK_DIR = args.outdir
GWAS_FILE = args.gwas
EQTL_FILE = args.eqtl
LD_FILE_B = args.ld
GTF_FILE = args.gtf
MPRA_RESULT_FILE = args.mpra
TARGET_GENES = ["ENSSSCG00000001396", "ENSSSCG00000001418", "ENSSSCG00000001395"]
BW_FILES = {
    "ATAC": str(args.atac),
    "CTCF": str(args.ctcf),
    "H3K27ac": str(args.h3k27ac),
    "H3K4me1": str(args.h3k4me1),
    "H3K4me3": str(args.h3k4me3),
}
for path in (GWAS_FILE, EQTL_FILE, LD_FILE_B, GTF_FILE, MPRA_RESULT_FILE, *BW_FILES.values()):
    if not Path(path).is_file():
        raise FileNotFoundError(f"Input file not found: {path}")
WORK_DIR.mkdir(parents=True, exist_ok=True)

TRACK_COLORS = {
    "ATAC": "#1B9E77",
    "CTCF": "#6A3D9A",
    "H3K27ac": "#E64B35",
    "H3K4me1": "#4DBBD5",
    "H3K4me3": "#D73027",
}

TRAIT = "M_LMDEP"
CHROM = "7"
CHROM_BW = "chr7"
TARGET_GENE = "ENSSSCG00000001396"
EQTL_LEAD_POS = 23_599_787   # Reference position for the DLC LD input.
EMVAR_POS = 23_797_335
EMVAR_LABEL = "chr7:23,797,335 C>T"
REGION_CENTER_POS = EMVAR_POS
FLANK_BP = 500_000
REGION_START = max(1, REGION_CENTER_POS - FLANK_BP)
REGION_END = REGION_CENTER_POS + FLANK_BP
GWAS_THRESHOLD = 5e-8
GWAS_SUGGESTIVE = 1e-5       # Suggestive GWAS threshold.
EQTL_FDR_THRESHOLD = 0.05

mpl.rcParams.update({"font.family": "Arial", "font.size": 9,
                     "pdf.fonttype": 42, "ps.fonttype": 42,
                     "svg.fonttype": "none"})

GREY = "#BDBDBD"
LEAD_BLUE = "#2166AC"
ORANGE = "#F39C12"
RED = "#D73027"
GREEN = "#238B45"
SUGGESTIVE = "#F4A582"       # Color for 5e-8 <= P < 1e-5.
EQTL_ONLY_BLUE = "#7FB3D5"
ASE_GREEN = "#2CA25F"        # Color reserved for ASE support.


# Utilities
def norm_chr(x):
    x = re.sub(r"^chr", "", str(x).strip(), flags=re.I)
    return x[:-2] if x.endswith(".0") else x


def find_col(cols, names, required=True):
    lookup = {str(c).strip().lower(): c for c in cols}
    for name in names:
        if name.lower() in lookup:
            return lookup[name.lower()]
    if required:
        raise RuntimeError(f"Expected one of {names}. Available columns: {list(cols)}")
    return None


def parse_snp(x):
    m = re.match(r"^(?:chr)?([^:_]+)[:_]([0-9]+)", str(x).strip(), flags=re.I)
    if m:
        return pd.Series([norm_chr(m.group(1)), int(m.group(2))])
    return pd.Series([pd.NA, pd.NA])


def gtf_attr(text, key):
    m = re.search(rf'{key} "([^"]+)"', str(text))
    return m.group(1) if m else pd.NA


def as_bool(series):
    return series.astype(str).str.strip().str.lower().isin(["true", "t", "1", "yes"])


# Use five discrete LD color bins.
LD_BREAKS = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]
LD_COLORS = ["#2166AC", "#67A9CF", "#FEE08B", "#FDAE61", "#D73027"]


def ld_color(r2):
    if pd.isna(r2):
        return GREY
    if r2 < 0.2:
        return LD_COLORS[0]
    if r2 < 0.4:
        return LD_COLORS[1]
    if r2 < 0.6:
        return LD_COLORS[2]
    if r2 < 0.8:
        return LD_COLORS[3]
    return LD_COLORS[4]


def require_nonempty(df, label):
    if df.empty:
        raise RuntimeError(f"No {label} records remain. Check chromosome, coordinates and column names.")


def load_bw_raw(bw_file, chrom, start, end):
    """Read unsmoothed per-base bigWig signal."""
    bw = pyBigWig.open(bw_file)
    try:
        values = bw.values(chrom, start, end)
    except Exception as e:
        print(f"  Warning: could not read {chrom}:{start}-{end} from {bw_file}: {e}")
        values = [0.0] * (end - start)
    finally:
        bw.close()
    xs = np.arange(start, end)
    ys = np.array([0.0 if v is None else float(v) for v in values], dtype=float)
    return xs, ys


def load_bw_signal(bw_file, chrom, start, end, bin_size=200):
    """Read bigWig signal averaged across bins."""
    bw = pyBigWig.open(bw_file)
    try:
        values = bw.values(chrom, start, end)
    except Exception as e:
        print(f"  Warning: could not read {chrom}:{start}-{end} from {bw_file}: {e}")
        values = [0] * (end - start)
    finally:
        bw.close()

    bins = np.arange(start, end, bin_size)
    score_vec = np.zeros(len(bins))
    for i, b in enumerate(bins):
        s = b - start
        e = min(s + bin_size, len(values))
        if s < len(values):
            chunk = values[s:e]
            chunk = [v for v in chunk if v is not None]
            if chunk:
                score_vec[i] = np.mean(chunk)
    return bins, score_vec


# 1. Load inputs
print("Loading inputs...")
print(f"  WORK_DIR = {WORK_DIR}")

# Load regional GWAS associations before filtering eQTLs.
print("Reading GWAS...")
gwas_parts = []
for chunk in pd.read_csv(GWAS_FILE, sep=r"\s+", compression="gzip", chunksize=1_000_000, low_memory=False):
    cc = find_col(chunk.columns, ["chromosome", "CHROM", "chr"])
    pc = find_col(chunk.columns, ["position", "POS", "bp"])
    pvc = find_col(chunk.columns, ["p_value", "pvalue", "P", "p.value"])
    rc = find_col(chunk.columns, ["reference_allele", "REF", "A1"], False)
    ac = find_col(chunk.columns, ["alternate_allele", "ALT", "A2"], False)
    nc = find_col(chunk.columns, ["sample_size", "N", "n"], False)
    chunk["chrom"] = chunk[cc].map(norm_chr)
    chunk["pos"] = pd.to_numeric(chunk[pc], errors="coerce")
    sub = chunk.loc[(chunk["chrom"] == CHROM) & (chunk["pos"] >= REGION_START) &
                    (chunk["pos"] <= REGION_END)].copy()
    if not sub.empty:
        sub["p"] = pd.to_numeric(sub[pvc], errors="coerce")
        sub["ref"] = sub[rc].astype(str).str.upper() if rc else "NA"
        sub["alt"] = sub[ac].astype(str).str.upper() if ac else "NA"
        sub["n"] = pd.to_numeric(sub[nc], errors="coerce") if nc else np.nan
        gwas_parts.append(sub)

gwas = pd.concat(gwas_parts, ignore_index=True)
gwas = gwas.loc[gwas["p"].notna() & (gwas["p"] > 0) & (gwas["p"] <= 1)].copy()
gwas["pos"] = gwas["pos"].astype(int)
gwas["mlogp"] = -np.log10(gwas["p"])
gwas = gwas.sort_values("p").drop_duplicates("pos", keep="first").sort_values("pos").reset_index(drop=True)
gwas_lead = gwas.loc[[gwas["p"].idxmin()]].copy()
gwas_target = gwas.loc[gwas["pos"] == EMVAR_POS].copy()
if not gwas_target.empty:
    allele_match = gwas_target.loc[gwas_target.apply(lambda r: {r["ref"], r["alt"]} == {"C", "T"}, axis=1)]
    if not allele_match.empty:
        gwas_target = allele_match
    gwas_target = gwas_target.loc[[gwas_target["p"].idxmin()]].copy()

n_gwas_sig = int((gwas["p"] < GWAS_THRESHOLD).sum())
n_gwas_sugg = int(((gwas["p"] >= GWAS_THRESHOLD) & (gwas["p"] < GWAS_SUGGESTIVE)).sum())
# Use suggestive GWAS sites to filter the eQTL records.
gwas_p1e5_pos = set(gwas.loc[gwas["p"] < GWAS_SUGGESTIVE, "pos"].astype(int))
print(f"  GWAS SNPs in region: {len(gwas):,}  P<5e-8: {n_gwas_sig:,}  5e-8≤P<1e-5: {n_gwas_sugg:,}")

# Load cis-eQTL associations for all genes in the region.
print("Reading cis-eQTLs...")
eqtl_parts = []
emvar_gene_parts = []  # Collect genes associated with the emVar position.
for chunk in pd.read_csv(EQTL_FILE, sep="\t", chunksize=500_000, low_memory=False):
    gene_col = find_col(chunk.columns, ["Gene.ID", "gene_id", "gene"])
    snp_col_tmp = find_col(chunk.columns, ["SNP_id", "snp", "variant_id"])
    chunk[["chrom_tmp", "pos_tmp"]] = chunk[snp_col_tmp].apply(parse_snp)
    sub = chunk.loc[(chunk["chrom_tmp"] == CHROM) & chunk["pos_tmp"].notna() &
                    (chunk["pos_tmp"] >= REGION_START) & (chunk["pos_tmp"] <= REGION_END)].copy()
    if not sub.empty:
        eqtl_parts.append(sub)

    emvar_sub = chunk.loc[(chunk["chrom_tmp"] == CHROM) & (chunk["pos_tmp"] == EMVAR_POS)].copy()
    if not emvar_sub.empty:
        emvar_gene_parts.append(emvar_sub[[gene_col]].copy())

eqtl = pd.concat(eqtl_parts, ignore_index=True)
snp_col = find_col(eqtl.columns, ["SNP_id", "snp", "variant_id"])
gene_col = find_col(eqtl.columns, ["Gene.ID", "gene_id", "gene"])
eqtl["chrom"] = eqtl["chrom_tmp"]
eqtl["pos"] = pd.to_numeric(eqtl["pos_tmp"], errors="coerce")
eqtl["gene_id_std"] = eqtl[gene_col].astype(str)
eqtl["p"] = pd.to_numeric(eqtl[find_col(eqtl.columns, ["p.value", "p_value", "pvalue", "P"])], errors="coerce")
eqtl["fdr"] = pd.to_numeric(eqtl[find_col(eqtl.columns, ["BBFDR", "FDR", "qvalue"])], errors="coerce")
eqtl["beta_std"] = pd.to_numeric(eqtl[find_col(eqtl.columns, ["beta", "effect", "slope"])], errors="coerce")
eqtl = eqtl.loc[eqtl["pos"].notna() & eqtl["p"].notna() &
                (eqtl["p"] > 0) & (eqtl["p"] <= 1)].copy()
require_nonempty(eqtl, "regional cis-eQTL")
eqtl["pos"] = eqtl["pos"].astype(int)

# Keep the target gene eQTL set for reference.
target_eqtl_positions = set(eqtl.loc[eqtl["gene_id_std"] == TARGET_GENE, "pos"].astype(int).tolist())

# Collect all gene IDs associated with the emVar position.
if emvar_gene_parts:
    emvar_gene_df = pd.concat(emvar_gene_parts, ignore_index=True)
    emvar_gene_col = emvar_gene_df.columns[0]
    emvar_gene_ids = set(emvar_gene_df[emvar_gene_col].astype(str).unique())
else:
    emvar_gene_ids = {TARGET_GENE}
print(f"  emVar-associated genes ({len(emvar_gene_ids)}): {sorted(emvar_gene_ids)}")

# Keep eQTLs at GWAS sites with P < 1e-5.
eqtl_plot = eqtl.loc[eqtl["pos"].isin(gwas_p1e5_pos)].copy()
require_nonempty(eqtl_plot, "eQTLs at GWAS positions with P < 1e-5")
eqtl_plot["mlogp"] = -np.log10(eqtl_plot["p"])
# Collect targets of the displayed eQTLs.
b_target_genes = set(eqtl_plot["gene_id_std"].unique())
print(f"  Panel C cis-eQTLs: {len(eqtl_plot):,} rows / {eqtl_plot['pos'].nunique():,} sites / {len(b_target_genes)} target genes")

# Read LD values for displayed eQTL positions.
plot_positions = set(eqtl_plot["pos"].astype(int).tolist())
ld_b = pd.read_csv(LD_FILE_B, sep=r"\s+")
ld_b.columns = [str(c).strip() for c in ld_b.columns]
ld_b["CHR_B"] = ld_b["CHR_B"].map(norm_chr)
ld_b["BP_B"] = pd.to_numeric(ld_b["BP_B"], errors="coerce")
ld_b["R2"] = pd.to_numeric(ld_b["R2"], errors="coerce")
ld_filtered_b = ld_b.loc[(ld_b["CHR_B"] == CHROM) & ld_b["BP_B"].isin(plot_positions) & ld_b["R2"].notna()].copy()
ld_filtered_b["BP_B"] = ld_filtered_b["BP_B"].astype(int)
ld_filtered_b = ld_filtered_b.sort_values("R2", ascending=False).drop_duplicates("BP_B", keep="first").sort_values("BP_B")

eqtl_plot = eqtl_plot.merge(ld_filtered_b[["BP_B", "R2"]], left_on="pos", right_on="BP_B", how="left")
eqtl_lead = eqtl_plot.loc[[eqtl_plot["p"].idxmin()]].copy()
eqtl_target = eqtl_plot.loc[eqtl_plot["pos"] == EMVAR_POS].copy()
if not eqtl_target.empty:
    # Prefer the target gene at the emVar position, then the smallest P value.
    focal = eqtl_target.loc[eqtl_target["gene_id_std"] == TARGET_GENE]
    eqtl_target = (focal if not focal.empty else eqtl_target)
    eqtl_target = eqtl_target.loc[[eqtl_target["p"].idxmin()]].copy()
n_emvar_genes_at_pos = int(eqtl_plot.loc[eqtl_plot["pos"] == EMVAR_POS, "gene_id_std"].nunique())

# Intersect the MPRA and displayed eQTL positions.
print("Reading MPRA results...")
mpra_all = pd.read_csv(MPRA_RESULT_FILE, sep="\t", low_memory=False)
mpra_all.columns = [str(c).strip() for c in mpra_all.columns]
chr_col = find_col(mpra_all.columns, ["chr", "chrom", "chromosome"])
pos_col = find_col(mpra_all.columns, ["variant_pos", "POS", "position"])
mpra_all["chr_std"] = mpra_all[chr_col].map(norm_chr)
mpra_all["pos_std"] = pd.to_numeric(mpra_all[pos_col], errors="coerce")
mpra_all["pk15_q"] = pd.to_numeric(mpra_all["PK15.q_BH_all"], errors="coerce").fillna(1)
mpra_all["psc_q"] = pd.to_numeric(mpra_all["PSC.q_BH_all"], errors="coerce").fillna(1)
mpra_all = mpra_all.loc[(mpra_all["chr_std"] == CHROM) & mpra_all["pos_std"].isin(plot_positions)].copy()
mpra_all["pos_std"] = mpra_all["pos_std"].astype(int)

pk15_best = mpra_all.loc[mpra_all.groupby("pos_std")["pk15_q"].idxmin()].copy()
psc_best = mpra_all.loc[mpra_all.groupby("pos_std")["psc_q"].idxmin()].copy()
mpra_merged = pd.merge(pk15_best[["pos_std", "chr_std", "pk15_q", "PK15.emVar_all_qc"]],
                        psc_best[["pos_std", "psc_q", "PSC.emVar_all_qc"]],
                        on="pos_std", how="outer").rename(columns={"pos_std": "variant_pos"})
mpra_merged["pk15_q"] = pd.to_numeric(mpra_merged["pk15_q"], errors="coerce").fillna(1)
mpra_merged["psc_q"] = pd.to_numeric(mpra_merged["psc_q"], errors="coerce").fillna(1)
mpra_merged["min_q"] = mpra_merged[["pk15_q", "psc_q"]].min(axis=1)
mpra_merged["min_q"] = mpra_merged["min_q"].clip(lower=np.nextafter(0, 1))
mpra_merged["neg_log10_q"] = -np.log10(mpra_merged["min_q"])
mpra_merged["is_emvar"] = as_bool(mpra_merged["PK15.emVar_all_qc"]) | as_bool(mpra_merged["PSC.emVar_all_qc"])

# Classify plotted MPRA positions by emVar status.
mpra_merged["evidence"] = np.where(mpra_merged["is_emvar"], "emVar", "eQTL")
print("  MPRA evidence counts:", mpra_merged["evidence"].value_counts().to_dict())

# Read gene annotations.
print("Reading GTF annotations...")
gtf_cols = ["seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"]
gtf_all = pd.read_csv(GTF_FILE, sep="\t", comment="#", names=gtf_cols,
                      usecols=["seqname", "feature", "start", "end", "strand", "attribute"],
                      low_memory=False)
gtf_all["chrom"] = gtf_all["seqname"].map(norm_chr)
gtf_all["gene_id"] = gtf_all["attribute"].map(lambda x: gtf_attr(x, "gene_id"))
gtf_all["gene_name"] = gtf_all["attribute"].map(lambda x: gtf_attr(x, "gene_name"))
genes = gtf_all.loc[(gtf_all["feature"] == "gene") & (gtf_all["chrom"] == CHROM) &
                    (gtf_all["end"] >= REGION_START) & (gtf_all["start"] <= REGION_END)].copy()
genes["display"] = genes["gene_name"].fillna(genes["gene_id"])
genes = genes.drop_duplicates(["gene_id", "start", "end", "strand"]).sort_values("start").reset_index(drop=True)

exons_raw = gtf_all.loc[(gtf_all["feature"] == "exon") & (gtf_all["chrom"] == CHROM) &
                        (gtf_all["end"] >= REGION_START) & (gtf_all["start"] <= REGION_END)].copy()
exons_raw["gene_id"] = exons_raw["attribute"].map(lambda x: gtf_attr(x, "gene_id"))
exons_raw["transcript_id"] = exons_raw["attribute"].map(lambda x: gtf_attr(x, "transcript_id"))
exons_raw = exons_raw.dropna(subset=["gene_id", "transcript_id"])
tx_len = exons_raw.groupby(["gene_id", "transcript_id"]).agg(
    tx_start=("start", "min"), tx_end=("end", "max")).reset_index()
tx_len["tx_len"] = tx_len["tx_end"] - tx_len["tx_start"]
best_tx = tx_len.loc[tx_len.groupby("gene_id")["tx_len"].idxmax()][["gene_id", "transcript_id"]]
exons = exons_raw.merge(best_tx, on=["gene_id", "transcript_id"])
exons = exons.groupby(["gene_id", "start", "end"]).first().reset_index()

# Color genes according to evidence at the displayed eQTL positions.
label_genes = set(b_target_genes) | set(emvar_gene_ids)
in_region_ids = set(genes["gene_id"].dropna())
n_label_in_region = int(genes["gene_id"].isin(label_genes).sum())
n_emvar_in_region = int(genes["gene_id"].isin(emvar_gene_ids).sum())
n_eqtl_only_in_region = int(len((b_target_genes - emvar_gene_ids) & in_region_ids))
print(f"  Genes labeled: {n_label_in_region} ({n_emvar_in_region} emVar targets and {n_eqtl_only_in_region} eQTL-only targets)")
print(f"  Target genes in the ±500 kb window: {len(b_target_genes & in_region_ids)} / {len(b_target_genes)}")

# Load the five chromatin tracks.
print("Reading chromatin signals...")
track_names = ["ATAC", "CTCF", "H3K27ac", "H3K4me1", "H3K4me3"]
track_signals = {}
for name in track_names:
    print(f"  {name}...")
    bins, scores = load_bw_signal(BW_FILES[name], CHROM_BW, REGION_START, REGION_END, bin_size=200)
    track_signals[name] = (bins, scores)

# Read base-level signal in the three 5-kb flanks for the inset tracks.
ZOOM_FLANK = 5_000
ZOOM_BIN = 50
ZOOM_AXIS_Y = -0.17        # Align zoom boxes with the MPRA x-axis.
ZOOM_BOX_H = 0.16
GWAS_LEAD_POS = int(gwas_lead.iloc[0]["pos"])   # GWAS lead position.
ZOOM_SITES = [
    ("eQTL lead", EQTL_LEAD_POS, LEAD_BLUE),
    ("GWAS lead", GWAS_LEAD_POS, "#7B3294"),
    ("emVar / FRV", EMVAR_POS, ORANGE),
]
print("Reading unsmoothed per-base signals for five tracks in the ±5 kb insets...")
zoom_signals = {}
for zname, zpos, _ in ZOOM_SITES:
    zs, ze = zpos - ZOOM_FLANK, zpos + ZOOM_FLANK
    for tn in track_names:
        zb, zsc = load_bw_raw(BW_FILES[tn], CHROM_BW, zs, ze)
        zoom_signals[(zname, tn)] = (zb, zsc)
    print(f"  {zname} @ {zpos:,}")

print("Inputs loaded.")

# 2. Plot the five panels





fig = plt.figure(figsize=(11.0, 9.8))
gs = fig.add_gridspec(5, 3, height_ratios=[2.0, 3.0, 2.0, 2.0, 3.0],
                      left=0.125, right=0.985, top=0.965, bottom=0.055, hspace=0.30)

ax_gwas = fig.add_subplot(gs[0, :])
ax_gene = fig.add_subplot(gs[1, :], sharex=ax_gwas)
ax_eqtl = fig.add_subplot(gs[2, :], sharex=ax_gwas)
ax_mpra = fig.add_subplot(gs[3, :], sharex=ax_gwas)

# Nest five chromatin tracks for each of three positions.
gs_zoom = gs[4, :].subgridspec(5, 3, hspace=0.28, wspace=0.12)
ax_zoom = [[fig.add_subplot(gs_zoom[r, c]) for r in range(5)] for c in range(3)]
ax_tracks = []

# A. Regional GWAS
sig = gwas.loc[gwas["p"] < GWAS_THRESHOLD]
sugg = gwas.loc[(gwas["p"] >= GWAS_THRESHOLD) & (gwas["p"] < GWAS_SUGGESTIVE)]
bg = gwas.loc[gwas["p"] >= GWAS_SUGGESTIVE]
ax_gwas.scatter(bg["pos"], bg["mlogp"], c=GREY, s=10, alpha=0.6, linewidth=0)
if not sugg.empty:
    ax_gwas.scatter(sugg["pos"], sugg["mlogp"], c=SUGGESTIVE, s=13, edgecolor="white",
                    linewidth=0.15, alpha=0.85, zorder=3)
if not sig.empty:
    ax_gwas.scatter(sig["pos"], sig["mlogp"], c=RED, s=18, edgecolor="white", linewidth=0.2,
                    zorder=4)

# Highlight the focal emVar.
if not gwas_target.empty:
    row = gwas_target.iloc[0]
    ax_gwas.scatter(row["pos"], row["mlogp"], marker="D", s=18, c=ORANGE, edgecolor="black", linewidth=0.8, zorder=7)
    ntxt = f"; n={int(row['n']):,}" if pd.notna(row["n"]) else ""
    ax_gwas.annotate(f"{EMVAR_LABEL}\nP={row['p']:.2e}{ntxt}",
                     (EMVAR_POS, row["mlogp"]), xytext=(14, 18), textcoords="offset points",
                     fontsize=8, color="#9C6500",
                     arrowprops={"arrowstyle": "-", "color": ORANGE, "lw": 0.9},
                     bbox=dict(boxstyle="round,pad=.15", fc="white", ec="none", alpha=0.82))

# Label the lead GWAS variant and its alleles.
gl = gwas_lead.iloc[0]
GWAS_LEAD_POS = int(gl["pos"])
ax_gwas.scatter(gl["pos"], gl["mlogp"], marker="D", s=18, c="#7B3294", edgecolor="black", linewidth=0.6, zorder=6)
if gl["pos"] != EMVAR_POS:
    ax_gwas.annotate(f"Regional lead\n{CHROM}:{gl['pos']:,} {gl['ref']}>{gl['alt']} P={gl['p']:.2e}",
                     (gl["pos"], gl["mlogp"]), xytext=(40, 10), textcoords="offset points",
                     fontsize=7.5, color="#4A1486",
                     arrowprops={"arrowstyle": "-", "color": "#7B3294", "lw": 0.9},
                     bbox=dict(boxstyle="round,pad=.15", fc="white", ec="none", alpha=0.82))

ax_gwas.axhline(-np.log10(GWAS_THRESHOLD), ls="--", lw=0.8, c="#555555")
ax_gwas.axhline(-np.log10(GWAS_SUGGESTIVE), ls=":", lw=0.8, c="#888888")
ax_gwas.text(0.01, 0.94, f"A  PigBiobank {TRAIT} GWAS", transform=ax_gwas.transAxes, va="top", fontweight="bold")
ax_gwas.text(0.99, 0.77,
             f"Regional SNPs: {len(gwas):,} | P<5×10⁻⁸: {len(sig):,} | 5×10⁻⁸≤P<10⁻⁵: {len(sugg):,}",
             transform=ax_gwas.transAxes, ha="right", va="top", fontsize=8,
             bbox=dict(boxstyle="round,pad=.20", fc="white", ec="none", alpha=0.88))
ax_gwas.set_ylabel(r"$-\log_{10}(P)$")
ax_gwas.legend(handles=[
    Line2D([0], [0], marker="D", color="none", markerfacecolor="#7B3294", markeredgecolor="black", markersize=6, label="Regional GWAS lead"),
    Line2D([0], [0], marker="D", color="none", markerfacecolor=ORANGE, markeredgecolor="black", markersize=6, label="Target MPRA emVar"),
    Line2D([0], [0], marker="o", color="none", markerfacecolor=RED, markeredgecolor="white", markersize=6, label=r"GWAS $P<5\times10^{-8}$"),
    Line2D([0], [0], marker="o", color="none", markerfacecolor=SUGGESTIVE, markeredgecolor="white", markersize=6, label=r"$5\times10^{-8}\leq P<10^{-5}$"),
], frameon=True, facecolor="white", edgecolor="none", framealpha=0.90, loc="lower right", fontsize=7.3)

# B. Gene structure
genes_plot = genes.sort_values(["start", "end"]).copy()
lane_ends = []
gene_lanes = {}
label_buffer_bp = 80_000
for idx, row in genes_plot.iterrows():
    visible_start = max(int(row["start"]), REGION_START)
    visible_end = min(int(row["end"]), REGION_END)
    assigned_lane = None
    for lane_idx, last_end in enumerate(lane_ends):
        if visible_start > last_end + label_buffer_bp:
            assigned_lane = lane_idx
            lane_ends[lane_idx] = visible_end
            break
    if assigned_lane is None:
        assigned_lane = len(lane_ends)
        lane_ends.append(visible_end)
    gene_lanes[idx] = assigned_lane

n_lanes = max(len(lane_ends), 1)

ax_gene.add_patch(plt.Rectangle((REGION_START, 0), REGION_END - REGION_START, n_lanes + 1.2,
                                facecolor='white', edgecolor='none', zorder=0))

for idx, row in genes_plot.iterrows():
    gid = row["gene_id"]
    y = n_lanes - gene_lanes[idx]
    start = int(row["start"])
    end = int(row["end"])
    strand = row["strand"]

    is_emvar_gene = gid in emvar_gene_ids
    is_eqtl_target = gid in b_target_genes
    is_labelled = is_emvar_gene or is_eqtl_target

    # Use evidence-dependent styling for gene structures.
    if is_emvar_gene:
        color = ORANGE
    elif is_eqtl_target:
        color = LEAD_BLUE
    else:
        color = "#666666"
    width = 2.8 if is_labelled else 1.3

    ax_gene.plot([start, end], [y, y], color=color, linewidth=width, solid_capstyle="butt", zorder=2, clip_on=True)
    for _, exon in exons.loc[exons["gene_id"] == gid].iterrows():
        exon_start = max(int(exon["start"]), REGION_START)
        exon_end = min(int(exon["end"]), REGION_END)
        if exon_end >= exon_start:
            ax_gene.add_patch(plt.Rectangle((exon_start, y - 0.18), exon_end - exon_start + 1, 0.36,
                                            facecolor=color, edgecolor="none", alpha=0.35, zorder=3, clip_on=True))


    midpoint = (max(start, REGION_START) + min(end, REGION_END)) / 2
    # Use gene names as labels.
    ax_gene.text(midpoint, y + 0.15, str(row["display"]),
                 ha="center", va="bottom",
                 fontsize=7.0 if is_labelled else 6.3,
                 fontstyle="italic",
                 fontweight="bold" if is_labelled else "normal",
                 color=color, clip_on=True)

ax_gene.set_ylim(-0.3, n_lanes + 1.5)
ax_gene.set_yticks([])
ax_gene.text(0.01, 0.94,
             f"B  Genes within ±500 kb (n={len(genes):,}; all labelled; coloured = eQTL/FRV targets, n={n_label_in_region})",
             transform=ax_gene.transAxes, va="top", fontweight="bold")
ax_gene.legend(handles=[
    Line2D([0], [0], color=ORANGE, lw=2.8, label=f"FRV 7:23,797,335 target gene (n={n_emvar_in_region})"),
    Line2D([0], [0], color=LEAD_BLUE, lw=2.2, label=f"GWAS-derived eQTL target gene (n={n_eqtl_only_in_region})"),
    Line2D([0], [0], color="#666666", lw=1.3, label=f"other genes (n={len(genes) - n_label_in_region})"),
], loc="lower right", fontsize=6.8, frameon=True, facecolor="white", edgecolor="none",
    framealpha=0.90, handlelength=1.6)
# Show genomic coordinates below the gene panel.
ax_gene.set_xlabel("Chromosome 7 position (Mb)")

# C. Cis-eQTL associations and LD
ax_eqtl.scatter(eqtl_plot["pos"], eqtl_plot["mlogp"],
                c=[ld_color(x) for x in eqtl_plot["R2"]],
                s=16, alpha=0.88, linewidth=0)
el = eqtl_lead.iloc[0]
ax_eqtl.scatter(el["pos"], el["mlogp"], marker="D", s=16, c=LEAD_BLUE, edgecolor="black", linewidth=0.65, zorder=6)
ax_eqtl.annotate(f"Lead cis-eQTL\n{el[snp_col]}\n{el[gene_col]}  P={el['p']:.2e}",
                 (el["pos"], el["mlogp"]), xytext=(-12, -8), textcoords="offset points",
                 fontsize=7.5, color=LEAD_BLUE, ha="right", va="top",
                 bbox=dict(boxstyle="round,pad=.15", fc="white", ec="none", alpha=0.80))
if not eqtl_target.empty:
    row = eqtl_target.iloc[0]
    target_ld = row["R2"]
    ld_text = f"{target_ld:.3f}" if pd.notna(target_ld) else "NA"
    ax_eqtl.scatter(row["pos"], row["mlogp"], marker="D", s=16, c=ORANGE, edgecolor="black", linewidth=0.8, zorder=7)
    genes_txt = f"; {n_emvar_genes_at_pos} genes" if n_emvar_genes_at_pos > 1 else ""
    ax_eqtl.annotate(f"{EMVAR_LABEL}{genes_txt}\n{row[gene_col]}  P={row['p']:.2e}; β={row['beta_std']:.3f}\n"
                     f"BBFDR={row['fdr']:.2e}; r²={ld_text}",
                     (EMVAR_POS, row["mlogp"]), xytext=(18, 20), textcoords="offset points",
                     fontsize=7.5, color="#9C6500",
                     arrowprops={"arrowstyle": "-", "color": ORANGE, "lw": 0.9},
                     bbox=dict(boxstyle="round,pad=.15", fc="white", ec="none", alpha=0.82))
ax_eqtl.text(0.01, 0.94,
             f"C  DLC muscle cis-eQTL (GWAS P<10⁻⁵ sites)",
             transform=ax_eqtl.transAxes, va="top", fontweight="bold")
ax_eqtl.text(0.99, 0.94,
             f"SNP-gene pairs: {len(eqtl_plot):,} | SNPs: {eqtl_plot['pos'].nunique():,} | target genes: {len(b_target_genes)}",
             transform=ax_eqtl.transAxes, ha="right", va="top", fontsize=7.5)
ax_eqtl.set_ylabel(r"$-\log_{10}(P)$")

# Draw the LD legend as a compact row of color swatches.
cax = ax_eqtl.inset_axes([0.882, 0.085, 0.098, 0.050])
cax.set_xlim(0, 1)
cax.set_ylim(0, 1)
cax.axis("off")
cax.add_patch(plt.Rectangle((0, 0), 1, 1, facecolor="white", edgecolor="#DDDDDD",
                            linewidth=0.5, alpha=0.94, zorder=0))
# Use gray for missing LD and five colors for available values.
SQ_Y, SQ_H, SQ_W = 0.42, 0.52, 0.135
cax.add_patch(plt.Rectangle((0.04, SQ_Y), SQ_W * 0.75, SQ_H,
                            facecolor=GREY, edgecolor="white", linewidth=0.4))
for i, col in enumerate(LD_COLORS):
    cax.add_patch(plt.Rectangle((0.04 + (i + 0.85) * SQ_W, SQ_Y), SQ_W, SQ_H,
                                facecolor=col, edgecolor="white", linewidth=0.4))
# Label the ends of the color scale.
cax.text(0.04 + 0.85 * SQ_W, 0.30, "0", fontsize=3.2, va="top", ha="left", color="#333333")
cax.text(0.04 + 5.85 * SQ_W, 0.30, "1", fontsize=3.2, va="top", ha="right", color="#333333")
# Label the LD measure as r².
cax.text(0.04 + 0.38 * SQ_W, 0.28, r"$r^2$", fontsize=3.0, va="top", ha="center", color="#666666")

# D. MPRA allelic-skew significance
EV_COLOR = {"eQTL": EQTL_ONLY_BLUE, "emVar": ORANGE}
EV_LABEL = {"eQTL": "eQTL (C panel)", "emVar": "MPRA emVar"}
for ev in ["eQTL", "emVar"]:
    sub = mpra_merged.loc[mpra_merged["evidence"] == ev]
    if sub.empty:
        continue
    size = 38 if ev == "emVar" else 28
    ax_mpra.scatter(sub["variant_pos"], sub["neg_log10_q"], s=size, c=EV_COLOR[ev],
                    marker="D" if ev == "emVar" else "o",
                    edgecolor="black" if ev == "emVar" else "#636363",
                    linewidth=0.65 if ev == "emVar" else 0.5, alpha=0.92, zorder=3 if ev == "emVar" else 2,
                    label=EV_LABEL[ev])
focal = mpra_merged.loc[mpra_merged["variant_pos"] == EMVAR_POS]
if not focal.empty:
    row = focal.iloc[0]
    ax_mpra.annotate(EMVAR_LABEL, xy=(EMVAR_POS, row["neg_log10_q"]), xytext=(12, 12),
                     textcoords="offset points", fontsize=8, color="#9C6500",
                     arrowprops={"arrowstyle": "-", "color": ORANGE, "lw": 0.8},
                     bbox=dict(boxstyle="round,pad=.15", fc="white", ec="none", alpha=0.84))
ax_mpra.axhline(-np.log10(0.1), ls="--", lw=0.8, c="#666666")
ax_mpra.text(0.01, 0.94, "D  MPRA $\\cap$ cis-eQTL sites (panel C)", transform=ax_mpra.transAxes, va="top", fontweight="bold")
ax_mpra.text(0.01, 0.76,
             f"cis-eqtl(C) $\\cap$ MPRA: {len(mpra_merged):,} | emVar: {int((mpra_merged['evidence']=='emVar').sum()):,}",
             transform=ax_mpra.transAxes, ha="left", va="top", fontsize=7.5,
             bbox=dict(boxstyle="round,pad=.18", fc="white", ec="none", alpha=0.88))
ax_mpra.set_ylabel(r"$-\log_{10}(\mathrm{padj-BH})$")
ax_mpra.set_ylim(bottom=ZOOM_AXIS_Y, top=2.5)
if not mpra_merged.empty:
    ax_mpra.legend(loc="upper right", frameon=True, facecolor="white", edgecolor="none",
                   framealpha=0.90, fontsize=7, ncol=2)

# E. Base-level chromatin tracks around three positions
# Scale tracks to their maxima and show row labels on the left.
for ci, (zname, zpos, zcol) in enumerate(ZOOM_SITES):
    zs, ze = zpos - ZOOM_FLANK, zpos + ZOOM_FLANK
    for ri, tn in enumerate(track_names):
        ax = ax_zoom[ci][ri]
        bins, scores = zoom_signals[(zname, tn)]
        color = TRACK_COLORS[tn]
        ax.fill_between(bins, scores, color=color, alpha=0.92, linewidth=0, clip_on=True)
        ax.plot(bins, scores, color=color, linewidth=0.25, clip_on=True)
        ymax = float(scores.max()) if scores.max() > 0 else 1.0
        ax.set_ylim(0, ymax * 1.15)
        ax.set_xlim(zs, ze)
        ax.set_yticks([])
        ax.spines[["top", "right"]].set_visible(False)
        # Mark the center of each inset.
        ax.axvline(zpos, color=zcol, ls="--", lw=0.8, alpha=0.6)
        # Share track labels across inset columns.
        if ci == 0:
            ax.text(-0.03, 0.5, tn, transform=ax.transAxes, va="center", ha="right",
                    fontsize=6.5, fontweight="bold", color=color)
        if ri == len(track_names) - 1:
            ax.set_xticks([zs, zpos, ze])
            ax.set_xticklabels([f"{zs/1e6:.4f}", f"{zpos/1e6:.4f}", f"{ze/1e6:.4f}"],
                               fontsize=6)
            ax.set_xlabel("Chr7 position (Mb)", fontsize=6.5)
        else:
            ax.tick_params(axis="x", labelbottom=False)

# Connect MPRA zoom boxes to the corresponding track panels.
for ci, (zname, zpos, zcol) in enumerate(ZOOM_SITES):
    xs, xe = zpos - ZOOM_FLANK, zpos + ZOOM_FLANK
    ax_mpra.add_patch(plt.Rectangle((xs, ZOOM_AXIS_Y), xe - xs, ZOOM_BOX_H,
                                    facecolor=zcol, alpha=0.35, edgecolor=zcol,
                                    linewidth=1.3, zorder=9, clip_on=False))
    for xa, xb in [(xs, 0.0), (xe, 1.0)]:
        con = ConnectionPatch(xyA=(xa, ZOOM_AXIS_Y), coordsA=ax_mpra.transData,
                              xyB=(xb, 1.02), coordsB=ax_zoom[ci][0].transAxes,
                              color=zcol, lw=0.9, alpha=0.8, zorder=6)
        fig.add_artist(con)

# Shared genomic x-axis
all_axes = [ax_gwas, ax_gene, ax_eqtl, ax_mpra] + ax_tracks
for ax in all_axes:
    ax.set_xlim(REGION_START, REGION_END)
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="x", color="#EEEEEE", lw=0.5)
    if REGION_START <= GWAS_LEAD_POS <= REGION_END:
        ax.axvline(GWAS_LEAD_POS, c="#7B3294", ls="--", lw=0.80, alpha=0.58)   # GWAS lead
    if REGION_START <= EMVAR_POS <= REGION_END:
        ax.axvline(EMVAR_POS, c=ORANGE, ls="--", lw=1.35, alpha=0.68)

# Draw the eQTL lead guide only below the gene panel.
for ax in [ax_eqtl, ax_mpra] + ax_tracks:
    if REGION_START <= EQTL_LEAD_POS <= REGION_END:
        ax.axvline(EQTL_LEAD_POS, c=LEAD_BLUE, ls="--", lw=0.80, alpha=0.58)


ticks = np.linspace(REGION_START, REGION_END, 5)
ax_gene.set_xticks(ticks)
ax_gene.set_xticklabels([f"{x / 1e6:.3f}" for x in ticks])
for ax in [ax_gwas, ax_eqtl, ax_mpra] + ax_tracks:
    ax.tick_params(axis="x", labelbottom=False)

fig.suptitle(f"Multi-layer evidence across chr7:23,797,335 ±500 kb",
             x=0.125, y=0.988, ha="left", fontsize=14, fontweight="bold")

# 3. Export figures
prefix = WORK_DIR / "chr7_23797335_500kb_5panel_ABCDE_update3"
fig.savefig(f"{prefix}.png", dpi=1200, bbox_inches="tight", facecolor="white")
fig.savefig(f"{prefix}.pdf", bbox_inches="tight", transparent=True)
plt.close(fig)

print("=" * 70)
print(f"Output: {prefix}.png / .pdf")
print(f"  GWAS   : {len(gwas):,} SNP  (P<5e-8: {len(sig):,} / 5e-8~1e-5: {len(sugg):,})")
print(f"  eQTL   : {len(eqtl_plot):,} rows / {eqtl_plot['pos'].nunique():,} sites / {len(b_target_genes)} target genes")
print(f"  MPRA   : {len(mpra_merged):,} points  {mpra_merged['evidence'].value_counts().to_dict()}")
print(f"  Genes  : {len(genes):,} total, {n_label_in_region} labeled ({n_emvar_in_region} emVar targets)")
print("=" * 70)


# 4. Summarize filtering and export Source Data

print("\n" + "=" * 70)
print("Filtering summary")
print("=" * 70)


n_eqtl_region_rows = len(eqtl)
n_eqtl_region_genes = eqtl["gene_id_std"].nunique()
n_genes_total = len(genes)
in_region_ids2 = set(genes["gene_id"].dropna())
all_target_ids = set(eqtl["gene_id_std"].unique())
genes_with_eqtl = in_region_ids2 & all_target_ids
genes_coloured = in_region_ids2 & (b_target_genes | emvar_gene_ids)
genes_out_of_window = b_target_genes - in_region_ids2
genes_frv = in_region_ids2 & emvar_gene_ids
genes_eqtl_only = in_region_ids2 & (b_target_genes - emvar_gene_ids)

# Summarize the absolute SNP-to-TSS distance.
if {"Start", "End"}.issubset(eqtl.columns):
    _tss = np.where(eqtl["Start"] < eqtl["End"], eqtl["Start"], eqtl["End"])
    d_tss = (eqtl["pos"] - _tss).abs()
else:
    d_tss = pd.Series([np.nan])
d_med, d_p99, d_max = d_tss.median(), d_tss.quantile(0.99), d_tss.max()

# Count all regional MPRA positions before the panel C filter.
_m = pd.read_csv(MPRA_RESULT_FILE, sep="\t", low_memory=False)
_m.columns = [str(c).strip() for c in _m.columns]
_mc = find_col(_m.columns, ["chr", "chrom", "chromosome"])
_mp = find_col(_m.columns, ["variant_pos", "POS", "position"])
_m = _m.loc[(_m[_mc].map(norm_chr) == CHROM) &
            (pd.to_numeric(_m[_mp], errors="coerce").between(REGION_START, REGION_END))]
n_mpra_region = pd.to_numeric(_m[_mp], errors="coerce").astype(int).nunique()

flow_rows = [
    ("A", "GWAS", "Regional SNPs",                       len(gwas),                    f"chr7:{REGION_START:,}-{REGION_END:,}"),
    ("A", "GWAS", "  P < 5×10⁻⁸",                     n_gwas_sig,                   "Genome-wide significant"),
    ("A", "GWAS", "  5×10⁻⁸ ≤ P < 10⁻⁵",              n_gwas_sugg,                  "Suggestive, highlighted separately"),
    ("A", "GWAS", "  Total P < 10⁻⁵",                  len(gwas_p1e5_pos),           "Positions used to filter panel C"),
    ("C", "eQTL", "Regional cis-eQTL associations before GWAS filtering", n_eqtl_region_rows,           f"{n_eqtl_region_genes} target genes"),
    ("C", "eQTL", "  Associations at GWAS sites with P < 10⁻⁵",  len(eqtl_plot),               f"{len(b_target_genes)} target genes"),
    ("C", "eQTL", "  Unique SNP positions",                 int(eqtl_plot['pos'].nunique()), "One SNP may have multiple target genes"),
    ("B", "Gene", "Regional genes",                        n_genes_total,                "GTF ±500kb"),
    ("B", "Gene", "  Genes with significant cis-eQTLs",             len(genes_with_eqtl),         "eQTL table BBFDR < 0.05"),
    ("B", "Gene", "  Genes with eQTLs at GWAS sites with P < 10⁻⁵", len(genes_coloured),          "Colored in the plot"),
    ("B", "Gene", "     Targets of chr7:23,797,335",      len(genes_frv),               "Orange"),
    ("B", "Gene", "     GWAS-associated eQTL target genes",     len(genes_eqtl_only),         "Blue"),
    ("B", "Gene", "  Regional genes without a significant eQTL",             n_genes_total - len(genes_with_eqtl), "Absent from the eQTL table"),
    ("B", "Gene", "  Panel C target genes outside the window",              len(genes_out_of_window),     "The cis window extends beyond the ±500 kb plot"),
    ("D", "MPRA", "Regional MPRA positions",                  n_mpra_region,                ""),
    ("D", "MPRA", "  Overlap with panel C eQTL positions",              len(mpra_merged),             "Displayed in the plot"),
    ("D", "MPRA", "     emVars",                   int((mpra_merged['evidence'] == 'emVar').sum()), "Orange"),
]
flow = pd.DataFrame(flow_rows, columns=["Panel", "Layer", "Step", "N", "Note"])

# Print the filtering summary with aligned columns.
_w = max(len(s) for s in flow["Step"])
print(f"{'Panel':<6}{'Layer':<7}{'Step':<{_w}}  {'N':>8}  Note")
print("-" * (6 + 7 + _w + 12 + 30))
for _, r in flow.iterrows():
    print(f"{r['Panel']:<6}{r['Layer']:<7}{r['Step']:<{_w}}  {r['N']:>8,}  {r['Note']}")
print("-" * (6 + 7 + _w + 12 + 30))
print(f"Cis |SNP−TSS| distance: median {d_med:,.0f} bp / 99% {d_p99:,.0f} bp / maximum {d_max:,.0f} bp")
print(f"Target genes: {len(b_target_genes & in_region_ids2)} inside and {len(genes_out_of_window)} outside the window")

# ---- Source data ----
gwas_src = gwas.copy()
gwas_src["category"] = np.where(gwas_src["p"] < GWAS_THRESHOLD, "P<5e-8",
                                np.where(gwas_src["p"] < GWAS_SUGGESTIVE, "5e-8~1e-5", "nonsig"))
gwas_src = gwas_src[["pos", "ref", "alt", "p", "mlogp", "n", "category"]].sort_values("pos")

genes_src = genes.copy()
genes_src["is_FRV_target"] = genes_src["gene_id"].isin(emvar_gene_ids)
genes_src["is_eQTL_target"] = genes_src["gene_id"].isin(b_target_genes)
genes_src["colour"] = np.where(genes_src["is_FRV_target"], "orange",
                               np.where(genes_src["is_eQTL_target"], "blue", "grey"))
genes_src["labelled"] = genes_src["is_FRV_target"] | genes_src["is_eQTL_target"]
genes_src = genes_src[["gene_id", "gene_name", "display", "start", "end", "strand",
                       "is_FRV_target", "is_eQTL_target", "colour", "labelled"]].sort_values("start")

eqtl_src = eqtl_plot[[snp_col, gene_col, "pos", "p", "fdr", "beta_std", "R2"]].copy()
eqtl_src.columns = ["SNP_id", "Gene.ID", "pos", "p.value", "BBFDR", "beta", "R2"]
eqtl_src = eqtl_src.sort_values(["pos", "p.value"])

mpra_src = mpra_merged.copy()
mpra_src = mpra_src[["variant_pos", "pk15_q", "psc_q", "min_q", "neg_log10_q", "is_emvar", "evidence"]]

ld_src = ld_filtered_b[["BP_B", "R2"]].copy()
ld_src.columns = ["pos", "R2"]

# Export the per-base values for each of the three zoom regions.
e_rows = []
e_info_rows = []
for zname, zpos, zcol in ZOOM_SITES:
    for tn in track_names:
        zb, zsc = zoom_signals[(zname, tn)]
        e_rows.append(pd.DataFrame({
            "Site": zname, "Track": tn, "Position": zb, "Value": zsc,
        }))
        ymax_used = float(zsc.max()) if zsc.max() > 0 else 1.0
        e_info_rows.append({
            "Site": zname, "Track": tn, "Centre_pos": zpos,
            "Window_start": zpos - ZOOM_FLANK, "Window_end": zpos + ZOOM_FLANK,
            "N_bp": len(zsc), "Max": ymax_used, "Ylim_top": ymax_used * 1.15,
            "Mean": float(np.mean(zsc)), "Colour": TRACK_COLORS[tn],
        })
e_src = pd.concat(e_rows, ignore_index=True)
e_info = pd.DataFrame(e_info_rows)

readme = pd.DataFrame({"Field": [
    "Figure", "Region", "Panel_A", "Panel_B", "Panel_C", "Panel_D", "Panel_E",
    "Filtering_flow", "Test_GWAS", "Test_eQTL", "MPRA_metric",
    "FRV_definition", "Colour_B", "LD_legend", "LD_reference", "LD_file",
    "Input_GWAS", "Input_eQTL", "Input_MPRA", "Input_GTF", "Input_bigwig",
], "Description": [
    "chr7_23797335_500kb_5panel_ABCDE_update3.pdf",
    f"chr7:{REGION_START:,}-{REGION_END:,} (chr7:23,797,335 ±500 kb)",
    "PigBiobank M_LMDEP GWAS; 5e-8 red, 5e-8~1e-5 light orange, else grey; 1e-5 and 5e-8 thresholds",
    f"All {n_genes_total} genes in window, all labelled; FRV targets orange, GWAS-derived eQTL targets blue",
    "cis-eQTL at GWAS P<1e-5 sites (all target genes); points coloured by LD r2 with lead cis-eQTL 7:23,599,787",
    "MPRA sites intersecting panel-C eQTL positions; eQTL vs MPRA emVar (no ASE)",
    "Three ±5 kb zoom panels (E1 lead 7:23,599,787 / E2 GWAS lead 7:23,619,470 / E3 emVar 7:23,797,335), "
    "each with ATAC / CTCF / H3K27ac / H3K4me1 / H3K4me3 from Muscle P350 RPKM bigwigs; "
    "per-base values (no binning/smoothing), ylim = max*1.15 so no peak is clipped; "
    "track names shared across the 3 columns, shown once at far left",
    "See sheet Filtering_flow for the full top-to-bottom numeric chain",
    "PigBiobank summary statistics, genome-wide threshold 5e-8",
    "eQTL_MAF005_kept.tsv pre-filtered to BBFDR<0.05; cis window ~1 Mb",
    "-log10(min PK15.q_BH_all, PSC.q_BH_all)",
    "FRV = functional regulatory variant = emVar (PK15.emVar_all_qc or PSC.emVar_all_qc)",
    "orange = FRV 7:23,797,335 target gene; blue = other GWAS-derived eQTL target gene; grey = other genes",
    "5 discrete bins (0-0.2-0.4-0.6-0.8-1.0) shown as stacked squares; grey = LD unavailable",
    "7:23,599,787 A>G — panel-C lead cis-eQTL after the GWAS P<1e-5 filter (plink --ld-snp 7:23599787:A:G)",
    "DLC_LD_chr7_23599787_lead_eQTL.ld — DLC population, n=189 (same cohort as the eQTL), plink 1.9 --r2, --ld-window-kb 1000, --ld-window-r2 0",
    "M_LMDEP.txt.gz",
    "eQTL_MAF005_kept.tsv",
    "RESULT_CROSSCELL_M3_cre_A0_F_eqtlase_filtered.tsv",
    "Sus_scrofa.Sscrofa11.1.109.gtf",
    "pig_NC_paper2 Muscle P350 RPKM bigwigs",
]})

source_base = WORK_DIR / "chr7_23797335_500kb_5panel_ABCDE_update3_SourceData"
sheets = {
    "README": readme,
    "Filtering_flow": flow,
    "A_GWAS": gwas_src,
    "B_Genes": genes_src,
    "C_eQTL": eqtl_src,
    "C_LD": ld_src,
    "D_MPRA": mpra_src,
    "E_tracks": e_src,
    "E_info": e_info,
}
try:
    with pd.ExcelWriter(f"{source_base}.xlsx", engine="openpyxl") as w:
        for nm, df_ in sheets.items():
            df_.to_excel(w, sheet_name=nm, index=False)
    print(f"\nSource data: {source_base}.xlsx  (sheets: {', '.join(sheets)})")
except Exception as e:
    print(f"\nWarning: XLSX export failed: {e}. Writing TSV files.")

for nm, df_ in sheets.items():
    df_.to_csv(f"{source_base}_{nm}.tsv", sep="\t", index=False)
print(f"Source data (tsv): {source_base}_*.tsv")
print("=" * 70)
