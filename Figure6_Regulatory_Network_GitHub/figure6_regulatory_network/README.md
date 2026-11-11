# Figure 6 — Multi-layered regulatory evidence

Plotting scripts for Figure 6 of *Multi-omic profiling and MPRA of hybrid pigs uncover extensive functional cis-regulatory variation undetectable in purebreds*. The scripts use processed DLC and PigBiobank results. Use the manuscript's Data Availability statement to locate study inputs.

## Setup

Use Python 3.10 or newer. Run the commands from this directory and replace the example paths with the files used for the manuscript.

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
```

## Fig. 6a — emVar–gene–trait network

| Input | Expected contents |
| --- | --- |
| Pathway TSV or TSV.gz | One row per `emVar`, `gene` and `trait` combination. Include `ASE`, `cis_eQTL`, `gene_based`, `TWAS_Muscle`, `SMR_Muscle`, `enloc_Muscle`, `SNP_GWAS_P_lt_5e8` and `gene_based_min_p`. Use the original position-and-sorted-allele emVar keys, e.g. `7:23797335:C/T`. |
| Gene metadata TSV or TSV.gz | Original v19 metadata table with `gene_id_key` and optional `Gene.Name`. |
| Trait category XLSX | `Main_Category` (or `Main Category`) and a trait key in `Trait`, `Study Symbol Raw` or `study_symbol`. |
| CRE and 3′UTR MPRA TSVs | Both tables require `chr`, `variant_pos`, `ref_allele`, `alt_allele`, `PSC.emVar_act_qc` and `PK15.emVar_act_qc`. |

```bash
python scripts/plot_fig6a_emvar_gene_trait_network.py \
  --pathways data/05_emVar_supported_pathways_with_enloc.tsv.gz \
  --v19-metadata data/01_all_nonredundant_MPRA_SNP_gene_trait.tsv.gz \
  --categories data/PigBiobank_trait264.list_3_cols.xlsx \
  --cre-mpra data/RESULT_CROSSCELL_M3_cre_A0_F_eqtlase_filtered.tsv \
  --utr-mpra data/RESULT_CROSSCELL_M3_utr_A0_F_eqtlase_filtered.tsv \
  --output results/fig6a --run
```

Omit `--run` for an input precheck. The plotted subset requires gene-based GWAS support and excludes Exterior and Adaptation. The upper panel labels the smallest raw gene-based GWAS P value per category, prioritizing paths with more external evidence types in case of ties. The command exports PDF and PNG figures and the single-sheet `figure_source_data_v34.xlsx`. Use a fresh output directory or explicitly pass `--force` to replace same-named outputs.

## Fig. 6b — Regulatory evidence at chr7:23,797,335 C>T

Use Sscrofa11.1 coordinates. The plot spans ±500 kb around the focal variant. The script joins plotted GWAS, cis-eQTL, LD and MPRA records by position within chromosome 7 as in the supplied plotting code.

| Input | Expected contents |
| --- | --- |
| PigBiobank M_LMDEP GWAS | Gzip-compressed whitespace-delimited summary table with chromosome, position and P value (e.g. `chromosome`, `position`, `p_value`). Optional REF, ALT and sample-size columns are retained when available. |
| DLC cis-eQTL TSV | `SNP_id`, `Gene.ID`, `p.value`, `BBFDR` and `beta` or supported aliases. |
| DLC LD table | PLINK whitespace-delimited LD with `CHR_B`, `BP_B` and `R2`, computed relative to chr7:23,599,787. |
| Gene models | Ensembl release 109 GTF for Sscrofa11.1. |
| CRE MPRA TSV | `chr`, `variant_pos`, `PK15.q_BH_all`, `PSC.q_BH_all`, `PK15.emVar_all_qc` and `PSC.emVar_all_qc`. |
| Chromatin tracks | ATAC, CTCF, H3K27ac, H3K4me1 and H3K4me3 muscle bigWigs. |

```bash
python scripts/plot_fig6b_chr7_regulatory_locus.py \
  --gwas data/M_LMDEP.txt.gz --eqtl data/eQTL_MAF005_kept.tsv \
  --ld data/DLC_LD_chr7_23599787_lead_eQTL.ld \
  --gtf data/Sus_scrofa.Sscrofa11.1.109.gtf \
  --mpra data/RESULT_CROSSCELL_M3_cre_A0_F_eqtlase_filtered.tsv \
  --atac data/ATAC_Muscle_P350_RPKM.bigwig --ctcf data/CTCF_Muscle_P350_RPKM.bigwig \
  --h3k27ac data/H3K27ac_Muscle_P350_RPKM.bigwig \
  --h3k4me1 data/H3K4me1_Muscle_P350_RPKM.bigwig \
  --h3k4me3 data/H3K4me3_Muscle_P350_RPKM.bigwig \
  --outdir results/fig6b
```

Outputs include PDF and PNG figures, a filtering summary and panel-level Source Data in Excel and TSV format. The locus view uses the CRE MPRA table and the displayed cis-eQTL positions.
