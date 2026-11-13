# Allele-specific expression analysis

Upstream ASE analysis scripts for the Duroc × Luchuan F1 study. Use Python with pandas and R with the poolr package. Run the commands from this directory.

## 1. Per-sample tests

Start from GATK ASEReadCounter TSV files named `<sample>_ASEReadCounter.table`. Sample IDs are the filename prefix before the first underscore. The script applies total depth ≥10, at least three reads per allele and at least 1% per allele. It performs a two-sided binomial test with null proportion 0.5 and applies BH adjustment within each sample.

```bash
Rscript scripts/01_sample_binomial_tests.R data/read_counts results/sample_tests
```

## 2. Merge sample results

Use the exonic DNA variant VCF (plain or gzip) and the per-sample tables with standard GATK columns `contig`, `position`, `refAllele`, `altAllele`, `refCount`, `altCount` and the generated `FDR`. The join uses the exact chromosome, position, REF and ALT strings in that order. The output contains `<sample>_info` (reference count, alternative count and FDR separated by underscores) and `<sample>_FDR`.

```bash
python -m pip install -r requirements.txt
python scripts/02_merge_sample_results.py --variants data/exonic_variants.vcf.gz \
  --sample-dir results/sample_tests --output results/ASE_result_189_sample.txt
```

## 3. Meta-analysis and annotation

Use the DNA-genotype-masked version of the merged table (`ASE_result_189_sample.True.txt` in the original workflow). The masking step must set sample–variant observations without a qualifying DNA heterozygous genotype to missing. The supplied scripts start from the already masked table for this step. Provide the comma-separated sample map with `Ear_ID` and the original tab-delimited annotation export with `Uploaded_variation` and the original column layout.

```bash
Rscript scripts/03_meta_analysis.R results/ASE_result_189_sample.True.txt \
  data/sample_map.csv data/ASE_result_189_sample_SNP_output.vcf results/meta
```

The script matches sample columns by name, combines per-sample binomial P values with the original poolr binomial procedure for sites with at least 19 qualifying samples, and summarizes allelic ratios after rounding to three decimals. It exports unannotated and annotated TSVs. The annotated output applies BH adjustment after the annotation join, following the supplied script. Annotation multiplicity therefore affects the correction universe.

The annotated TSV provides an upstream input for Fig. 4b when the annotation export supplies `Gene`. The Fig. 4 module applies the downstream annotation, MAF and Mratio filters. GATK counting, DNA genotype masking and generation of the annotation table are external steps.
