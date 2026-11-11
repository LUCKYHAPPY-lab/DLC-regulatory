#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""V5: population-specific site-level Hp, retaining fixed sites as zero.

Same counting rule as V1 except that monomorphic sites are retained.
No cross-population site intersection, MAF filter or new AN>=100 filter.
Missing/invalid counts are excluded, never imputed as homozygous reference.
"""
import argparse
import csv
import gzip
import hashlib
import json
import os
import shlex
import shutil
import signal
import subprocess
import sys
from array import array
from datetime import datetime
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

POPS = ('DU', 'LC', 'DLC')
EXPECTED = {'DU': 100, 'LC': 100, 'DLC': 189}
COLORS = {'DU': '#D94841', 'LC': '#4A98BD', 'DLC': '#F5B85C'}


def table(path, rows):
    with open(path, 'w', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter='\t')
        w.writeheader()
        w.writerows(rows)


def validate(vcf, sample_dir, binary, out):
    names = subprocess.check_output([binary, 'query', '-l', str(vcf)], text=True).splitlines()
    if len(names) != len(set(names)):
        raise ValueError('Duplicate VCF sample IDs')
    assigned = {}
    files = {}
    manifest = []
    dst = out / 'sample_lists'
    dst.mkdir()
    for pop in POPS:
        ids = (sample_dir / (pop + '.txt')).read_text().split()
        if len(ids) != EXPECTED[pop] or len(set(ids)) != EXPECTED[pop]:
            raise ValueError(f'{pop}: expected {EXPECTED[pop]} unique sample IDs')
        for name in ids:
            if name in assigned:
                raise ValueError(f'Sample occurs in multiple groups: {name}')
            assigned[name] = pop
            manifest.append({'Sample_ID': name, 'Population': pop})
        files[pop] = dst / (pop + '.txt')
        files[pop].write_text('\n'.join(ids) + '\n')
    if set(names) != set(assigned):
        raise ValueError(f'Lists do not partition VCF: {len(set(assigned)-set(names))} list-only IDs; {len(set(names)-set(assigned))} VCF-only IDs')
    table(out / 'sample_manifest.tsv', manifest)
    return files


def group(pop, vcf, samples, binary, out):
    fmt = '%CHROM\\t%POS\\t%REF\\t%ALT\\t%INFO/AN\\t%INFO/AC\\n'
    q = shlex.quote
    cmd = ('set -o pipefail; '
           f'{q(binary)} view -S {q(str(samples))} -Ou {q(str(vcf))} | '
           f'{q(binary)} +fill-tags -Ou -- -t AC,AN | '
           f'{q(binary)} query -f {q(fmt)}')
    source = out / f'{pop}_Hp_by_site_including_fixed.tsv.gz'
    log_path = out / f'{pop}_bcftools.log'
    print(f'[{pop}] Starting; see {log_path}', flush=True)
    print(cmd, flush=True)
    values = array('d')
    counts = {'Population': pop, 'N_samples': EXPECTED[pop],
              'N_VCF_records_scanned': 0, 'N_valid_sites_retained': 0,
              'N_polymorphic_sites_retained': 0, 'N_fixed_sites_retained_Hp_zero': 0,
              'N_missing_or_invalid_sites_excluded': 0}
    with open(log_path, 'w') as log, gzip.open(source, 'wt', compresslevel=3) as dst:
        log.write(cmd + '\n')
        log.flush()
        proc = subprocess.Popen(cmd, shell=True, executable='/bin/bash', stdout=subprocess.PIPE,
                                stderr=log, text=True, start_new_session=True)
        try:
            dst.write('Population\tCHROM\tPOS\tREF\tALT\tAN\tAC\tN_observed_alleles\tSite_status\tHp\n')
            for line in proc.stdout:
                counts['N_VCF_records_scanned'] += 1
                fields = line.rstrip('\n').split('\t')
                if len(fields) != 6:
                    raise ValueError(f'{pop}: malformed bcftools query row')
                chrom, pos, ref, alt, an_txt, ac_txt = fields
                valid = True
                try:
                    an = int(an_txt)
                    ac = [int(x) for x in ac_txt.split(',')]
                    ref_count = an - sum(ac)
                    valid = an > 0 and all(x >= 0 for x in ac) and ref_count >= 0
                except ValueError:
                    valid = False
                if not valid:
                    counts['N_missing_or_invalid_sites_excluded'] += 1
                else:
                    allele_counts = [ref_count] + ac
                    observed = sum(x > 0 for x in allele_counts)
                    fixed = observed == 1
                    hp = 0.0 if fixed else max(0.0, min(1.0, 1.0-sum((x/an)**2 for x in allele_counts)))
                    values.append(hp)
                    counts['N_valid_sites_retained'] += 1
                    counts['N_fixed_sites_retained_Hp_zero' if fixed else 'N_polymorphic_sites_retained'] += 1
                    dst.write('\t'.join([pop, chrom, pos, ref, alt, str(an), ac_txt,
                                          str(observed), 'Fixed' if fixed else 'Polymorphic', format(hp, '.17g')])+'\n')
                if counts['N_VCF_records_scanned'] % 1000000 == 0:
                    print(f"[{pop}] scanned {counts['N_VCF_records_scanned']:,}; retained {counts['N_valid_sites_retained']:,}; fixed {counts['N_fixed_sites_retained_Hp_zero']:,}", flush=True)
            rc = proc.wait()
            if rc:
                raise subprocess.CalledProcessError(rc, cmd)
        finally:
            proc.stdout.close()
            if proc.poll() is None:
                os.killpg(proc.pid, signal.SIGTERM)
                try:
                    proc.wait(timeout=10)
                except subprocess.TimeoutExpired:
                    os.killpg(proc.pid, signal.SIGKILL)
                    proc.wait()
    if not values:
        raise RuntimeError(f'{pop}: no valid sites')
    if counts['N_VCF_records_scanned'] != counts['N_valid_sites_retained'] + counts['N_missing_or_invalid_sites_excluded']:
        raise AssertionError('Record accounting mismatch')
    x = np.frombuffer(values, dtype=np.float64)
    n = len(x)
    mean = float(np.mean(x))
    sd = float(np.std(x, ddof=1)) if n > 1 else float('nan')
    # In-place sorting keeps memory bounded to one population; quartiles exact.
    x.sort()
    q1, median, q3 = (float(z) for z in np.quantile(x, [.25, .5, .75]))
    iqr = q3 - q1
    lo = int(np.searchsorted(x, q1 - 1.5*iqr, side='left'))
    hi = int(np.searchsorted(x, q3 + 1.5*iqr, side='right')) - 1
    whislo, whishi = float(x[min(lo, n-1)]), float(x[max(hi, 0)])
    summary = {**counts, 'Fixed_site_fraction': counts['N_fixed_sites_retained_Hp_zero']/n,
               'Mean_Hp': mean, 'SD_Hp': sd, 'Median_Hp': median, 'Q1_Hp': q1, 'Q3_Hp': q3,
               'Min_Hp': float(x[0]), 'Max_Hp': float(x[-1]), 'Whisker_low': whislo, 'Whisker_high': whishi}
    plot = dict(label=pop, q1=q1, med=median, q3=q3, whislo=whislo, whishi=whishi, fliers=[])
    table(out / f'{pop}_summary.tsv', [summary])
    print(f"[{pop}] Finished: retained={n:,}; fixed={counts['N_fixed_sites_retained_Hp_zero']:,}; median={median:.9g}", flush=True)
    return summary, plot


def main():
    if hasattr(sys.stdout, 'reconfigure'):
        sys.stdout.reconfigure(line_buffering=True)
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--vcf', type=Path, required=True)
    parser.add_argument('--sample-dir', type=Path, required=True)
    parser.add_argument('--out-root', type=Path, required=True)
    parser.add_argument('--bcftools', default='bcftools')
    a = parser.parse_args()
    if not a.vcf.is_file():
        raise FileNotFoundError(a.vcf)
    if not shutil.which(a.bcftools):
        raise FileNotFoundError(a.bcftools)
    out = a.out_root / datetime.now().strftime('Hp_V5_independent_fixed_zero_%Y%m%d_%H%M%S_%f')
    out.mkdir(parents=True, exist_ok=False)
    shutil.copyfile(__file__, out / Path(__file__).name)
    metadata = {'parameters': {k: str(v) for k,v in vars(a).items()}, 'version': 'V5',
                'script_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                'bcftools_version': subprocess.check_output([a.bcftools, '--version'], text=True),
                'rule': 'Independent population subsets; retain every record with valid AC/AN and AN>0. No MAF, no cross-population callability, no added minimum-called-individual filter. Fixed sites retained as Hp=0.',
                'formula': 'Hp = 1 - sum(p_i^2); for biallelic records Hp = 2pq.',
                'scope': 'All valid variant records in this input VCF, including SNPs, indels and multiallelic records. Same variant-type scope as V1. Does not infer genotypes at records absent from the original population callset.',
                'limits': 'GT-derived allele counts, not per-individual observed heterozygosity or per-base genomic diversity. No added DP/GQ filter. Different site ascertainment across populations limits direct ranking. No SNP-wise significance test because linked sites are not independent.'}
    (out / 'run_metadata.json').write_text(json.dumps(metadata, indent=2)+'\n')
    print(f'Output: {out}', flush=True)
    print('V5: independent valid site sets; fixed retained; missing excluded; no intersection; no MAF filter.', flush=True)
    sample_files = validate(a.vcf, a.sample_dir, a.bcftools, out)
    summaries, plot_stats = [], []
    for pop in POPS:
        summary, plot = group(pop, a.vcf, sample_files[pop], a.bcftools, out)
        summaries.append(summary)
        plot_stats.append(plot)
        table(out / 'Fig2e_Hp_V5_summary.tsv', summaries)
    plt.rcParams['pdf.fonttype'] = 42
    plt.rcParams['ps.fonttype'] = 42
    fig, ax = plt.subplots(figsize=(5.5, 5))
    bp = ax.bxp(plot_stats, showfliers=False, patch_artist=True, widths=.6,
                medianprops=dict(color='black', linewidth=1.5))
    for box, pop in zip(bp['boxes'], POPS):
        box.set_facecolor(COLORS[pop])
        box.set_alpha(.9)
    ax.set_ylabel(r'Expected heterozygosity ($H_p$)', fontsize=13)
    ax.set_xticklabels([f"{r['Population']}\n(n={r['N_valid_sites_retained']:,} sites)" for r in summaries], fontsize=9)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    fig.tight_layout()
    for suffix in ('pdf', 'png'):
        fig.savefig(out / ('Fig2e_Hp_V5_independent_fixed_zero.'+suffix), dpi=600, bbox_inches='tight')
    plt.close(fig)
    (out / 'README.txt').write_text(
        'V5: independent population-specific valid-site Hp distributions, including fixed sites.\n'
        'This is V1 with monomorphic records retained, without adding a shared-site or AN>=100 filter.\n'
        'AN>0 plus valid nonnegative AC counts is required independently within each population.\n'
        'No MAF filter. Missing/invalid allele counts excluded, never changed to Hp=0.\n'
        'SNPs, indels and multiallelic records are included, matching the V1 record-type scope.\n'
        'Multiallelic Hp can exceed 0.5. AN is the number of called alleles, not individual sample size.\n'
        'Each source table row is one retained VCF record, with allele counts, status and Hp.\n'
        'Summary includes exact quartiles, boxplot whiskers, fixed counts and missing/invalid counts.\n'
        'Site sets differ across groups; their distributions are descriptive and subject to ascertainment.\n'
        'These are variant-record estimates, not individual Ho or diversity across all genomic bases.\n'
        'No site-level significance test is performed.\n', encoding='utf-8')
    (out / 'COMPLETED.txt').write_text(datetime.now().isoformat()+'\n')
    (a.out_root / 'latest_completed.txt').write_text(str(out.resolve())+'\n')
    print(f'COMPLETED: {out}', flush=True)


if __name__ == '__main__':
    main()
