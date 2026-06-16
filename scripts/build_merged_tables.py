#!/usr/bin/env python3
"""
build_merged_tables.py
----------------------
Consolidates processed data tables into two clean master files:

  data/sample_metadata.tsv   — one row per sample (41 rows)
  data/arg_contigs.tsv       — one row per ARG-carrying contig (791 rows)

Run once after any upstream update to regenerate these files.
"""

from pathlib import Path
import pandas as pd

BASE = Path(__file__).parent.parent
D    = BASE / "data"

# ── 1. sample_metadata.tsv ────────────────────────────────────────────────────
summary  = pd.read_csv(D / "merged_summary.tsv",            sep="\t")
sylph    = pd.read_csv(D / "sylph_ciment_summary.tsv",      sep="\t")
libsize  = pd.read_csv(D / "libsize_pa_continuous.tsv",     sep="\t")
cftr     = pd.read_csv(D / "cftr_status.tsv",               sep="\t")
mac_rpkm = pd.read_csv(D / "macrolide_rpkm_mapped.tsv",     sep="\t")
mobile   = pd.read_csv(D / "genomad_mobile_per_sample.tsv", sep="\t")
libmap   = pd.read_csv(D / "mapping_library_sizes.tsv",     sep="\t")

# Start from summary (has Leeds, AZM, PA_sylph_pct, ARG richness counts)
meta = summary.copy()

# Sylph: add species fractions (PA_sylph_pct and AZM already in summary)
sylph_cols = ["sample_id", "Ach_sylph_pct", "Steno_sylph_pct", "Haem_sylph_pct",
              "Rothia_sylph_pct", "Strep_sylph_pct", "PA_pos_sylph"]
meta = meta.merge(sylph[sylph_cols], on="sample_id", how="left")

# Library size (reads_M from libsize_pa_continuous; total_mapped_reads from idxstats)
meta = meta.merge(libsize[["sample_id", "reads_M"]], on="sample_id", how="left")
meta = meta.merge(libmap[["sample_id", "total_mapped_reads"]], on="sample_id", how="left")

# CFTR status
meta = meta.merge(cftr[["sample_id", "illumina_id", "ivacaftor", "orkambi", "cftr_pos"]],
                  on="sample_id", how="left")

# Macrolide RPKM (AZM and Leeds already present)
meta = meta.merge(mac_rpkm[["sample_id", "total_mac_rpkm_mapped"]], on="sample_id", how="left")

# Mobile ARG counts (AZM and PA_pct already present)
mobile_cols = ["sample_id", "n_total_arg", "n_mobile_arg", "n_chrom_arg", "pct_mobile"]
meta = meta.merge(mobile[mobile_cols], on="sample_id", how="left")

meta = meta.sort_values("sample_id").reset_index(drop=True)
meta.to_csv(D / "sample_metadata.tsv", sep="\t", index=False)
print(f"[OK] sample_metadata.tsv — {len(meta)} rows × {len(meta.columns)} columns")
print(f"     columns: {', '.join(meta.columns)}\n")

# ── 2. arg_contigs.tsv ────────────────────────────────────────────────────────
mge  = pd.read_csv(D / "arg_mge_colocalisation.tsv", sep="\t")
rpkm = pd.read_csv(D / "arg_rpkm_mapped.tsv",        sep="\t")

# rpkm has: sample_id, contig, length, mapped_reads, library_size_M, rpkm
# mge  has: sample_id, contig, status, taxid, organism, is_PA, broad_group,
#           PA_pct, on_plasmid, plasmid_score, n_hallmarks, conjugation_genes
contigs = mge.merge(
    rpkm[["sample_id", "contig", "length", "mapped_reads", "library_size_M", "rpkm"]],
    on=["sample_id", "contig"],
    how="left"
)
contigs = contigs.sort_values(["sample_id", "contig"]).reset_index(drop=True)
contigs.to_csv(D / "arg_contigs.tsv", sep="\t", index=False)
print(f"[OK] arg_contigs.tsv — {len(contigs)} rows × {len(contigs.columns)} columns")
print(f"     columns: {', '.join(contigs.columns)}")
