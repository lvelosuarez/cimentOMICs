#!/usr/bin/env python3
"""
compute_arg_rpkm.py
-------------------
Compute proper RPKM for all ARG-carrying contigs using pre-summarised
read count tables (generated from raw idxstats by summarize_idxstats.py).

RPKM = (mapped_reads_to_contig / contig_length_kb) / library_size_M

Inputs
------
  data/sample_metadata.tsv  — total_mapped_reads column
  data/arg_contigs.tsv      — contig list with length, mapped_reads, plasmid status

Outputs
-------
  data/arg_rpkm_mapped.tsv      — per-(sample, contig, ARO) RPKM values
  data/macrolide_rpkm_mapped.tsv — macrolide subset, matches fig2 input

To regenerate the input tables from raw idxstats:
  python scripts/summarize_idxstats.py
"""

import re
from pathlib import Path
import pandas as pd

# --- Paths -------------------------------------------------------------------
BASE      = Path(__file__).parent.parent
ARO_TABLE = BASE / "data" / "merged_aro_per_sample.tsv"
OUT_ALL   = BASE / "data" / "arg_rpkm_mapped.tsv"
OUT_MAC   = BASE / "data" / "macrolide_rpkm_mapped.tsv"

# --- Load pre-summarised mapping tables --------------------------------------
meta        = pd.read_csv(BASE / "data" / "sample_metadata.tsv", sep="\t")
lib_sizes   = meta.set_index("sample_id")["total_mapped_reads"]
contigs_tbl = pd.read_csv(BASE / "data" / "arg_contigs.tsv", sep="\t")
counts      = contigs_tbl[["sample_id", "contig", "length", "mapped_reads"]].drop_duplicates()
arg_contigs = contigs_tbl[["sample_id", "contig"]].drop_duplicates()

# --- Load ARO annotations (to get drug class / macrolide flag) ---------------
aro = pd.read_csv(ARO_TABLE, sep="\t")[
    ["sample_id", "ARO_id", "gene_label", "drug_class", "is_mac", "detected_by"]
].drop_duplicates()

# --- Compute RPKM ------------------------------------------------------------
records = []

for sample, df_arg in counts.groupby("sample_id"):
    if sample not in lib_sizes.index:
        print(f"[WARN] {sample}: no library size — skipping")
        continue

    library_size_M = lib_sizes[sample] / 1e6
    if library_size_M == 0:
        print(f"[WARN] {sample}: zero mapped reads — skipping")
        continue

    df_arg = df_arg.copy()
    df_arg["length_kb"] = df_arg["length"] / 1000.0
    df_arg["rpkm"] = df_arg["mapped_reads"] / (df_arg["length_kb"] * library_size_M)
    df_arg["library_size_M"] = library_size_M

    records.append(df_arg[["sample_id", "contig", "length", "mapped_reads",
                            "library_size_M", "rpkm"]])

    print(f"[OK] {sample}: {len(df_arg)} ARG contigs, "
          f"library={library_size_M:.2f}M reads, "
          f"median RPKM={df_arg['rpkm'].median():.1f}")

if not records:
    raise RuntimeError(
        "No data found. Regenerate input tables with:\n"
        "  python scripts/summarize_idxstats.py"
    )

result = pd.concat(records, ignore_index=True)

# --- Merge with ARO annotations ----------------------------------------------
# arg_mge_colocalisation has one row per contig (not per ARO);
# merged_aro_per_sample has one row per (sample, ARO).
# We need contig → ARO mapping. Build it from the funcscan output path or
# approximate: join on sample_id + contig via arg_mge_colocalisation,
# then bring in ARO via sample_id only (some contigs carry multiple AROs).

# Load full contig→ARO mapping from the RGI/AMRFinder combined table
# (arg_mge_colocalisation has contig but not ARO; we need the contig_id from RGI)
# Fallback: merge result with aro on sample_id only for drug class annotation,
# then aggregate RPKM per sample per drug_class.

result_full = result.merge(
    arg_contigs,   # sample_id, contig
    on=["sample_id", "contig"],
    how="left"
)

# --- Save per-contig RPKM ----------------------------------------------------
result_full.to_csv(OUT_ALL, sep="\t", index=False)
print(f"\n[SAVED] {OUT_ALL}")

# --- Macrolide subset --------------------------------------------------------
# Use the is_mac flag from the ARO table joined on sample_id.
# Since we don't have per-contig ARO here, use the contig list from the
# arg_mge_colocalisation table filtered to macrolide AROs via merged_aro table.

# Get macrolide ARG contig list
mac_aro = aro[aro["is_mac"] == True][["sample_id", "ARO_id", "gene_label"]].drop_duplicates()

# arg_mge_colocalisation doesn't have ARO_id — need to join via the original
# RGI/AMRFinder output. For now, flag contigs that belong to macrolide-detected
# samples and match drug class in a separate pass.
# Simpler: load the contig-level drug class from the colocalisation table
# (it has broad_group but not drug_class per ARO).
# Use the known macrolide gene labels to identify macrolide contigs.

# Load a richer contig table if available
arg_full = pd.read_csv(ARG_CONTIGS, sep="\t")

# Check if drug_class column exists (it may not — depends on pipeline step)
if "drug_class" in arg_full.columns:
    mac_contigs = arg_full[arg_full["drug_class"].str.contains(
        "macrolide|MLSB", case=False, na=False
    )][["sample_id", "contig"]].drop_duplicates()
else:
    # Fall back to: any contig in a sample that has macrolide AROs in aro table
    # This is conservative — we include the contig if the sample has any mac ARO
    # and the contig is an ARG contig. Better to keep all and note limitation.
    print("[NOTE] No drug_class column in arg_mge_colocalisation — "
          "macrolide subset uses sample-level flag from merged_aro_per_sample")
    mac_samples = mac_aro["sample_id"].unique()
    mac_contigs = arg_contigs[arg_contigs["sample_id"].isin(mac_samples)]

mac_result = result_full.merge(mac_contigs, on=["sample_id", "contig"], how="inner")

# Per-sample total macrolide RPKM (sum across all macrolide ARG contigs)
mac_per_sample = (
    mac_result
    .groupby("sample_id")["rpkm"]
    .sum()
    .reset_index()
    .rename(columns={"rpkm": "total_mac_rpkm_mapped"})
)

mac_per_sample.to_csv(OUT_MAC, sep="\t", index=False)
print(f"[SAVED] {OUT_MAC}")

# --- Summary stats -----------------------------------------------------------
print("\n--- RPKM Summary (all ARG contigs) ---")
print(result_full.groupby("sample_id")["rpkm"].agg(["count", "median", "sum"])
      .rename(columns={"count": "n_arg_contigs", "median": "median_rpkm",
                       "sum": "total_rpkm"})
      .round(2)
      .to_string())

print("\n--- Macrolide RPKM per sample ---")
print(mac_per_sample.sort_values("total_mac_rpkm_mapped", ascending=False).to_string(index=False))
