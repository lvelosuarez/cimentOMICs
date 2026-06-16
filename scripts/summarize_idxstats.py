#!/usr/bin/env python3
"""
summarize_idxstats.py
---------------------
Collapses per-sample samtools idxstats files into two compact tables:

  data/mapping_library_sizes.tsv      — total mapped reads per sample
  data/mapping_arg_contig_counts.tsv  — read counts for ARG contigs only

Run once after fetching raw idxstats from the sequencing server.
The two output tables are sufficient to recompute all RPKM values via
compute_arg_rpkm.py and are much smaller than the full idxstats (~28 MB → ~200 KB).
"""

from pathlib import Path
import pandas as pd

BASE         = Path(__file__).parent.parent
IDXSTATS_DIR = BASE / "data" / "mapping_idxstats"
OUT_LIBSIZES = BASE / "data" / "mapping_library_sizes.tsv"
OUT_COUNTS   = BASE / "data" / "mapping_arg_contig_counts.tsv"

arg_contigs = pd.read_csv(BASE / "data" / "arg_contigs.tsv", sep="\t")[["sample_id", "contig"]].drop_duplicates()

lib_rows, count_rows = [], []

for f in sorted(IDXSTATS_DIR.glob("*.idxstats.tsv")):
    sample = f.stem.replace(".idxstats", "")

    df = pd.read_csv(f, sep="\t", header=None,
                     names=["contig", "length", "mapped_reads", "unmapped_reads"])

    total_mapped = df.loc[df["contig"] != "*", "mapped_reads"].sum()
    lib_rows.append({"sample_id": sample, "total_mapped_reads": total_mapped})

    sample_arg = arg_contigs.loc[arg_contigs["sample_id"] == sample, "contig"].tolist()
    df_arg = df[df["contig"].isin(sample_arg)][["contig", "length", "mapped_reads"]].copy()
    df_arg.insert(0, "sample_id", sample)
    count_rows.append(df_arg)

    print(f"[OK] {sample}: library={total_mapped/1e6:.2f}M reads, {len(df_arg)} ARG contigs")

pd.DataFrame(lib_rows).sort_values("sample_id").to_csv(OUT_LIBSIZES, sep="\t", index=False)
print(f"\n[SAVED] {OUT_LIBSIZES}")

pd.concat(count_rows, ignore_index=True).to_csv(OUT_COUNTS, sep="\t", index=False)
print(f"[SAVED] {OUT_COUNTS}")
