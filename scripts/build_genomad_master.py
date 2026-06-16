"""
build_genomad_master.py
-----------------------
Collects geNomad plasmid_summary.tsv files from all 40 samples on tachyon
and builds a single master CSV analogous to plasx_master.csv.

Filters applied (informed by TNEG QC):
  - contig length >= 1500 bp   (removes TNEG false-positive noise at 206-351 bp)
  - contig coverage > 2x        (removes 1x assembly artifact calls)
  - plasmid_score >= 0.5        (permissive; sensitivity analysis downstream)

Output: data/genomad_master.csv
Columns: ID_meta, node_num, contig_length, contig_coverage,
         plasmid_score, n_hallmarks, marker_enrichment, conjugation_genes, amr_genes

Usage:
  # On local machine after rsync (see below)
  python3 scripts/build_genomad_master.py

Rsync command to fetch results from tachyon:
  rsync -av --progress \\
    lourdes@tachyon:/mnt/san/microbio/projects/ciment/genomad/CO*_summary/ \\
    data/genomad/
"""

import pandas as pd
import re
from pathlib import Path
import warnings
warnings.filterwarnings("ignore")

BASE        = Path(__file__).parent.parent
GENOMAD_DIR = BASE / "data" / "genomad"
OUT_FILE    = BASE / "data" / "genomad_master.csv"

# QC filters (from TNEG analysis)
MIN_LENGTH   = 1500   # bp — TNEG artifacts were 241/351 bp
MIN_COVERAGE = 2.0    # x  — TNEG artifacts at ~1x
MIN_SCORE    = 0.5    # plasmid_score threshold (same as PlasX primary threshold)

def parse_contig_header(seq_name):
    """
    Extract node_num, length, and coverage from metaSPAdes header.
    Format: NODE_<num>_length_<len>_cov_<cov>
    geNomad may append |provirus_start_end for proviruses — strip that.
    """
    seq_name = seq_name.split("|")[0]  # strip provirus suffix
    m = re.match(r"NODE_(\d+)_length_(\d+)_cov_([\d.]+)", seq_name)
    if m:
        return int(m.group(1)), int(m.group(2)), float(m.group(3))
    return None, None, None

def load_sample(sample_dir: Path) -> pd.DataFrame:
    """Load plasmid_summary.tsv for one sample."""
    sample = sample_dir.name  # e.g. "CO1_summary" → ID_meta = "CO1"
    id_meta = sample.replace("_summary", "")

    tsv = sample_dir / f"{id_meta}_plasmid_summary.tsv"
    if not tsv.exists():
        print(f"  MISSING: {tsv}")
        return pd.DataFrame()

    df = pd.read_csv(tsv, sep="\t", na_values=["NA"])

    # Parse contig metadata from seq_name
    parsed = df["seq_name"].apply(parse_contig_header)
    df["node_num"]       = [p[0] for p in parsed]
    df["contig_length"]  = [p[1] for p in parsed]
    df["contig_coverage"]= [p[2] for p in parsed]
    df["ID_meta"]        = id_meta

    return df

def main():
    if not GENOMAD_DIR.exists():
        print(f"ERROR: {GENOMAD_DIR} not found.")
        print("Run rsync first:")
        print("  rsync -av lourdes@tachyon:/mnt/san/microbio/projects/ciment/genomad/CO*_summary/ data/genomad/")
        return

    sample_dirs = sorted(GENOMAD_DIR.glob("CO*_summary"))
    print(f"Found {len(sample_dirs)} sample directories")

    frames = []
    for d in sample_dirs:
        df = load_sample(d)
        if not df.empty:
            frames.append(df)

    if not frames:
        print("No data loaded.")
        return

    master = pd.concat(frames, ignore_index=True)
    print(f"\nTotal contigs before filtering: {len(master):,}")

    # Apply QC filters
    before = len(master)
    master = master[
        (master["contig_length"]   >= MIN_LENGTH)   &
        (master["contig_coverage"] >= MIN_COVERAGE) &
        (master["plasmid_score"]   >= MIN_SCORE)
    ]
    print(f"After QC filters (len≥{MIN_LENGTH}, cov≥{MIN_COVERAGE}, score≥{MIN_SCORE}): {len(master):,} ({before-len(master):,} removed)")

    # Select and rename output columns
    keep = ["ID_meta", "node_num", "contig_length", "contig_coverage",
            "plasmid_score", "n_hallmarks", "marker_enrichment",
            "conjugation_genes", "amr_genes"]
    master = master[[c for c in keep if c in master.columns]]

    # Check expected samples
    found = sorted(master["ID_meta"].unique())
    expected = [f"CO{i}" for i in range(1, 41)]
    missing = [s for s in expected if s not in found]
    if missing:
        print(f"\nWARNING: Missing samples: {missing}")
    print(f"Samples with data: {len(found)}/40")

    # Summary stats
    print(f"\nPlasmid-classified contigs per sample:")
    print(master.groupby("ID_meta")["node_num"].count().describe())
    print(f"\nScore distribution:")
    print(master["plasmid_score"].describe())

    # Check for AMR-carrying plasmids
    amr = master[master["amr_genes"].notna()]
    print(f"\nPlasmids with AMR genes annotated by geNomad: {len(amr)}")
    if len(amr) > 0:
        print(amr[["ID_meta", "node_num", "contig_length", "plasmid_score", "amr_genes"]].to_string(index=False))

    master.to_csv(OUT_FILE, index=False)
    print(f"\nSaved → {OUT_FILE}")
    print(f"Shape: {master.shape}")

if __name__ == "__main__":
    main()
