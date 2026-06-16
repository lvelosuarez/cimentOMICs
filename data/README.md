# Data directory

Processed summary tables from the CIMENTOMICS metagenomics pipeline (n=40 CF sputum samples).

---

## Files in this repository

### `sample_metadata.tsv` — 40 rows × 29 columns

One row per patient. Consolidates all per-sample metadata.

| Column | Description |
|--------|-------------|
| `sample_id` | Patient identifier (CO1–CO40) |
| `Leeds` | Leeds microbiology subgroup (LeedsA / LeedsB) |
| `AZM` | Azithromycin prophylaxis status (AZMpos / AZMneg) |
| `PA_sylph_pct` | *P. aeruginosa* relative abundance (%, sylph) |
| `PA_pos_sylph` | Boolean: PA detected (sylph) |
| `Ach_sylph_pct` | *Achromobacter* spp. relative abundance (%) |
| `Steno_sylph_pct` | *Stenotrophomonas* spp. relative abundance (%) |
| `Haem_sylph_pct` | *Haemophilus* spp. relative abundance (%) |
| `Rothia_sylph_pct` | *Rothia* / *Actinomyces* relative abundance (%) |
| `Strep_sylph_pct` | *Streptococcus* / oral Firmicutes relative abundance (%) |
| `n_rgi_strict` | Unique ARGs detected by RGI (STRICT) |
| `n_amrf` | Unique ARGs detected by AMRFinderPlus |
| `n_consensus` | ARGs detected by both tools |
| `n_rgi_only` / `n_amrf_only` | ARGs detected by one tool only |
| `n_mac_rgi` / `n_mac_amrf` / `n_mac_consensus` | Macrolide ARG counts per tool |
| `reads_M` | Total library size (millions of reads) |
| `total_mapped_reads` | Reads mapped to assembly (Bowtie2) |
| `illumina_id` | Internal sequencing ID |
| `ivacaftor` / `orkambi` | CFTR modulator use (1 = yes) |
| `cftr_pos` | Boolean: on any CFTR modulator |
| `total_mac_rpkm_mapped` | Total macrolide ARG abundance (sum RPKM, Bowtie2) |
| `n_total_arg` | Total ARG-carrying contigs (geNomad analysis) |
| `n_mobile_arg` | Plasmid-associated ARG contigs (geNomad score ≥ 0.7) |
| `n_chrom_arg` | Chromosomal ARG contigs |
| `pct_mobile` | % of ARG contigs that are plasmid-associated |

Used by: Fig 1, Fig 3, Fig 4, Fig 5, `compute_arg_rpkm.py`

---

### `arg_contigs.tsv` — 791 rows × 16 columns

One row per ARG-carrying contig per sample. Consolidates contig taxonomy, plasmid classification, and mapping RPKM.

| Column | Description |
|--------|-------------|
| `sample_id` | Patient identifier |
| `contig` | metaSPAdes contig name |
| `status` | Contig assembly status |
| `taxid` | NCBI taxonomy ID (KrakenUnique) |
| `organism` | Organism name |
| `is_PA` | Boolean: contig assigned to *P. aeruginosa* |
| `broad_group` | Organism group (e.g. "Pseudomonas aeruginosa", "Haemophilus spp.") |
| `PA_pct` | Patient-level PA% (repeated per contig for convenience) |
| `on_plasmid` | Boolean: geNomad plasmid score ≥ 0.7 |
| `plasmid_score` | geNomad plasmid probability (0–1) |
| `n_hallmarks` | geNomad plasmid hallmark gene count |
| `conjugation_genes` | geNomad conjugation gene annotations |
| `length` | Contig length (bp) |
| `mapped_reads` | Reads mapped to contig (Bowtie2) |
| `library_size_M` | Sample library size in millions (for RPKM) |
| `rpkm` | Reads per kilobase per million mapped reads |

Used by: Fig 3, Fig 4, `compute_arg_rpkm.py`, `summarize_idxstats.py`

---

### `merged_aro_per_sample.tsv` — 1278 rows

One row per (sample, ARO accession). Union of RGI STRICT and AMRFinderPlus detections, harmonised via argNorm. Columns include `ARO_id`, `gene_label`, `drug_class`, `is_mac`, `detected_by`.

Used by: Fig 1, Fig 5, `compute_arg_rpkm.py`

---

### `aro_prevalence.tsv` — 227 rows

One row per unique ARO accession detected in the cohort. Columns: `ARO_id`, `gene_label`, `broad_class`, `Resistance Mechanism`, `n_samples` (number of patients with detection). Sorted by prevalence descending.

Used by: Fig 1

---

### `resistome_class_per_sample.tsv` — 40 rows

Wide-format pivot: one row per sample, one column per drug class (unique ARO count). Also contains `PA_pct` and `AZM`. Used directly for the Fig 1 heatmap.

Used by: Fig 1

---

### `assembly_arg_richness_taxonomy.tsv` — 40 rows

Wide-format pivot: one row per sample, one column per organism group (unique ARO count attributed to that group). Also contains `n_total`, `n_PA`, `n_nonPA`, `PA_pct`, `AZM`, `Leeds`.

Used by: Fig 2

---

### `macrolide_depth_hits.tsv` — 100 rows

One row per macrolide gene detection per sample. Columns: `sample_id`, `gene`, `Contig id`, `depth_contig`, `gene_len`, `pct_id`, `depth`, `reads_M`, `RPKM`. Used to assign RPKM values to individual macrolide genes.

Used by: Fig 3, Fig 4

---

## External data

| Dataset | Location |
|---------|----------|
| Raw metagenomic reads | ENA: [PRJEB111121](https://www.ebi.ac.uk/ena/browser/view/PRJEB111121) |
| `hamronization_combined_report.tsv` (~338 MB) | Zenodo DOI: [pending] |
| `depth.tsv` (~59 MB) | Zenodo DOI: [pending] |
