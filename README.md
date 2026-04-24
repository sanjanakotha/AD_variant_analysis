# AD Variant Analysis

A Python toolkit for analyzing single nucleotide variants (SNVs) in **human** transcription factor (TF) protein domains, including activation domains (ADs), DNA-binding domains (DBDs), repression domains (RDs), bifunctional domains (Bif), and intrinsically disordered regions (IDRs).

> **Note:** The CDS BED files (`data/cds_beds/`), protein sequences (`data/raw_files/proteins.dat`), and DNA transcript sequences (`data/raw_files/dna_transcripts.dat`) provided in this repository are specific to **human TFs**.

## Features

- Intersect genomic variant BED files with CDS regions
- Classify CDS SNVs as synonymous, non-synonymous, or nonsense
- Map protein domain coordinates to genomic BED format
- Intersect domain regions with variants to identify domain-specific mutations

## Requirements

- Python ≥ 3.8
- [bedtools](https://bedtools.readthedocs.io/) (must be installed separately and available in `PATH`)

## Installation

Install directly from GitHub:

```bash
pip install git+https://github.com/sanjanakotha/AD_variant_analysis.git
```

Or clone the repository and install locally:

```bash
git clone https://github.com/sanjanakotha/AD_variant_analysis.git
cd AD_variant_analysis
pip install .
```

Or in editable mode for development:

```bash
git clone https://github.com/sanjanakotha/AD_variant_analysis.git
cd AD_variant_analysis
pip install -e .
```

## Repository layout

```
data/
├── cds_beds/                          # CDS BED files (one per ENST ID, 1,559 human TFs)
├── domain_beds/                       # Pre-generated genomic domain BED files
├── raw_files/
│   ├── proteins.dat                   # Human TF protein sequences (pickle)
│   └── dna_transcripts.dat            # Human TF CDS nucleotide sequences (pickle)
├── all_TFs_table_proteins_with_IDR.txt  # TF domain annotation table (1,590 human TFs)
└── SFARI_TF_de_novo_variants_gpf.bed    # Example variants (117 SFARI de novo SNVs)

output/                                # Pipeline results (created at runtime)
├── intersections/
├── sorted/
├── classified/
└── domain_snvs/
```

## Input data formats

### TF domain annotation table

A tab-separated file with one row per transcript. Required columns:

| Column | Description |
|---|---|
| `uniprotID` | UniProt accession |
| `ENSG` | Ensembl gene ID |
| `ENST` | Ensembl transcript ID |
| `DBD_coords` | DNA-binding domain coordinates (e.g. `10-50,60-90`) |
| `AD_coords` | Activation domain coordinates |
| `RD_coords` | Repression domain coordinates |
| `Bif_coords` | Bifunctional domain coordinates |
| `IDR_coords` | Intrinsically disordered region coordinates |

Coordinates are 1-based amino acid positions formatted as comma-separated `start-end` ranges.

An example file (`data/all_TFs_table_proteins_with_IDR.txt`) covering 1,590 human TFs is included in the repository.

### Variants BED file

A 6-column tab-separated BED file (no header):

```
chrom  start  end  ref_allele  alt_allele  score
```

An example file (`data/SFARI_TF_de_novo_variants_gpf.bed`) with 117 de novo variants from the SFARI database is included.

### CDS BED files

One BED file per transcript (named `<ENST_ID>`), listing the CDS exon intervals with strand information in column 5. These are available from Ensembl BioMart or can be generated with tools such as `gffread`.

Human TF CDS BED files are provided in the `data/cds_beds/` directory (one file per ENST ID, covering 1,559 human TF transcripts).

### Reference sequence files (human TFs)

Two Python pickle files (`.dat`) containing reference sequences for human TFs are provided in `data/raw_files/`:

| File | Description |
|---|---|
| `data/raw_files/proteins.dat` | Protein (amino acid) sequences keyed by ENST ID |
| `data/raw_files/dna_transcripts.dat` | CDS nucleotide sequences keyed by ENST ID |

## Pipeline overview

```
TF domain annotation + CDS BED files
        │
        ▼
  map-domains                 ← step 0: map AA domain coords to genomic BED

CDS BED files  +  Variants BED
        │                │
        ▼                ▼
  intersect-variants          ← step 1: find variants in CDS regions
        │
        ▼
  classify-snvs               ← step 2: Syn / Non-Syn / Nonsense
        │
        ▼
Domain BED files ──► classify-domain-snvs  ← step 3: per-domain results
```

## CLI commands

All commands accept `--help` for full argument details.

### `map-domains`

Map protein domain amino acid coordinates to genomic BED format using CDS BED files.
This generates the domain BED files required by subsequent steps.

```bash
map-domains \
  --input            data/all_TFs_table_proteins_with_IDR.txt \
  --cds-directory    data/cds_beds/ \
  --output-directory data/domain_beds/
```

### `intersect-variants`

Intersect a variants BED file with CDS BED files. Automatically sorts inputs and caches sorted files for reuse.

```bash
intersect-variants \
  --cds-dir      data/cds_beds/ \
  --variants     data/SFARI_TF_de_novo_variants_gpf.bed \
  --output-dir   output/
```

Output structure:

```
output/
├── intersections/          # one .bed per transcript with overlapping variants
└── sorted/                 # sorted copies of inputs (reused on subsequent runs)
```

### `classify-snvs`

Classify CDS variants as synonymous, non-synonymous, or nonsense using protein and nucleotide reference sequences.

```bash
classify-snvs \
  --input-dir      output/intersections/ \
  --output-dir     output/classified/ \
  --proteins       data/raw_files/proteins.dat \
  --dna-transcripts data/raw_files/dna_transcripts.dat \
  --sorted-cds-dir  output/sorted/cds/
```

### `classify-domain-snvs`

Intersect classified SNVs with protein domain regions to identify which variants fall in specific domains.

```bash
classify-domain-snvs \
  --domain-dir         data/domain_beds/ \
  --classified-snv-dir output/classified/ \
  --output-dir         output/domain_snvs/ \
  --mapping            data/all_TFs_table_proteins_with_IDR.txt
```

### `intersect-domain-variants`

Intersect domain BED files with CDS variant BED files (without SNV classification).

```bash
intersect-domain-variants \
  --domain-dir      data/domain_beds/ \
  --cds-variants-dir output/intersections/ \
  --output-dir       output/domain_variants/ \
  --mapping          data/all_TFs_table_proteins_with_IDR.txt
```

Use `--domain-type DBD` (or `AD`, `RD`, `Bif`, `IDR`) to restrict to a single domain type.

## Python API

The CLI tools can also be imported and called from Python:

```python
from pathlib import Path
from AD_variant_analysis.intersect_variants import intersect_variants, sort_bed_files
from AD_variant_analysis.classify_snvs import SNVClassifier
from AD_variant_analysis.intersect_domains_variants import intersect_domains_with_variants
```

See the `example.ipynb` notebook for an interactive walkthrough.

## Example notebook

`example.ipynb` demonstrates:

1. Loading and inspecting the included TF annotation table
2. Loading and inspecting the included variants BED file
3. Running the full analysis pipeline using the CLI commands
4. Interpreting the output

## License

MIT
