# FETA — Functional Enrichment Tool for maize Annotations

FETA is a locally-hosted web application for functional enrichment analysis of maize gene lists. Paste a set of gene model IDs, select your annotation types, and FETA runs Fisher's Exact Tests across a wide range of functional databases to identify over- and under-enriched terms in your query relative to the whole-genome background.

---

## Features

- **B73 enrichment** across 17 annotation databases (GO, KEGG, CornCyc/MetaCyc, Reactome, Enzyme Commission, Pfam, scATAC, loci, traits, DeepLoc subcellular localization, InterProScan)
- **Non-B73 enrichment** for all 25 NAM founder genotypes (B97, CML52, CML69, CML103, CML228, CML247, CML277, CML322, CML333, HP301, Il14H, Ki3, Ki11, Ky21, M37W, M162W, Mo18W, Ms71, NC350, NC358, Oh7B, Oh43, P39, Tx303, Tzi8) across GO, Enzyme Commission, Pfam, DeepLoc, and InterProScan annotations
- **Multiple testing correction** — Benjamini-Hochberg (FDR) or Bonferroni, configurable p-value cutoff
- **Auto-detection** of genotype from gene model ID prefixes
- **Three visualizations** per ontology result:
  - Bar chart (−log₁₀ adjusted p-value, colored by enrichment direction)
  - Dot plot (gene ratio × term, dot size = hit count, dual green/orange color scales for over/under enrichment)
  - CNET plot (bipartite force-directed network of terms and genes)
- **Toggle** to show/hide underenriched terms independently in each visualization
- **Other Analytics** tab with PSAURON protein-coding score distributions and Phylostratr evolutionary strata analysis, including Mann-Whitney U, Kolmogorov-Smirnov, and Chi-squared tests
- **Download** results as TSV, PNG, or SVG

---

## Project Structure

```
feta/
├── flask_render.py          # Flask application and all API routes
├── requirements.txt
├── templates/
│   ├── index.html           # B73 enrichment page
│   ├── nonb73.html          # Non-B73 enrichment page
│   ├── analytics.html       # PSAURON & Phylostratr analytics page
│   ├── faq.html
│   └── data_sources.html
└── data/
    ├── genotypes.txt        # PREFIX → GENOTYPE mapping table
    ├── B73/
    │   ├── B73.tsv                        # Background gene model list
    │   ├── background_genes.txt
    │   ├── gomap_go_v5.tsv
    │   ├── entap_go_v5.tsv
    │   ├── pannzer_go_v5.tsv
    │   ├── expanded.tsv                   # GO expanded/propagated
    │   ├── entap_kegg_v5.tsv
    │   ├── corncyc_v5.tsv
    │   ├── reactome_v5.tsv
    │   ├── pannzer_enzyme_commission_v5.tsv
    │   ├── e2p2_enzyme_commission_v5.tsv
    │   ├── expanded_ec.tsv
    │   ├── pfam_v5.tsv
    │   ├── scatac_v5.tsv
    │   ├── loci_v5.tsv
    │   ├── wallace_traits_v5.tsv
    │   ├── atlas_traits_v5.tsv
    │   ├── deeploc_v5.tsv
    │   ├── interproscan_v5.tsv
    │   ├── psauron_v5.tsv
    │   └── phylostratr_v5.tsv
    └── nonB73/
        ├── {GENOTYPE}.tsv                           # Background gene model list per genotype
        ├── {GENOTYPE}_background.txt
        ├── {GENOTYPE}_pannzer_go.tsv
        ├── {GENOTYPE}_entap_go.tsv
        ├── {GENOTYPE}_expanded.tsv
        ├── {GENOTYPE}_pannzer_enzyme_commission.tsv
        ├── {GENOTYPE}_pfam.tsv
        ├── {GENOTYPE}_deeploc.tsv
        ├── {GENOTYPE}_interproscan.tsv
        ├── {GENOTYPE}_psauron.tsv
        └── {GENOTYPE}_phylostratr.tsv
```

All annotation files are tab-separated with at minimum the columns `GENE_MODEL`, `LABEL`, and `NAME`. `GENE_MODEL` fields may contain semicolon-delimited lists of gene IDs for terms shared across multiple models.

---

## Installation

### 1. Clone the repository

```bash
git clone https://github.com/olivia-c-haley/MaizeFETA.git
cd MaizeFETA
```

### 2. Create conda environment (recommended)

```bash
conda create --name feta_env python=3.11
conda activate feta_env
```

### 3. Install dependencies

```bash
pip install -r requirements.txt
```

### 4. Add annotation data

FETA requires pre-built annotation TSV files under `data/B73/` and `data/nonB73/`. These are not bundled in the repository due to their size. Place your data files following the directory structure above. The `data/genotypes.txt` file should be a tab-separated table with columns `PREFIX` and `GENOTYPE`, e.g.:

```
PREFIX          GENOTYPE
Zm00001eb       B73
Zm00018ab       B97
Zm00019ab       CML52
...
```

---

## Running the App

```bash
python flask_render.py
```

Then open your browser to [http://127.0.0.1:5000](http://127.0.0.1:5000).

By default Flask runs in debug mode (set in `flask_render.py`). 

---

## Usage

### B73 Enrichment

1. Navigate to **B73 Enrichment**
2. Paste your B73 gene model IDs (one per line, e.g. `Zm00001eb000001`)
3. Click **Load Genes** to validate your list against the background
4. Select one or more enrichment types (GO, KEGG, Pathways, etc.) and choose an annotation source per type
5. Set your p-value correction method and cutoff
6. Click **Run Enrichment**

### Non-B73 Enrichment

1. Navigate to **Non-B73 Enrichment**
2. Paste your gene model IDs — the genome is auto-detected from the ID prefix, and the dropdown will update automatically
3. Confirm or manually select the NAM line background
4. Follow the same steps as B73

### Other Analytics

1. Navigate to **Other Analytics**
2. Paste your gene list and select PSAURON, Phylostratr, or both
3. Click **Run Analysis**
4. Results include score distribution plots, summary statistics, and statistical test results, downloadable as PNG, SVG, or TSV

---


## Dependencies

| Package | Purpose |
|---------|---------|
| Flask | Web framework and routing |
| pandas | TSV parsing and data manipulation |
| scipy | Fisher's Exact Test, Mann-Whitney U, KS test, KDE, Chi-squared |
| statsmodels | Benjamini-Hochberg and Bonferroni multiple testing correction |
| numpy | Numerical arrays for KDE grid computation |

Frontend libraries are loaded from CDN and require no installation:
- [D3.js](https://d3js.org/) — dot plots and CNET force graphs
- [Chart.js](https://www.chartjs.org/) — bar charts and KDE line charts

---

## License

See `LICENSE` for details.
