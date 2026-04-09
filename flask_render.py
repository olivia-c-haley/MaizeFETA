from flask import Flask, request, render_template, jsonify
import pandas as pd
from math import isnan
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests


app = Flask(__name__)

# ── B73 paths ────────────────────────────────────────────────────────────────
B73_DIR = "./data/B73"

B73_GENE_MODEL_FILE = f"{B73_DIR}/B73.tsv"

B73_ONTOLOGY_FILES = {
    "gomap":      f"{B73_DIR}/gomap_go_v5.tsv",
    "entapgo":    f"{B73_DIR}/entap_go_v5.tsv",
    "pannzergo":  f"{B73_DIR}/pannzer_go_v5.tsv",
    "goexpanded": f"{B73_DIR}/expanded_go.tsv",
    "entapkegg":  f"{B73_DIR}/entap_kegg_v5.tsv",
    "metacyc":    f"{B73_DIR}/corncyc_v5.tsv",
    "reactome":   f"{B73_DIR}/reactome_v5.tsv",
    "pannzerec":  f"{B73_DIR}/pannzer_enzyme_commission_v5.tsv",
    "e2p2ec":     f"{B73_DIR}/e2p2_enzyme_commission_v5.tsv",
    "ecexpanded": f"{B73_DIR}/expanded_ec.tsv",
    "pfam":       f"{B73_DIR}/pfam_v5.tsv",
    "scatac":     f"{B73_DIR}/scatac_v5.tsv",
    "loci":       f"{B73_DIR}/loci_v5.tsv",
    "wallace_traits": f"{B73_DIR}/wallace_traits_v5.tsv",
    "atlas_traits":   f"{B73_DIR}/atlas_traits_v5.tsv",
    "deeploc":    f"{B73_DIR}/deeploc_v5.tsv",
    "interproscan": f"{B73_DIR}/interproscan_v5.tsv",
}

# ── Non-B73 paths ─────────────────────────────────────────────────────────────
NON_B73_DIR = "./data/nonB73"

NON_B73_GENOTYPES = [
    "B97","CML52","CML69","CML103","CML228","CML247","CML277",
    "CML322","CML333","HP301","Il14H","Ki3","Ki11","Ky21",
    "M37W","M162W","Mo18W","Ms71","NC350","NC358",
    "Oh7B","Oh43","P39","Tx303","Tzi8"
]

# Filename templates for non-B73 — {genotype} is substituted at runtime
# All files live flat in NON_B73_DIR (no per-genotype subfolder)
NON_B73_ONTOLOGY_TEMPLATES = {
    "pannzergo":  "{genotype}_pannzer_go.tsv",
    "entapgo":    "{genotype}_entap_go.tsv",
    "goexpanded": "{genotype}_expanded_go.tsv",
    "pannzerec":  "{genotype}_pannzer_enzyme_commission.tsv",
    "pfam":       "{genotype}_pfam.tsv",
    "deeploc":    "{genotype}_deeploc.tsv",
    "interproscan": "{genotype}_interproscan.tsv",
}


# ── Genotype prefix lookup (loaded lazily from data/genotypes.txt) ───────────
_GENOTYPE_PREFIXES = None

def get_genotype_prefixes():
    """Load prefix→genotype mapping on first use. Returns {} if file missing."""
    global _GENOTYPE_PREFIXES
    if _GENOTYPE_PREFIXES is None:
        try:
            df = pd.read_csv("./data/genotypes.txt", sep="\t")
            _GENOTYPE_PREFIXES = dict(zip(df["PREFIX"], df["GENOTYPE"]))
        except Exception as e:
            print(f"WARNING: could not load genotypes.txt: {e}")
            _GENOTYPE_PREFIXES = {}
    return _GENOTYPE_PREFIXES

def detect_genotype(gene_list):
    """Detect genotype from the prefix of the first recognisable gene ID.
    Returns (genotype, prefix) or (None, None) if unrecognised."""
    prefixes = get_genotype_prefixes()
    for gene in gene_list:
        for prefix, genotype in prefixes.items():
            if gene.startswith(prefix):
                return genotype, prefix
    return None, None

def psauron_file_for(genotype):
    """Return the correct PSAURON file path for a given genotype."""
    if genotype == "B73":
        return f"{B73_DIR}/psauron_v5.tsv"
    else:
        return f"{NON_B73_DIR}/{genotype}_psauron.tsv"


def phylostratr_file_for(genotype):
    """Return the correct Phylostratr file path for a given genotype."""
    if genotype == "B73":
        return f"{B73_DIR}/phylostratr_v5.tsv"
    else:
        return f"{NON_B73_DIR}/{genotype}_phylostratr.tsv"


# ── Core enrichment logic ─────────────────────────────────────────────────────
def calculate_enrichment(gene_model_list, query_genes, background_genes):
    genes_in_pathway     = set(gene_model_list)
    genes_in_query       = set(query_genes)
    genes_not_in_pathway = set(background_genes) - genes_in_pathway
    genes_not_in_query   = set(background_genes) - genes_in_query

    contingency_table = [
        [len(genes_in_pathway & genes_in_query),     len(genes_in_pathway & genes_not_in_query)],
        [len(genes_not_in_pathway & genes_in_query), len(genes_not_in_pathway & genes_not_in_query)],
    ]

    observed      = len(genes_in_pathway & genes_in_query) / len(genes_in_query) if genes_in_query else 0
    expected      = len(genes_in_pathway) / len(background_genes) if background_genes else 0
    expected_hits = round(expected * len(genes_in_query))
    etype         = "Overenriched" if observed > expected else "Underenriched"

    _, p_value = fisher_exact(contingency_table)

    return {
        "p_value":               float(p_value),
        "type":                  etype,
        "observed":              observed,
        "expected":              expected,
        "expected_hits":         expected_hits,
        "hits":                  len(genes_in_pathway & genes_in_query),
        "genes":                 "; ".join(genes_in_pathway & genes_in_query),
    }


def run_enrichment(classification_file, query_genes, background_genes, p_value_type, p_value_cutoff):
    try:
        df = pd.read_csv(classification_file, sep="\t")
    except FileNotFoundError:
        raise FileNotFoundError(f"Annotation file not found: {classification_file}")
    results = []

    for _, row in df.iterrows():
        gene_models = row["GENE_MODEL"].split(";")
        if not set(gene_models) & set(query_genes):
            continue
        result = calculate_enrichment(gene_models, query_genes, background_genes)
        result["label"]                = str(row["LABEL"]) if pd.notna(row["LABEL"]) else ""
        result["name"]                 = str(row["NAME"])  if pd.notna(row["NAME"])  else ""
        result["total_genes_in_label"] = len(gene_models)
        results.append(result)

    if not results:
        return results

    # Multiple testing correction
    if p_value_type == "Bonferroni":
        n = len(results)
        for r in results:
            r["adj_p_value"] = float(min(r["p_value"] * n, 1.0))
    elif p_value_type == "Benjamini-Hochberg":
        pvals     = [r["p_value"] for r in results]
        corrected = multipletests(pvals, method="fdr_bh")[1]
        for r, p_corr in zip(results, corrected):
            r["adj_p_value"] = 1.0 if isnan(p_corr) else float(p_corr)
    else:
        # Fallback — treat unrecognised method as BH
        pvals     = [r["p_value"] for r in results]
        corrected = multipletests(pvals, method="fdr_bh")[1]
        for r, p_corr in zip(results, corrected):
            r["adj_p_value"] = 1.0 if isnan(p_corr) else float(p_corr)

    return [r for r in results if r["adj_p_value"] <= p_value_cutoff]


# ── Routes ────────────────────────────────────────────────────────────────────
@app.route("/", methods=["GET"])
def index():
    return render_template("index.html")


@app.route("/nonb73", methods=["GET"])
def nonb73():
    return render_template("nonb73.html")


@app.route("/rnaseq", methods=["GET"])
def rnaseq():
    return render_template("rnaseq.html")


@app.route("/rnaseq_studies", methods=["GET"])
def rnaseq_studies():
    import os, csv
    studies_file = "./data/rnaseq_studies.txt"
    try:
        studies = []
        try:
            fh = open(studies_file, newline="", encoding="utf-8")
        except UnicodeDecodeError:
            fh = open(studies_file, newline="", encoding="latin-1")
        with fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                studies.append({
                    "filename":     row.get("Filename",     "").strip(),
                    "display_name": row.get("Display Name", "").strip(),
                    "genotype":     row.get("Genotype",     "").strip(),
                    "curation":     row.get("Curation Set", "").strip(),
                    "tissue":       row.get("Tissue",       "").strip(),
                    "condition":    row.get("Condition",    "").strip(),
                })
        return jsonify({"studies": studies})
    except FileNotFoundError:
        return jsonify({"error": f"rnaseq_studies.txt not found at {studies_file}"}), 404
    except Exception as e:
        return jsonify({"error": str(e)}), 500


@app.route("/rnaseq_file/<path:filename>", methods=["GET"])
def rnaseq_file(filename):
    import os
    # Sanitise — only allow filenames with no directory separators
    safe = os.path.basename(filename)
    filepath = os.path.join("./data/rnaseq", safe)
    if not os.path.isfile(filepath):
        return jsonify({"error": f"File not found: {safe}"}), 404
    try:
        with open(filepath, encoding="utf-8") as f:
            content = f.read()
    except UnicodeDecodeError:
        with open(filepath, encoding="latin-1") as f:
            content = f.read()
    return content, 200, {"Content-Type": "text/plain; charset=utf-8"}


@app.route("/analytics", methods=["GET"])
def analytics():
    return render_template("analytics.html")

@app.route("/faq", methods=["GET"])
def faq():
    return render_template("faq.html")

@app.route("/data-sources", methods=["GET"])
def data_sources():
    return render_template("data_sources.html")


@app.route("/detect_genotype", methods=["POST"])
def detect_genotype_route():
    try:
        genes_text  = request.form.get("genes", "").strip()
        gene_list   = [g.strip().replace("\r", "") for g in genes_text.split("\n") if g.strip()]
        if not gene_list:
            return jsonify({"error": "No genes provided."}), 400
        genotype, prefix = detect_genotype(gene_list)
        if genotype is None:
            return jsonify({"error": "Could not detect genotype. Please verify your gene model IDs."}), 400
        return jsonify({"genotype": genotype, "prefix": prefix})
    except Exception as e:
        import traceback
        traceback.print_exc()
        return jsonify({"error": f"Server error in genotype detection: {str(e)}"}), 500


@app.route("/psauron", methods=["POST"])
def psauron():
    import numpy as np
    import traceback
    from scipy.stats import mannwhitneyu, ks_2samp, chi2_contingency, gaussian_kde

    try:
        genes_text  = request.form.get("genes", "").strip()
        query_genes = [g.strip().replace("\r", "") for g in genes_text.split("\n") if g.strip()]
        if not query_genes:
            return jsonify({"error": "No genes provided."}), 400

        genotype, prefix = detect_genotype(query_genes)
        if genotype is None:
            return jsonify({"error": "Could not detect genotype from gene IDs. Please check that your gene models use a recognised NAM prefix (e.g. Zm00001eb for B73, Zm00018ab for B97)."}), 400

        psauron_file = psauron_file_for(genotype)
        try:
            df = pd.read_csv(psauron_file, sep="\t")
        except FileNotFoundError:
            return jsonify({"error": f"PSAURON data file not found for {genotype}: {psauron_file}"}), 500

        query_set  = set(query_genes)
        pop_df     = df.copy()
        query_df   = df[df["GENE_MODEL"].isin(query_set)]

        pop_scores    = pop_df["SCORE"].dropna().values
        query_scores  = query_df["SCORE"].dropna().values
        query_missing = int(len(query_set) - len(query_df))

        if len(query_scores) < 2:
            return jsonify({"error": f"Fewer than 2 query genes found in PSAURON data for {genotype}. Please check your gene IDs."}), 400

        # Compute KDE server-side on fixed [0, 1] grid
        x_grid      = np.linspace(0, 1, 300)
        pop_kde_y   = gaussian_kde(pop_scores,   bw_method="silverman")(x_grid).tolist()
        query_kde_y = gaussian_kde(query_scores, bw_method="silverman")(x_grid).tolist()

        pop_true_n   = int(pop_df["NAME"].astype(bool).sum())
        query_true_n = int(query_df["NAME"].astype(bool).sum())

        mw_stat, mw_p = mannwhitneyu(query_scores, pop_scores, alternative="two-sided")
        ks_stat, ks_p = ks_2samp(query_scores, pop_scores)

        q_true  = query_true_n
        q_false = int(len(query_scores)) - q_true
        p_true  = pop_true_n
        p_false = int(len(pop_scores)) - p_true
        contingency = [[q_true, q_false], [p_true - q_true, p_false - q_false]]
        try:
            _, chi2_p, _, _ = chi2_contingency(contingency)
        except Exception:
            chi2_p = 1.0

        return jsonify({
            "genotype":      genotype,
            "x_grid":        x_grid.tolist(),
            "pop_kde_y":     pop_kde_y,
            "query_kde_y":   query_kde_y,
            "pop_n":         int(len(pop_scores)),
            "query_n":       int(len(query_scores)),
            "query_missing": query_missing,
            "stats": {
                "query_mean":     float(np.mean(query_scores)),
                "pop_mean":       float(np.mean(pop_scores)),
                "query_median":   float(np.median(query_scores)),
                "pop_median":     float(np.median(pop_scores)),
                "query_true_pct": q_true / len(query_scores) if len(query_scores) else 0,
                "pop_true_pct":   p_true / len(pop_scores)   if len(pop_scores)   else 0,
                "mw_pvalue":      float(mw_p),
                "ks_pvalue":      float(ks_p),
                "chi2_pvalue":    float(chi2_p),
            }
        })

    except Exception as e:
        traceback.print_exc()
        return jsonify({"error": f"Server error: {str(e)}"}), 500


@app.route("/phylostratr", methods=["POST"])
def phylostratr():
    import numpy as np
    import traceback
    from scipy.stats import chi2_contingency, fisher_exact
    from statsmodels.stats.multitest import multipletests

    try:
        genes_text  = request.form.get("genes", "").strip()
        query_genes = [g.strip().replace("\r", "") for g in genes_text.split("\n") if g.strip()]
        if not query_genes:
            return jsonify({"error": "No genes provided."}), 400

        genotype, prefix = detect_genotype(query_genes)
        if genotype is None:
            return jsonify({"error": "Could not detect genotype from gene IDs."}), 400

        phyfile = phylostratr_file_for(genotype)
        try:
            df = pd.read_csv(phyfile, sep="\t")
        except FileNotFoundError:
            return jsonify({"error": f"Phylostratr file not found for {genotype}: {phyfile}"}), 500

        df = df.dropna(subset=["LABEL"])
        df["LABEL"] = df["LABEL"].astype(float)

        strata = (df[["LABEL", "NAME"]]
                  .drop_duplicates()
                  .sort_values("LABEL")
                  .values.tolist())

        query_set     = set(query_genes)
        query_df      = df[df["GENE_MODEL"].isin(query_set)]
        pop_total     = len(df)
        query_total   = len(query_df)
        query_missing = len(query_set) - query_total

        pop_counts   = df.groupby("LABEL").size().to_dict()
        query_counts = query_df.groupby("LABEL").size().to_dict()

        # Per-stratum Fisher's exact test
        raw_pvals = []
        strata_out = []
        for label, name in strata:
            pc = int(pop_counts.get(label, 0))
            qc = int(query_counts.get(label, 0))
            # 2x2: in-stratum vs out-of-stratum, query vs population
            table = [
                [qc,               query_total - qc],
                [pc - qc,          (pop_total - query_total) - (pc - qc)],
            ]
            try:
                _, p = fisher_exact(table)
            except Exception:
                p = 1.0
            raw_pvals.append(p)
            direction = "over" if (qc / query_total if query_total else 0) > (pc / pop_total if pop_total else 0) else "under"
            strata_out.append({
                "label":       label,
                "name":        name,
                "pop_count":   pc,
                "query_count": qc,
                "pop_pct":     round(pc / pop_total   * 100, 2) if pop_total   else 0,
                "query_pct":   round(qc / query_total * 100, 2) if query_total else 0,
                "direction":   direction,
                "raw_pvalue":  float(p),
            })

        # BH correction across strata
        if raw_pvals:
            adj = multipletests(raw_pvals, method="fdr_bh")[1]
            for s, ap in zip(strata_out, adj):
                s["adj_pvalue"] = float(ap)
        else:
            for s in strata_out:
                s["adj_pvalue"] = 1.0

        # Global chi-squared
        contingency = [
            [s["query_count"] for s in strata_out],
            [s["pop_count"]   for s in strata_out],
        ]
        try:
            _, chi2_p, _, _ = chi2_contingency(contingency)
        except Exception:
            chi2_p = 1.0

        return jsonify({
            "genotype":      genotype,
            "strata":        strata_out,
            "pop_total":     pop_total,
            "query_total":   query_total,
            "query_missing": query_missing,
            "chi2_pvalue":   float(chi2_p),
        })

    except Exception as e:
        traceback.print_exc()
        return jsonify({"error": f"Server error: {str(e)}"}), 500


@app.route("/validate_genes", methods=["POST"])
def validate_genes():
    import traceback
    try:
        genes_text  = request.form.get("genes", "").strip()
        query_genes = [g.strip().replace("\r", "") for g in genes_text.split("\n") if g.strip()]
        if not query_genes:
            return jsonify({"error": "No genes provided."}), 400

        genotype, prefix = detect_genotype(query_genes)

        # Load background gene list to validate against
        if genotype == "B73":
            bg_file = f"{B73_DIR}/background_genes.txt"
        elif genotype is not None:
            bg_file = f"{NON_B73_DIR}/{genotype}_background.txt"
        else:
            bg_file = None

        known_genes = set()
        if bg_file:
            try:
                with open(bg_file) as f:
                    known_genes = {l.strip() for l in f if l.strip()}
            except FileNotFoundError:
                # Fall back to any ontology file to get known gene models
                pass

        # If we still have no background, try loading from a reference ontology file
        if not known_genes and genotype:
            ref_file = (f"{B73_DIR}/deeploc_v5.tsv" if genotype == "B73"
                        else f"{NON_B73_DIR}/{genotype}_deeploc.tsv")
            try:
                import pandas as pd
                df = pd.read_csv(ref_file, sep="\t")
                for row in df["GENE_MODEL"]:
                    known_genes.update(str(row).split(";"))
            except Exception:
                pass

        not_found = [g for g in query_genes if g not in known_genes] if known_genes else query_genes

        return jsonify({
            "genotype": genotype,
            "total": len(query_genes),
            "found": len(query_genes) - len(not_found),
            "not_found": not_found,
        })
    except Exception as e:
        traceback.print_exc()
        return jsonify({"error": str(e)}), 500


@app.route("/upload", methods=["POST"])
def upload():
    print("B73 upload route called")
    for key, value in request.form.items():
        print(f"  {key}: {value}")

    genes_text  = request.form.get("genes", "").strip()
    query_genes = [g.strip().replace("\r", "") for g in genes_text.split("\n") if g.strip()]
    if not query_genes:
        return jsonify({"error": "No genes provided."}), 400

    p_value_type   = request.form.get("evaluation", "Benjamini-Hochberg")
    p_value_cutoff = float(request.form.get("pvalue", 0.05))

    gene_model_df    = pd.read_csv(B73_GENE_MODEL_FILE, sep="\t")
    background_genes = set(gene_model_df["GENE_MODEL"].values)

    selected_ontologies = {}
    for key in ["go-ontology", "pathway-ontology", "ec-ontology", "domain-ontology",
                "scatac-ontology", "loci-ontology", "traits-ontology", "deeploc-ontology"]:
        ontology_name = request.form.get(key)
        if ontology_name and ontology_name in B73_ONTOLOGY_FILES:
            selected_ontologies[ontology_name] = B73_ONTOLOGY_FILES[ontology_name]

    if not selected_ontologies:
        return jsonify({"error": "No ontologies selected."}), 400

    all_results = {}
    for ontology_name, filepath in selected_ontologies.items():
        try:
            all_results[ontology_name] = run_enrichment(
                filepath, query_genes, background_genes, p_value_type, p_value_cutoff
            )
        except FileNotFoundError as e:
            return jsonify({"error": str(e)}), 500
        except Exception as e:
            import traceback; traceback.print_exc()
            return jsonify({"error": f"Error running enrichment for {ontology_name}: {str(e)}"}), 500

    return jsonify(all_results)


@app.route("/upload_nonb73", methods=["POST"])
def upload_nonb73():
    print("Non-B73 upload route called")
    for key, value in request.form.items():
        print(f"  {key}: {value}")

    background = request.form.get("background", "").strip()
    if background not in NON_B73_GENOTYPES:
        return jsonify({"error": f"Unknown background: '{background}'."}), 400

    genes_text  = request.form.get("genes", "").strip()
    query_genes = [g.strip().replace("\r", "") for g in genes_text.split("\n") if g.strip()]
    if not query_genes:
        return jsonify({"error": "No genes provided."}), 400

    p_value_type   = request.form.get("evaluation", "Benjamini-Hochberg")
    p_value_cutoff = float(request.form.get("pvalue", 0.05))

    bg_file = f"{NON_B73_DIR}/{background}.tsv"
    try:
        gene_model_df    = pd.read_csv(bg_file, sep="\t")
        background_genes = set(gene_model_df["GENE_MODEL"].values)
    except FileNotFoundError:
        return jsonify({"error": f"Background file not found: {bg_file}"}), 500

    selected_ontologies = {}
    for key in ["go-ontology", "ec-ontology", "domain-ontology", "deeploc-ontology"]:
        ontology_name = request.form.get(key)
        if ontology_name and ontology_name in NON_B73_ONTOLOGY_TEMPLATES:
            filename = NON_B73_ONTOLOGY_TEMPLATES[ontology_name].format(genotype=background)
            selected_ontologies[ontology_name] = f"{NON_B73_DIR}/{filename}"

    if not selected_ontologies:
        return jsonify({"error": "No ontologies selected."}), 400

    all_results = {}
    for ontology_name, filepath in selected_ontologies.items():
        try:
            all_results[ontology_name] = run_enrichment(
                filepath, query_genes, background_genes, p_value_type, p_value_cutoff
            )
        except FileNotFoundError as e:
            return jsonify({"error": str(e)}), 500
        except Exception as e:
            import traceback; traceback.print_exc()
            return jsonify({"error": f"Error running enrichment for {ontology_name}: {str(e)}"}), 500

    return jsonify(all_results)


if __name__ == "__main__":
    app.run(debug=True)
