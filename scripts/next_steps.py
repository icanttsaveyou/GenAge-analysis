#!/usr/bin/env python3
"""
Next-steps analysis for GenAge results.
Does all 5 things:
1. Filter to classic aging pathways
2. Show which genes drive insulin/mTOR/FoxO
3. Compare worm vs fly vs mouse
4. Highlight well-known aging genes
5. Make a simple bar plot of top pathways
"""

import pandas as pd
from pathlib import Path
from gprofiler import GProfiler
try:
    import matplotlib.pyplot as plt
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False
    print("matplotlib not installed — skipping plot (run: pip install matplotlib)")

results_dir = Path("results")
data_dir = Path("data")

# ---------------------------------------------------------------------------
# Load the files we already generated
# ---------------------------------------------------------------------------
print("Loading previous results...")
mapping = pd.read_csv(results_dir / "model_genes_with_human_homologs.csv")
enrichment = pd.read_csv(results_dir / "pathway_enrichment_human_orthologs.csv")

print(f"Mapping rows: {len(mapping)}")
print(f"Enrichment rows: {len(enrichment)}")

# ---------------------------------------------------------------------------
# 1. Filter to classic pathways only
# ---------------------------------------------------------------------------
print("\n" + "="*60)
print("1. CLASSIC AGING PATHWAYS")
print("="*60)

keywords = ["insulin", "igf", "mtor", "foxo", "ampk", "longevity", "growth hormone"]
mask = enrichment["name"].str.lower().apply(lambda x: any(k in x for k in keywords))
classic = enrichment[mask].sort_values("p_value")

print(classic[["source", "name", "p_value", "intersection_size"]].head(20).to_string(index=False))
classic.to_csv(results_dir / "classic_aging_pathways.csv", index=False)
print(f"\nSaved → results/classic_aging_pathways.csv")

# ---------------------------------------------------------------------------
# 2. Which human genes drive insulin / mTOR / FoxO
# ---------------------------------------------------------------------------
print("\n" + "="*60)
print("2. GENES DRIVING KEY PATHWAYS (re-running with evidences)")
print("="*60)

# Get unique human symbols
human_genes = mapping["human_symbol"].dropna().unique().tolist()

gp = GProfiler(return_dataframe=True)
detailed = gp.profile(
    organism="hsapiens",
    query=human_genes,
    sources=["KEGG", "REAC", "GO:BP"],
    no_evidences=False,          # this gives us the gene lists
    user_threshold=0.05,
)

# Keep only the strongest classic ones
key_terms = detailed[
    detailed["name"].str.contains("insulin|mTOR|FoxO|longevity|AMPK", case=False, na=False)
].sort_values("p_value").head(8)

for _, row in key_terms.iterrows():
    genes = row.get("intersections", [])
    if isinstance(genes, str):
        genes = eval(genes) if genes.startswith("[") else genes.split(",")
    print(f"\n{row['name']}  (p={row['p_value']:.2e})")
    print("  Genes:", ", ".join(genes[:12]) + ("..." if len(genes) > 12 else ""))

# ---------------------------------------------------------------------------
# 3. Compare worm vs fly vs mouse
# ---------------------------------------------------------------------------
print("\n" + "="*60)
print("3. WORM vs FLY vs MOUSE")
print("="*60)

org_map = {
    "Caenorhabditis elegans": "worm",
    "Drosophila melanogaster": "fly",
    "Mus musculus": "mouse",
}

for org_full, short in org_map.items():
    genes = mapping[mapping["organism"] == org_full]["human_symbol"].dropna().unique().tolist()
    print(f"\n--- {short.upper()} ({len(genes)} human orthologs) ---")
    if len(genes) < 5:
        print("  Too few genes, skipping enrichment")
        continue
    res = gp.profile(organism="hsapiens", query=genes, sources=["KEGG", "GO:BP"], no_evidences=True)
    top = res[res["name"].str.contains("insulin|mTOR|FoxO|longevity|AMPK|growth", case=False, na=False)]
    if len(top) > 0:
        print(top[["name", "p_value"]].head(5).to_string(index=False))
    else:
        print("  No classic aging pathways in top hits")

# ---------------------------------------------------------------------------
# 4. Look at well-known aging genes
# ---------------------------------------------------------------------------
print("\n" + "="*60)
print("4. WELL-KNOWN AGING GENES IN THE DATA")
print("="*60)

famous = ["daf-2", "daf-16", "age-1", "let-363", "akt-1", "akt-2",
          "InR", "foxo", "Akt1", "Mtor", "Igf1r", "Insr", "Foxo1", "Foxo3"]

found = mapping[mapping["symbol"].str.lower().isin([f.lower() for f in famous]) |
                mapping["human_symbol"].str.lower().isin([f.lower() for f in famous])]
print(found[["organism", "symbol", "human_symbol", "name"]].to_string(index=False))

# ---------------------------------------------------------------------------
# 5. Simple bar plot of top classic pathways
# ---------------------------------------------------------------------------
print("\n" + "="*60)
print("5. MAKING BAR PLOT")
print("="*60)

if not HAS_MATPLOTLIB:
    print("Skipped (matplotlib not available). Install with: pip install matplotlib")
else:
    plot_data = classic.head(12).copy()
    plot_data["-log10(p)"] = -plot_data["p_value"].apply(lambda p: __import__("math").log10(max(p, 1e-300)))
    plot_data = plot_data.sort_values("-log10(p)")

    plt.figure(figsize=(10, 6))
    plt.barh(plot_data["name"], plot_data["-log10(p)"], color="steelblue")
    plt.xlabel("-log10(p-value)")
    plt.title("Top classic aging pathways (human orthologs of worm/fly/mouse GenAge genes)")
    plt.tight_layout()
    plot_path = results_dir / "top_pathways_barplot.png"
    plt.savefig(plot_path, dpi=150)
    print(f"Plot saved → {plot_path}")

print("\nAll 5 steps finished.")
