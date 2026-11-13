#!/usr/bin/env python3
"""
GenAge analysis tutorial script.

1. Download the files manually (or use the URLs below) into a folder called genage_data/
   - https://genomics.senescence.info/genes/human_genes.zip  → extract to genage_human.csv
   - https://genomics.senescence.info/genes/models_genes.zip → extract to genage_models.csv
   - Human homologs: https://genomics.senescence.info/genes/search.php?organism=Homo+sapiens&show=5&export=true&homologs=true
     (save as human_homologs.tsv or .csv)

2. Then run this script.
"""

import pandas as pd
from pathlib import Path
from gprofiler import GProfiler

# ---------------------------------------------------------------------------
# Point to the folder containing the downloaded CSVs
# ---------------------------------------------------------------------------
data_dir = Path("data")          # change this if your files are elsewhere
# or data_dir = Path("genage_data")

human_path = data_dir / "genage_human.csv"
models_path = data_dir / "genage_models.csv"
homologs_path = data_dir / "human_homologs.csv"   # or .tsv

print("Loading data from", data_dir.resolve())
human = pd.read_csv(human_path)
models = pd.read_csv(models_path)
# The homologs file is tab-separated and has a trailing tab (extra empty column)
homologs = pd.read_csv(
    homologs_path,
    sep="\t",
    header=0,
    names=["Species", "Symbol", "Entrez_ID", "Model_Organism", "Model_Symbol", "Model_Entrez", "extra"],
    usecols=range(6),  # ignore the extra empty column
)

print(f"Human GenAge genes: {len(human)}")
print(f"Model organism genes: {len(models)}")
print(f"Human homolog mappings: {len(homologs)}")

# ---------------------------------------------------------------------------
# Focus on worm, fly, mouse
# ---------------------------------------------------------------------------
focus_orgs = [
    "Caenorhabditis elegans",
    "Drosophila melanogaster",
    "Mus musculus",
]
models_focus = models[models["organism"].isin(focus_orgs)].copy()
print(f"\nGenes in worm / fly / mouse: {len(models_focus)}")
print(models_focus["organism"].value_counts())

# ---------------------------------------------------------------------------
# Find human orthologs
# ---------------------------------------------------------------------------
homologs_renamed = homologs.rename(columns={
    "Model_Organism": "organism",
    "Model_Symbol": "symbol",
    "Symbol": "human_symbol",
    "Entrez_ID": "human_entrez",
})

merged = models_focus.merge(
    homologs_renamed[["organism", "symbol", "human_symbol", "human_entrez"]],
    on=["organism", "symbol"],
    how="left",
)

has_human = merged.dropna(subset=["human_symbol"])
unique_human_symbols = has_human["human_symbol"].unique().tolist()

print(f"\nModel genes with at least one human homolog: {has_human['symbol'].nunique()}")
print(f"Unique human genes (orthologs): {len(unique_human_symbols)}")

out_map = Path("results") / "model_genes_with_human_homologs.csv"
has_human.to_csv(out_map, index=False)
print(f"Saved detailed mapping → {out_map}")

# ---------------------------------------------------------------------------
# Pathway enrichment with g:Profiler
# ---------------------------------------------------------------------------
print("\nRunning g:Profiler on the unique human orthologs...")
gp = GProfiler(return_dataframe=True)

results = gp.profile(
    organism="hsapiens",
    query=unique_human_symbols,
    sources=["KEGG", "REAC", "WP", "GO:BP"],
    no_evidences=True,
    user_threshold=0.05,
    significance_threshold_method="g_SCS",
)

print(f"Significant terms: {len(results)}")

aging_keywords = [
    "insulin", "igf", "mtor", "longevity", "aging", "ageing",
    "ampk", "foxo", "growth hormone"
]

def is_aging_related(name):
    n = str(name).lower()
    return any(kw.lower() in n for kw in aging_keywords)

results["aging_related"] = results["name"].apply(is_aging_related)
aging_hits = results[results["aging_related"]].sort_values("p_value")

print("\n=== Top aging-related pathways ===")
cols = ["source", "native", "name", "p_value", "intersection_size", "term_size"]
print(aging_hits[cols].head(20).to_string(index=False))

results.to_csv(Path("results") / "pathway_enrichment_human_orthologs.csv", index=False)
print("\nFull results saved to results/pathway_enrichment_human_orthologs.csv")

print("\n--- Bonus: enrichment on the 307 curated human GenAge genes ---")
human_results = gp.profile(
    organism="hsapiens",
    query=human["symbol"].tolist(),
    sources=["KEGG", "REAC", "GO:BP"],
    no_evidences=True,
)
human_results["aging_related"] = human_results["name"].apply(is_aging_related)
print(human_results[human_results["aging_related"]]
      .sort_values("p_value")[["source", "name", "p_value"]].head(10).to_string(index=False))

human_results.to_csv(Path("results") / "human_genage_pathway_enrichment.csv", index=False)
print("Human GenAge enrichment saved.")

print("\nAll done!")
