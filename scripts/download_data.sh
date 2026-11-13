#!/bin/bash
# Download GenAge data into the data/ folder
set -e
cd "$(dirname "$0")/../data"

echo "Downloading human genes..."
curl -sL -o human_genes.zip https://genomics.senescence.info/genes/human_genes.zip
unzip -o human_genes.zip
rm -f human_genes.zip release.html

echo "Downloading model organism genes..."
curl -sL -o models_genes.zip https://genomics.senescence.info/genes/models_genes.zip
unzip -o models_genes.zip
rm -f models_genes.zip release.html

echo "Downloading human homologs..."
curl -sL -o human_homologs.csv "https://genomics.senescence.info/genes/search.php?organism=Homo+sapiens&show=5&export=true&homologs=true"

echo "Done. Files:"
ls -lh
