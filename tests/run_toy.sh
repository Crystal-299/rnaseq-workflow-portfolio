#!/usr/bin/env bash
# Requires Linux/macOS and the tools in envs/rnaseq.yaml.
set -euo pipefail
cd "$(dirname "$0")/.."
python tests/make_toy_data.py
python workflow/scripts/validate_config.py tests/config.toy.yaml tests/samples.toy.tsv
snakemake --configfile tests/config.toy.yaml -n --cores 2
snakemake --configfile tests/config.toy.yaml --cores 2
python workflow/scripts/check_outputs.py tests/config.toy.yaml tests/samples.toy.tsv
echo "Synthetic upstream workflow passed. DESeq2 is an optional separate target."
