#!/usr/bin/env bash
sbatch \
  --mail-type=END \
  --mem-per-cpu=2000 \
  --time=120:00:00 \
  -o snake_dcv.out -e snake_dcv.err \
snakemake \
--profile profile_simple/ \
-s workflow/Snakefile_dcv \
--rerun-incomplete \
--keep-incomplete \
-pr \
--cores 200 \
--use-conda \
--latency-wait 30 \
--show-failed-logs \
"$@"
