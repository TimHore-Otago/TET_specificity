# TET_specificity
This repository contains the scripts and datasets for the paper Pronounced sequence specificity guides cellular function of TET enzymes, Ravichandran <em>et al</em>.

## Components of the analysis

### Scripts
- Python script to quantitate methylation by hexamer (call_hexamers.py).
- Python script to sort methylation according to CGI location (classify_CGI.py).
- Python script to merge results and obtain summary tables (summarize.py).
- R Script to perform the Intra-motif positional preference analysis and plot the results  in R (IMPP_analysis.R).
- R Scripts to determine sequence logos using linear regression for logarithmic rates (seqlogo_NgTET_weighted_in_vitro.html, seqlogo_mTET[1-3]_weighted_in_vitro.Rmd, TET3_in_vivo.Rmd)

### Datasets
- CpG Illingworth bed file converted to GRCm38 (CpG_Il_mm10.bed).
- Demethylation velocity table (Slope_summary.txt).

Last update: 12:57 26/01/2021 EST
