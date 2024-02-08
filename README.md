# WID SMK code repository

This code accompanies the manuscript '**Cigarette smoking and e-cigarette use induce shared DNA methylation changes linked to carcinogenesis**', by Herzog *et al.* (2024) and provides the necessary scripts to perform analyses and reproduce figures and display items.


## Guide to contents and structure

* `0-source`:
	+ Where applicable, this directory contains datasets and information required (e.g., lists of CpGs) for scripts used in this repository.
	+ This directory also contains functions that are used in the manuscript, e.g., for plotting.
* `1-analysis-pipeline`:
	+ This directory contains, in ordered number, the main files to reproduce the analyses and major computations. Note: scripts may be nested in subfolders (1-delta-beta) or provided as singular R files. 
* `2-markdown`: 
	+ This directory contains code to reproduce the main display items (`WID_SMK.Rmd`) or Supplementary Figures (`Supplementary_WID_SMK.Rmd`).

If you have any questions, please get in touch.

To compute the new smoking indices in your dataset, please take a look at the [WID.smk](https://github.com/chiaraherzog/WID.smk) package.