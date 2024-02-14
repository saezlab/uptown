# UPTOWN: A comprehensive workflow embedding baseline network contextualisation methods
This repository is part of my Master's Thesis, conducted at the Institute for Computational Biomedicine, Universitätsklinikum Heidelberg, between 01.08.2023 and 15.02.2024. For a frozen version of this repo, please consult the appropiate release.

## Introduction
 UPTOWN (Unified PlaTform fOr netWork iNference) is a unified and standardised environment to run six well-established, baseline methods for context-specific network inference. All methods aim to connect a source node, origin of the perturbation, to a set of several downstream nodes containing molecular measurements of the perturbation effects. The six methods contain very elementary approaches, but they are, to some degree, represented in many widely used tools for network inference. The shortest-paths and all-paths modules are used to connect sources to targets, while PPR acts as a heat-diffusion-like algorithm to prioritise well-connected nodes in the solution networks. The combination of these modules, coupled with a sign consistency check, provides a good framework for network inference using baseline, well established methods. 

 [picture]

 ## Contents
 This repository contains the two main files from UPTOWN, evaluation.py and solver.py. In addition, it contains files necessary for the analyses detailed in my Master's Thesis. Below is a brief overview of each file:
* decryptm_downloader.py: Downloads the raw files from the PRIDE database. decryptm_downloader.sh is a wrapper to be run in a cluster.
* decryptm_fit.py: It takes the dose-response data points and fits a logistic model, extracted from the original DecryptM study. It provides a file containing all the fitted parameters and an optional file containing the data points in a long format.
* decryptm_exploranalysis.R: provides the proteomics_targets.tsv file to be used as measurements for the Solver, and ec50_info.tsv containing the ec50 for each detected protein (whose filtering pased quality controls)
* fragpipe folder: contains the necessary files to reproduce the processing of the raw files from DecryptM. Fragpipe needs to be downloaded. These files need to be in the same folder as the raw files. It needs a proteome (Fasta file), which can be donwloaded by Fragpipe itself. Please check their documentation https://fragpipe.nesvilab.org/.
* mthesis_env.yaml: contains the necessary dependencies to run UPTOWN.
* panacea_proc.R: provides panacea_sources.tsv (file that contains the drug targets) and panacea_offtargets.tsv (for evaluation)
* pkn_preprocessing.py: provides the context-agnostic PKN (network_collectri.sif) and a list of TF regulons to be used with an enrichment analysis (collectri.tsv).
* plots.R: code for the plots contained in the manuscript
* uptown_analysis_cluster: uses the two datasets to run network inference in parallel. Uptown_network_calc.nf is a wrapper for a cluster.
* uptown_evaluation_cluster: uses the two datasets to run network inference in parallel. Uptown_network_eval.nf is a wrapper for a cluster. 
* uptown_network.config: file containing the Nextflow config.

