This directory is messy. The cause is that initially, the plan was to use the 29 mammals dataset, which gives TFBS motifs conserved in 29 mammals. However, we later
decided to use motif enrichment analysis of a dataset with all TFBS motifs in the whole human genome. This was done because MHC genes themselves are some of the most
rapidly evolving genes in the human genome. Hence, conservation was less important. Therefore, only the data and scripts in the folder HumanWG_motifs are important.
These use a dataset by Pouya Kheradpour and M Kellis (Kheradpour, P. & Kellis, M. Systematic discovery and characterization of regulatory motifs
in ENCODE TF binding experiments. Nucleic Acids Res. 42, 2976 2987 (2014).) In particular, HumanWG_motifs/Scripts/EnrichmentTFBSRefactored.R does the heavy lifting
of finding cross-validated motifs enriched in positive set genes. See HumanWG_motifs/Documentation for
information about inputs and workings of the TFBS enrichment scripts. Note that these scripts were run remotely on Linux pcs and were not publication-ready, i.e. there are
hardcoded directories and hence none of the scripts work out of the box.