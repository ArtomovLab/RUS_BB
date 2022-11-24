# **Russian Biobank paper**
Here is list of scripts which is used for pilot study of Russian Biobank

1. `1000G.sh and 1000G.py` merge Russian samples with 1000 Genomes samples and perform PCA
2. `CSQ.R` parses VEP annotation and collects the most severe effect
3. `ECCE_gwas.sh` and ECCE_gwas.py run GWAS analysis 
4. `ECCE_pca_paper.R` makes PCA plots
5. `gprior.sh GPrior_paper.R` run gene prioritization
6. `IBD_graphs.R` makes IBD heatmap 
7. `IBD_merge.sh` merges IBD segments
8. `IBD_unite.R` collects agregated IBD data
9. `IBDne.R` makes plot of populational size
10. `IBDne.sh` calculates populational size
11. `LDSC.nf` calculates genetic correlations for pheweb platform (modified from https://github.com/statgen/pheweb-rg-pipeline)
12. `RUS_variants_paper.R` makes plots for Finnish and Russian enriched variants
13. `Replication_paper.R` collects Finnish PheWas for Finnish enriched variants
14. `SEX_masmatch_paper.R` collects sex mismatch
15. `admixture1.sh` runs ADMIXTURE analysis
16. `admixture_graph.R` makes plot of ADMIXTURE analysis
17. `af.sh` and af.py collect alleles frequencies from gnomAD
18. `beagle.sh` runs IBD analysis
19. `clustering.R` runs clustering analysis of Russian population
20. `clusters.py` collects allele frequencies for each Russian cluster
21. `fst.sh` runs Fst analysis
22. `fst_plot_paper.R` makes Fst plots
23. `pca_filter.R` performs filters to Russian samples PCA
24. `pheweb.R` merges VEP and PheWeb annotations
25. `pheweb.sh` reformats summary statistics for PheWeb
26. `pheweb_pipline.txt` list of steps for running PheWeb
27. `qc_plink1.sh` perform QC for Russian samples
28. `relatives_paper.R` makes relatives plots
29. `treemix.R, treemix.sh and treemix.py` run TreeMix analysis
30. `vep.parse.sh and vep.sh` run VEP analysis

