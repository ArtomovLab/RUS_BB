# **Russian Biobank paper**
Here is list of scripts which is used for pilot study of Russian Biobank

## **Quality control**

1. `SEX_masmatch_paper.R` collects sex mismatch
2. `qc_plink1.sh` performs QC for Russian samples
3. `relatives_paper.R` makes relatives plots

## **PCA analysis and clustering**

4. `1000G.sh and 1000G.py` merge Russian samples with 1000 Genomes samples and perform PCA
5. `pca_filter.R` performs filters to Russian samples PCA
6. `ECCE_pca_paper.R` makes PCA plots
7. `clustering.R` runs clustering analysis of the Russian population

## **Populational structure analysis**

8. `fst.sh` runs Fst analysis
9. `fst_plot_paper.R` makes Fst plots
10. `beagle.sh` runs IBD analysis
11. `IBD_merge.sh` merges IBD segments
12. `IBD_unite.R` collects agregated IBD data
13. `IBD_hraphs.R` makes median IBD heatmap
14. `admixture1.sh` runs ADMIXTURE analysis
15. `admixture_graph.R` makes plot of ADMIXTURE analysis

## **GWAS and construction PheWeb database**

16. `ECCE_gwas.sh` and ECCE_gwas.py run GWAS analysis
17. `pheweb.sh` reformats summary statistics for PheWeb
18. `pheweb_pipline.txt` list of steps for running PheWeb
19. `vep.parse.sh and vep.sh` run VEP analysis
20. `CSQ.R` parses VEP annotation and collects the most severe effect
21. `pheweb.R` merges VEP and PheWeb annotations
22. `LDSC.nf` calculates genetic correlations for pheweb platform (modified from https://github.com/statgen/pheweb-rg-pipeline) 
23. `postgap.sh, postgap_unite.R, gprior.sh and GPrior_paper.R` run gene prioritization
24. `Replication_paper.R` collects Finnish PheWas for Finnish enriched variants

## **Finnish and Russian Enriched variants**

25. `af.sh and af.py` collect alleles frequencies from gnomAD
26. `clusters.py` collects allele frequencies for each Russian cluster
27. `RUS_variants_paper.R` makes plots for Finnish and Russian enriched variants`
28. `IBDne.sh` calculates a populational size
29. `IBDne.R` makes plot of a populational size
30. `treemix.R, treemix.sh and treemix.py` run TreeMix analysis
