# **Russian Biobank paper**
List of scripts which is used in the Russian Biobank study. Preprint available at (https://www.biorxiv.org/content/10.1101/2023.03.23.534000v1)

## **Quality Control**

1. `SEX_masmatch_paper.R` collects sex mismatch
2. `qc_plink1.sh` performs QC for Russian samples
3. `relatives_paper.R` generates relatives plots
4. `map.py` generates map at Fig. 1A and Sup. Fig. S14

## **PCA Analysis and Clustering**

5. `1000G.py` merges Russian samples with 1000 Genomes samples and performs PCA
6. `pca_filter.R` performs filters to Russian samples PCA
7. `ECCE_pca_paper.R` creats PCA plots
8. `clustering.R` runs clustering analysis of the Russian population

## **Populational Structure Analysis**

9. `fst.sh` runs Fst analysis
10. `fst_plot_paper.R` generates Fst plots
11. `beagle.sh` runs IBD analysis
12. `IBD_merge.sh` merges IBD segments
13. `IBD_unite.R` collects agregated IBD data
14. `IBD_graphs.R` creates median IBD heatmap
15. `admixture1.sh` runs ADMIXTURE analysis
16. `admixture_graph.R` generates plot of ADMIXTURE analysis

## **GWAS and Construction PheWeb database**

17. `ECCE_gwas.py` runs GWAS analysis
18. `pheweb.sh` reformats summary statistics for PheWeb
19. `pheweb_pipline.txt` lists steps for running PheWeb
20. `vep.parse.sh and vep.sh` run VEP analysis
21. `CSQ.R` parses VEP annotation and collects the most severe effect
22. `pheweb.R` merges VEP and PheWeb annotations
23. `LDSC.nf` calculates genetic correlations for pheweb platform (modified from https://github.com/statgen/pheweb-rg-pipeline) 
24. `postgap.sh, postgap_unite.R, gprior.sh and GPrior_paper.R` run gene prioritization
25. `Replication_paper.R` collects Finnish PheWas for Finnish enriched variants

## **Finnish and Russian Enriched variants**

26. `af.py` collect alleles frequencies from gnomAD
27. `clusters.py` collects allele frequencies for each Russian cluster
28. `RUS_variants_paper.R` generates plots for Finnish and Russian enriched variants`
29. `IBDne.sh` calculates a populational size
30. `IBDne.R` generates a plot of a populational size
31. `HGDP_EGDP.py` runs PCA for 1000G, HGDP, EGDP, Russians dataset
32. `treemix.R, treemix.sh and treemix.py` run TreeMix analysis
