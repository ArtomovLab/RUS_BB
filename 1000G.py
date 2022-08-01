import hail as hl
import math
from bokeh.io import output_notebook, show
import pandas as pd
from bokeh.layouts import row, column
from bokeh.io import output_file
from bokeh.plotting import figure
from bokeh.io import save
import os
from bokeh.io import export_png
import argparse

def main(args):

        hl.init(tmp_dir='/broad/hptmp')
        rg = hl.get_reference('GRCh38')
        rg37 = hl.get_reference('GRCh37')
        recode = {f"{i}":f"chr{i}" for i in (list(range(1, 23)) + ['X', 'Y'])}
        recode['MT'] = 'chrM'
        rg37.add_liftover('/humgen/atgu1/methods/dusoltsev/biobank/references_grch37_to_grch38.over.chain.gz', rg)
        
        #Intersection with 1000G
        
        ECCE = hl.read_matrix_table('/humgen/atgu1/methods/dusoltsev/biobank/HRC/esse.concat.hrc_qc_duprem_all.ht')
        G1000 = hl.read_matrix_table(f'/humgen/atgu1/methods/dusoltsev/biobank/HRC/1000G_WGS/1000G_WGS_chr1_filtered.ht')
        for i in range(2,23):
             G1000_ = hl.read_matrix_table(f'/humgen/atgu1/methods/dusoltsev/biobank/HRC/1000G_WGS/1000G_WGS_chr{i}_filtered.ht')
             G1000 = G1000.union_rows(G1000_)
             #G1000 = G1000.filter_rows(hl.is_defined(ECCE.rows()[G1000.row_key]))
             #G1000 = hl.import_vcf(f'/humgen/atgu1/methods/nkolosov/all_1000G/vcf/ALL.chr{i}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz',force_bgz = True, reference_genome='GRCh37',skip_i$        
        G1000.write(f'/humgen/atgu1/methods/dusoltsev/biobank/HRC/1000G_WGS/1000G_WGS_all_filtered.ht',overwrite = True)
        
        G_1000 = hl.read_matrix_table('/humgen/atgu1/methods/dusoltsev/biobank/HRC/1000G_WGS/1000G_WGS_all_filtered.ht')
        mt = hl.read_matrix_table('/humgen/atgu1/methods/dusoltsev/biobank/HRC/esse.concat.hrc_qc_duprem_all.ht')

        mt_intersect = mt.filter_rows(hl.is_defined(G_1000.rows()[mt.row_key]))
        G_1000_intersect = G_1000.filter_rows(hl.is_defined(mt.rows()[G_1000.row_key]))
        
        ECCE_1000G = hl.experimental.full_outer_join_mt(mt_intersect,G_1000_intersect)
        ECCE_1000G = ECCE_1000G.select_entries(GT=hl.or_else(ECCE_1000G.left_entry.GT, ECCE_1000G.right_entry.GT))
        ECCE_1000G.write('/humgen/atgu1/methods/dusoltsev/biobank/HRC/esse.concat.hrc_qc_duprem_all_1000G_WGS.ht',overwrite = True)
        print(ECCE_1000G.count())
        
        #make VCF        
        mt = hl.read_matrix_table('/humgen/atgu1/methods/dusoltsev/biobank/HRC/esse.concat.hrc_qc_duprem_all.ht')
        hl.export_vcf(mt,f'/humgen/atgu1/methods/dusoltsev/biobank/HRC/esse.concat.hrc_qc_duprem_all1.vcf.bgz',tabix=True)
        
        #LD
        mt = hl.read_matrix_table('/humgen/atgu1/methods/dusoltsev/biobank/HRC/esse.concat.hrc_qc_duprem_all_1000G_WGS.ht')
        
        #EUR
        POP_keep = hl.import_table('/humgen/atgu1/methods/dusoltsev/biobank/HRC/fst/samples_WGS_RUS_EUR.txt', impute=True, no_header=True)
        POP_keep = POP_keep.key_by(POP_keep.f0)
        mt = mt.filter_cols(hl.is_defined(POP_keep[mt.col_key]), keep=True)
        mt = hl.sample_qc(mt)
        mt = hl.variant_qc(mt)
        mt = mt.filter_cols(mt.sample_qc.call_rate >= 0.9)
        mt = mt.filter_rows(mt.variant_qc.call_rate >= 0.9)
        mt = mt.filter_rows(hl.min(mt.variant_qc.AF) > 0.01)
        mt = mt.filter_rows(mt.variant_qc.p_value_hwe > 0.0001)
        mt.write('/humgen/atgu1/methods/dusoltsev/biobank/HRC/esse.concat.hrc_qc_duprem_all_1000G_WGS_EUR_postqc.ht',overwrite = True)
        
        mt = hl.read_matrix_table('/humgen/atgu1/methods/dusoltsev/biobank/HRC/esse.concat.hrc_qc_duprem_all_1000G_WGS_EUR_postqc.ht')
        hl.export_vcf(mt,f'/humgen/atgu1/methods/dusoltsev/biobank/HRC/esse.concat.hrc_qc_duprem_all_1000G_WGS_postqc.vcf.bgz',tabix=True)

        biallelic_mt = mt.filter_rows(hl.len(mt.alleles) == 2)
        pruned_variant_table = hl.ld_prune(biallelic_mt.GT,r2=0.2, bp_window_size = 250000)
        mt_pruned = mt.filter_rows(hl.is_defined(pruned_variant_table[mt.row_key]))
        mt_pruned.write('/humgen/atgu1/methods/dusoltsev/biobank/HRC/esse.concat.hrc_qc_duprem_all_pruned_1000G_WGS_EUR.ht')

        #PCA with outliers
        ECCE_1000G = hl.read_matrix_table('/humgen/atgu1/methods/dusoltsev/biobank/HRC/esse.concat.hrc_qc_duprem_all_pruned_1000G_WGS.ht')
        sa_1000G = hl.import_table('/humgen/atgu1/methods/dusoltsev/biobank/1kg_annotations.txt', impute=True, key='Sample')

        intervals = hl.parse_locus_interval('1-22',reference_genome=rg37)
        ECCE_1000G_auto = hl.filter_intervals(ECCE_1000G, intervals=[intervals], keep = True)
        pca_eigenvalues, pca_scores, pca_loadings = hl.hwe_normalized_pca(ECCE_1000G_auto.GT,compute_loadings = True)

        pca_scores.write('/humgen/atgu1/methods/dusoltsev/biobank/HRC/esse.concat.hrc_qc_duprem_all_pruned_1000G_WGS_EUR_pca.ht', overwrite = True)
        pca_scores.export('/humgen/atgu1/methods/dusoltsev/biobank/HRC/esse.concat.hrc_qc_duprem_all_pruned_1000G_WGS_EUR_pca.txt')
        pca_loadings.export('/humgen/atgu1/methods/dusoltsev/biobank/HRC/esse.concat.hrc_qc_duprem_all_pruned_1000G_WGS_EUR_pca_loadings.txt')
        print(pca_eigenvalues)

        #PCA without outliers
        PCA_exclude = hl.import_table('/humgen/atgu1/methods/dusoltsev/biobank/HRC/samples_to_remove2.txt', impute=True, no_header=True)
        PCA_exclude = PCA_exclude.key_by(PCA_exclude.f0)
        ECCE_1000G_exclude = ECCE_1000G.filter_cols(hl.is_defined(PCA_exclude[ECCE_1000G.col_key]), keep=False) 
        intervals = hl.parse_locus_interval('1-22',reference_genome=rg37)
        ECCE_1000G_exclude_auto = hl.filter_intervals(ECCE_1000G_exclude, intervals=[intervals], keep = True)
        pca_eigenvalues, pca_scores, pca_loadings = hl.hwe_normalized_pca(ECCE_1000G_exclude_auto.GT,compute_loadings = True)

        pca_scores.write('/humgen/atgu1/methods/dusoltsev/biobank/HRC/esse.concat.hrc_qc_duprem_all_pruned_1000G_WGS_pca_exclude.ht', overwrite = True)
        pca_scores.export('/humgen/atgu1/methods/dusoltsev/biobank/HRC/esse.concat.hrc_qc_duprem_all_pruned_1000G_WGS_pca_exclude.txt')
        pca_loadings.export('/humgen/atgu1/methods/dusoltsev/biobank/HRC/esse.concat.hrc_qc_duprem_all_pruned_1000G_WGS_pca_exclude_loadings.txt')
        print(pca_eigenvalues)       
        
        #Original PCA
        mt = hl.read_matrix_table('/humgen/atgu1/methods/dusoltsev/biobank/HRC/esse.concat.hrc_qc_duprem_all.ht')
        rsids = hl.import_table('/humgen/atgu1/methods/dusoltsev/biobank/HRC/esse.hrc.dr08_pruned_PLINK2_final.prune.in', impute=True, no_header=True)
        rsids = rsids.key_by(rsids.f0)
        rsids = rsids.f0.collect()
        rsids = hl.literal(rsids)
        mt_pruned = mt.filter_rows(rsids.contains(mt.rsid),keep=True)

        PCA_exclude = hl.import_table('/humgen/atgu1/methods/dusoltsev/biobank/HRC/samples_to_remove2.txt', impute=True, no_header=True)
        PCA_exclude = PCA_exclude.key_by(PCA_exclude.f0)

        intervals = hl.parse_locus_interval('1-22',reference_genome=rg37)
        mt_pruned_auto = hl.filter_intervals(mt_pruned, intervals=[intervals], keep = True)
        pca_eigenvalues, pca_scores, pca_loadings = hl.hwe_normalized_pca(mt_pruned_auto.GT,compute_loadings = True)

        pca_scores.write('/humgen/atgu1/methods/dusoltsev/biobank/HRC/esse.concat.hrc_qc_duprem_all_pruned_pca_exclude.ht')
        pca_scores.export('/humgen/atgu1/methods/dusoltsev/biobank/HRC/esse.concat.hrc_qc_duprem_all_pruned_pca_exclude.txt')
        pca_loadings.export('/humgen/atgu1/methods/dusoltsev/biobank/HRC/esse.concat.hrc_qc_duprem_all_pruned_pca_exclude_loadings.txt')
        print(pca_eigenvalues)
    

if __name__=='__main__':

        parser = argparse.ArgumentParser(
                formatter_class=argparse.ArgumentDefaultsHelpFormatter)

        parser.add_argument('-ph', type=int, default = 1, help='GWAS phenotype id')
        args = parser.parse_args()
        main(args)
