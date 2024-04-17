import hail as hl

def main():

        esse = 'hail mt file of Russian genotypes'
        G1K  = 'hail mt file of 1000G genotypes by chromosomes'
        esse_G1K = 'merged esse and 1000G hail mt file'
        RUS_EUR = 'Russian and 1000G Europeans samples txt file'
        esse_G1K_EUR = 'merged esse and 1000G EUR hail mt file'
        esse_G1K_EUR_pruned = 'merged esse and 1000G EUR pruned hail mt file'
        G_annotations = 'annotations of 1000G populations txt file'
        esse_G1K_EUR_pruned_pca = 'PCA merged esse and 1000G EUR pruned hail mt file'
        outliers = 'IDs of esse outliers txt file'
        esse_G1K_EUR_pruned_pca_exclude = 'PCA merged esse and 1000G EUR pruned hail mt file, outliers excluded'
        pruned_variants = 'esse pruned variants txt file'
        esse_pca_exclude = 'PCA esse pruned hail mt file, outliers excluded'

        hl.init(tmp_dir='tmp')
        rg = hl.get_reference('GRCh38')
        rg37 = hl.get_reference('GRCh37')
        recode = {f"{i}":f"chr{i}" for i in (list(range(1, 23)) + ['X', 'Y'])}
        recode['MT'] = 'chrM'
        
        #Intersection with 1000G
        
        ECCE = hl.read_matrix_table(f'{esse}')
        G1000 = hl.read_matrix_table(f'{G1K}_chr1.ht')
        for i in range(2,23):
             G1000_ = hl.read_matrix_table(f'{G1K}_chr{i}.ht')
             G1000 = G1000.union_rows(G1000_)
             #G1000 = G1000.filter_rows(hl.is_defined(ECCE.rows()[G1000.row_key]))

        G1000.write(f'{G1K}_all.ht',overwrite = True)
        
        G_1000 = hl.read_matrix_table(f'{G1K}_all.ht')
        mt = hl.read_matrix_table(f'{esse}')

        mt_intersect = mt.filter_rows(hl.is_defined(G_1000.rows()[mt.row_key]))
        G_1000_intersect = G_1000.filter_rows(hl.is_defined(mt.rows()[G_1000.row_key]))
        
        ECCE_1000G = hl.experimental.full_outer_join_mt(mt_intersect,G_1000_intersect)
        ECCE_1000G = ECCE_1000G.select_entries(GT=hl.or_else(ECCE_1000G.left_entry.GT, ECCE_1000G.right_entry.GT))
        ECCE_1000G.write(f'{esse_G1K}',overwrite = True)
        print(ECCE_1000G.count())
        
        #make VCF        
        mt = hl.read_matrix_table(f'{esse}')
        hl.export_vcf(mt,f'{esse}.vcf.bgz',tabix=True)
        
        #LD
        mt = hl.read_matrix_table(f'{esse_G1K}')
        
        #EUR
        POP_keep = hl.import_table(f'{RUS_EUR}', impute=True, no_header=True)
        POP_keep = POP_keep.key_by(POP_keep.f0)
        mt = mt.filter_cols(hl.is_defined(POP_keep[mt.col_key]), keep=True)
        mt = hl.sample_qc(mt)
        mt = hl.variant_qc(mt)
        mt = mt.filter_cols(mt.sample_qc.call_rate >= 0.9)
        mt = mt.filter_rows(mt.variant_qc.call_rate >= 0.9)
        mt = mt.filter_rows(hl.min(mt.variant_qc.AF) > 0.01)
        mt = mt.filter_rows(mt.variant_qc.p_value_hwe > 0.0001)
        mt.write(f'{esse_G1K_EUR}',overwrite = True)
        
        mt = hl.read_matrix_table(f'{esse_G1K_EUR}')
        hl.export_vcf(mt,f'{esse_G1K_EUR}.vcf.bgz',tabix=True)

        biallelic_mt = mt.filter_rows(hl.len(mt.alleles) == 2)
        pruned_variant_table = hl.ld_prune(biallelic_mt.GT,r2=0.2, bp_window_size = 250000)
        mt_pruned = mt.filter_rows(hl.is_defined(pruned_variant_table[mt.row_key]))
        mt_pruned.write(f'{esse_G1K_EUR_pruned}')

        #PCA with outliers
        ECCE_1000G = hl.read_matrix_table(f'{esse_G1K_EUR_pruned}')
        sa_1000G = hl.import_table(f'{G_annotations}', impute=True, key='Sample')

        intervals = hl.parse_locus_interval('1-22',reference_genome=rg37)
        ECCE_1000G_auto = hl.filter_intervals(ECCE_1000G, intervals=[intervals], keep = True)
        pca_eigenvalues, pca_scores, pca_loadings = hl.hwe_normalized_pca(ECCE_1000G_auto.GT,compute_loadings = True)

        pca_scores.write(f'{esse_G1K_EUR_pruned_pca}', overwrite = True)
        pca_scores.export(f'{esse_G1K_EUR_pruned_pca}.txt')
        pca_loadings.export(f'{esse_G1K_EUR_pruned_pca}_loadings.txt')
        print(pca_eigenvalues)

        #PCA without outliers
        PCA_exclude = hl.import_table(f'{outliers}', impute=True, no_header=True)
        PCA_exclude = PCA_exclude.key_by(PCA_exclude.f0)
        ECCE_1000G_exclude = ECCE_1000G.filter_cols(hl.is_defined(PCA_exclude[ECCE_1000G.col_key]), keep=False) 
        intervals = hl.parse_locus_interval('1-22',reference_genome=rg37)
        ECCE_1000G_exclude_auto = hl.filter_intervals(ECCE_1000G_exclude, intervals=[intervals], keep = True)
        pca_eigenvalues, pca_scores, pca_loadings = hl.hwe_normalized_pca(ECCE_1000G_exclude_auto.GT,compute_loadings = True)

        pca_scores.write(f'{esse_G1K_EUR_pruned_pca_exclude}', overwrite = True)
        pca_scores.export(f'{esse_G1K_EUR_pruned_pca_exclude}.txt')
        pca_loadings.export(f'{esse_G1K_EUR_pruned_pca_exclude}_loadings.txt')
        print(pca_eigenvalues)       
        
        #Original PCA
        mt = hl.read_matrix_table(f'{esse}')
        rsids = hl.import_table(f'{pruned_variants}', impute=True, no_header=True)
        rsids = rsids.key_by(rsids.f0)
        rsids = rsids.f0.collect()
        rsids = hl.literal(rsids)
        mt_pruned = mt.filter_rows(rsids.contains(mt.rsid),keep=True)

        PCA_exclude = hl.import_table(f'{outliers}', impute=True, no_header=True)
        PCA_exclude = PCA_exclude.key_by(PCA_exclude.f0)

        intervals = hl.parse_locus_interval('1-22',reference_genome=rg37)
        mt_pruned_auto = hl.filter_intervals(mt_pruned, intervals=[intervals], keep = True)
        pca_eigenvalues, pca_scores, pca_loadings = hl.hwe_normalized_pca(mt_pruned_auto.GT,compute_loadings = True)

        pca_scores.write(f'{esse_pca_exclude}')
        pca_scores.export(f'{esse_pca_exclude}.txt')
        pca_loadings.export(f'{esse_pca_exclude}_loadings.txt')
        print(pca_eigenvalues)
 
if __name__=='__main__':

        main()
