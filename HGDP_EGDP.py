import hail as hl

def main():

        HGDP_EGDP = 'hail mt file of HGDP + EGDP genotypes'
        esse_G1K = 'hail mt file of Russians + 1000 Genomes genotypes'
        esse_G1K_HGDP_EGDP = 'merged hail mt file of  Russians + 1000 Genomes + HGDP + EGDP'
        pruned_variants = 'pruned variants of Russians + 1000 Genomes'
        esse_G1K_HGDP_EGDP_pruned = 'merged pruned hail mt file of  Russians + 1000 Genomes + HGDP + EGDP'
        outliers = 'txt file with IDs of esse outliers'
        G_annotations = 'txt file with annotations of 1000G populations'
        HGDP_samples = 'HGDP samples txt file'
        EGDP_samples = 'EGDP samples txt file'
        esse_G1K_HGDP_EGDP_pruned_pca = 'PCA merged pruned hail mt file of Russians + 1000 Genomes + HGDP + EGDP'

        hl.init(tmp_dir='tmp')
        rg = hl.get_reference('GRCh38')
        rg37 = hl.get_reference('GRCh37')
        recode = {f"{i}":f"chr{i}" for i in (list(range(1, 23)) + ['X', 'Y'])}
        recode['MT'] = 'chrM'

        rg37.add_liftover('references_grch37_to_grch38.over.chain.gz', rg)
        #intersection with HGDP_EGDP
        HGDP_EGDP = hl.read_matrix_table(f'{HGDP_EGDP}')
        mt = hl.read_matrix_table(f'{esse_G1K}')
        mt_intersect = mt.filter_rows(hl.is_defined(HGDP_EGDP.rows()[mt.row_key]))
        HGDP_EGDP_intersect = HGDP_EGDP.filter_rows(hl.is_defined(mt.rows()[HGDP_EGDP.row_key]))
        ECCE_1000G = hl.experimental.full_outer_join_mt(mt_intersect,HGDP_EGDP_intersect)
        ECCE_1000G = ECCE_1000G.select_entries(GT=hl.or_else(ECCE_1000G.left_entry.GT, ECCE_1000G.right_entry.GT))
        ECCE_1000G.write(f'{esse_G1K_HGDP_EGDP}',overwrite = True)

        mt = hl.read_matrix_table(f'{esse_G1K_HGDP_EGDP}')
        rsids = hl.import_table(f'{pruned_variants}', impute=True, no_header=True)
        rsids = rsids.key_by(rsids.f0)
        rsids = rsids.f0.collect()
        rsids = hl.literal(rsids)
        mt_pruned = mt.filter_rows(rsids.contains(mt.left_row.left_row.rsid),keep=True)
        mt_pruned.write(f'{esse_G1K_HGDP_EGDP_pruned}')

        #PCA without outliers
        HGDP_EGDP = hl.read_matrix_table(f'{esse_G1K_HGDP_EGDP_pruned}')
        PCA_exclude = hl.import_table(f'{outliers}', impute=True, no_header=True)
        PCA_exclude = PCA_exclude.key_by(PCA_exclude.f0)
        mt = HGDP_EGDP.filter_cols(hl.is_defined(PCA_exclude[HGDP_EGDP.col_key]), keep=False)
        
        #annotation and selection
        sa_1000G = hl.import_table(f'{G_annotations}', impute=True, key='Sample')
        mt = mt.annotate_cols(pheno1 = sa_1000G[mt.s])
        mt = mt.annotate_cols(SUPERPOP = hl.if_else((mt.s.startswith('HG0') | mt.s.startswith('NA')), mt.pheno1.SuperPopulation, 'NaN'))
        mt = mt.annotate_cols(POP = hl.if_else((mt.s.startswith('HG0') | mt.s.startswith('NA')), mt.pheno1.Population, 'NaN'))
        mt = mt.annotate_cols(SUPERPOP = hl.if_else(regions_difinition)
        mt = mt.annotate_cols(POP = hl.if_else(regions_difinition)
        mt = mt.annotate_cols(POP = hl.if_else(regions_difinition)
        mt = mt.annotate_cols(POP = hl.if_else(regions_difinition)
        sa_1000G = hl.import_table(f'{HGDP_samples}', impute=True, key='ID')
        mt = mt.annotate_cols(pheno2 = sa_1000G[mt.s])
        mt = mt.annotate_cols(POP = hl.if_else(mt.s.startswith('HGDP'), mt.pheno2.Population, mt.POP))

        sa_1000G = hl.import_table(f'{EGDP_samples}', impute=True, key='ID')
        sa_1000G = sa_1000G.annotate(x = sa_1000G.ID + '_' + sa_1000G.ID)
        sa_1000G = sa_1000G.key_by(sa_1000G.x)
        mt = mt.annotate_cols(pheno3 = sa_1000G[mt.s])
        mt = mt.annotate_cols(POP = hl.if_else(mt.s.startswith('GS0'), mt.pheno3.Population, mt.POP))
        POPs = ['SAMARA','ORENBURG','SPB',
        #'IBS','TSI','CEU','GBR','FIN',
        "Abkhazians_EGDP","Armenians_EGDP","Avars_EGDP","Azerbaijanis_EGDP","Balkars_EGDP",
        "Circassians_EGDP","Georgians_EGDP","Kabardins_EGDP","Kumyks_EGDP","Lezgins_EGDP",
        "North_Ossetians_EGDP","Tabasarans_EGDP","Ishkasim_EGDP","Kazakhs_EGDP","Kyrgyz_EGDP",          
        "Kyrgyz_Tdj_EGDP","Rushan_Vanch_EGDP","Shugnan_EGDP","Tajiks_EGDP","Turkmens_EGDP" ,       
        "Uygurs_EGDP","Uzbek_EGDP","Yaghnobi_EGDP","Albanians_EGDP","Bashkirs_EGDP",        
        "Belarusians_EGDP","Chuvashes_EGDP","Cossacks_EGDP","Cossacks_Kuban_EGDP","Croats_EGDP",          
        "Estonians_EGDP","Finnish_EGDP" ,"Germans_EGDP","Hungarians_EGDP" ,"Ingrians_EGDP",        
        "Karelians_EGDP","Komis_EGDP","Kryashen_Tatars_EGDP","Latvians_EGDP","Lithuanians_EGDP",     
        "Maris_EGDP","Mishar_Tatars_EGDP","Moldavians_EGDP","Mordvins_EGDP","Poles_EGDP",           
        "Roma_EGDP","Russians_EGDP","Russians_Central_EGDP","Russians_North_EGDP","Russians_West_EGDP",   
        "Saami_EGDP","Swedes_EGDP","Tatars_EGDP","Udmurds_EGDP","Ukrainians_east_EGDP", 
        "Ukrainians_north_EGDP","Ukrainians_west_EGDP","Vepsas_EGDP",
        "Altaians_EGDP","Buryats_EGDP",         
        "Chukchis_EGDP","Eskimo_EGDP","Evenks_EGDP","Evens_Magadan_EGDP","Evens_Sakha_EGDP",     
        "Forest_Nenets_EGDP","Kets_EGDP","Khantys_EGDP","Koryaks_EGDP","Mansis_EGDP",          
        "Mongolians_EGDP","Nganasans_EGDP","Sakha_EGDP","Selkups_EGDP","Shor_EGDP",            
        "Tundra_Nenets_EGDP","Tuvinians_EGDP","Yakuts_EGDP",
        #'French_HGDP','Sardinian_HGDP','Orcadian_HGDP','Russian_HGDP','BergamoItalian_HGDP','Tuscan_HGDP','Basque_HGDP','Adygei_HGDP',
        'Russian_HGDP']
        #'Balochi_HGDP','Bedouin_HGDP','Brahui_HGDP','Burusho_HGDP','Cambodian_HGDP','Dai_HGDP','Daur_HGDP','Druze_HGDP','Han_HGDP','Hazara_HGDP','Hezhen_HGDP','Japanese_HGDP','Kalash_HGDP','Lahu_HGDP','Makrani_HGDP','Miao_HGDP','Mongolian_HGDP','Mozabite_HGDP','Naxi_HGDP','NorthernHan_HGDP','Oroqen_HGDP','Palestinian_HGDP','Pathan_HGDP','She_HGDP','Sindhi_HGDP','Tu_HGDP','Tujia_HGDP','Uygur_HGDP','Xibo_HGDP','Yakut_HGDP','Yi_HGDP']
        POPs = hl.literal(POPs)
        mt =mt.filter_cols(POPs.contains(mt.POP),keep=True)

        intervals = hl.parse_locus_interval('1-22',reference_genome=rg37)
        ECCE_1000G_exclude_auto = hl.filter_intervals(mt, intervals=[intervals], keep = True)
        pca_eigenvalues, pca_scores, pca_loadings = hl.hwe_normalized_pca(ECCE_1000G_exclude_auto.GT,compute_loadings = True,k=10)

        #pca_scores.write(f'{esse_G1K_HGDP_EGDP_pruned_pca}', overwrite = True)
        pca_scores.export(f'{esse_G1K_HGDP_EGDP_pruned_pca}.txt')
        pca_loadings.export(f'{esse_G1K_HGDP_EGDP_pruned_pca}_loadings.txt')
        #print(pca_eigenvalues)

if __name__=='__main__':

        main()
