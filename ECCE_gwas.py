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

#(7142159, 4589)
#(6816775, 4572)

def main(args):
        hl.init(tmp_dir='/broad/hptmp')
        hl.set_global_seed(0)
        rg = hl.get_reference('GRCh38')
        rg37 = hl.get_reference('GRCh37')
        recode = {f"{i}":f"chr{i}" for i in (list(range(1, 23)) + ['X', 'Y'])}
        recode['MT'] = 'chrM'
        rg37.add_liftover('/humgen/atgu1/methods/dusoltsev/biobank/references_grch37_to_grch38.over.chain.gz', rg)

        #read geno
        mt = hl.read_matrix_table('/humgen/atgu1/methods/dusoltsev/biobank/HRC/esse.concat.hrc_qc_duprem_all.ht')

        sa = hl.import_table('/humgen/atgu1/methods/dusoltsev/biobank/Phen_gwas_f.txt',key='ID',impute=True)
        pca = hl.import_table('/humgen/atgu1/methods/dusoltsev/biobank/HRC/esse.hrc.dr08_pruned_PLINK2_final.eigenvec',key='IID',impute=True)
        related_samples_to_remove = hl.import_table('/humgen/atgu1/methods/dusoltsev/biobank/HRC/plink2.king.cutoff.out.id', impute=True, no_header=True)
        
        mt = mt.annotate_cols(pca = pca[mt.s])
        mt = mt.annotate_cols(pheno = sa[mt.s])
        related_samples_to_remove =  related_samples_to_remove.key_by(related_samples_to_remove.f0)
        mt = mt.filter_cols(hl.is_defined(related_samples_to_remove[mt.col_key]), keep=False)
        mt = mt.filter_cols(hl.is_defined(mt.pca),keep=True) 
        
        variants = pd.read_table('/humgen/atgu1/methods/dusoltsev/biobank/HRC/variants_to_exclude.txt')
        variants = list(variants['rsid'])
        variants = hl.literal(variants)
        mt = mt.filter_rows(variants.contains(mt.rsid),keep=False)

        df = pd.read_csv(f'{args.list}', sep ='\t')
        cols = df['pheno'].tolist()
        mods = df['mod'].tolist()
        covs = df['cov'].tolist()
        mod = mods[args.ph-1]
        cov = covs[args.ph-1]
        pheno = cols[args.ph-1]
        print(f'GWAS {pheno} start with mod {mod} amd add cov {cov}')
        col = f'{args.ph}'
        
        mt = mt.filter_cols(hl.is_nan(eval(f'mt.pheno.{pheno}')),keep=False)
        mt = mt.filter_cols(hl.is_defined(mt.pca),keep=True)
        
        if f'{pheno}' in ['BMI_under','BMI_over','BMI_obes']:
                 mt = mt.filter_cols((eval(f'mt.pheno.{pheno}') == 1) | (mt.pheno.BMI_norm == 1))
        if f'{pheno}' in ['PHT']:
                 mt = mt.filter_cols((eval(f'mt.pheno.{pheno}') == 1) | (mt.pheno.OptimBP == 1))
        if f'{pheno}' in ['HyperGLU12']:
                 mt = mt.filter_cols((eval(f'mt.pheno.{pheno}') == 1) | (mt.pheno.HyperGLU1 == 0))
       
        mt = hl.variant_qc(mt)
        mt = mt.filter_rows(hl.min(mt.variant_qc.AF) > 0.01)
        mt = mt.filter_rows(mt.variant_qc.p_value_hwe > 0.0001)
        Nsamp = mt.cols().count()
        
        if not pd.isnull(cov) and int(mod) == 1:
                  mt = mt.filter_cols(hl.is_nan(eval(f'mt.pheno.{cov}')),keep=False)
                  gwas = hl.linear_regression_rows(y=eval(f'mt.pheno.{pheno}'),
                  x=mt.GT.n_alt_alleles(),
                  covariates=[1.0,
                  mt.pheno.AGE,
                  mt.pheno.SEX,
                  mt.pca.PC1,
                  mt.pca.PC2,
                  mt.pca.PC3,
                  mt.pca.PC4 ] + [eval(f'mt.pheno.{cov}')], pass_through=['RSID',mt.variant_qc.AF])
                  gwas = gwas.filter(hl.is_nan(gwas.beta), keep = False)
        if pd.isnull(cov) and int(mod) == 1:       
                  gwas = hl.linear_regression_rows(y=eval(f'mt.pheno.{pheno}'),
                  x=mt.GT.n_alt_alleles(),
                  covariates=[1.0,
                  mt.pheno.AGE,
                  mt.pheno.SEX,
                  mt.pca.PC1,
                  mt.pca.PC2,
                  mt.pca.PC3,
                  mt.pca.PC4 ], pass_through=['RSID',mt.variant_qc.AF])
                  gwas = gwas.filter(hl.is_nan(gwas.beta), keep = False)
        if not pd.isnull(cov) and int(mod) == 2:

                  f = mt.aggregate_cols(hl.agg.counter(eval(f'mt.pheno.{pheno}')))
                  control = f[0]
                  case = f[1]

                  gwas = hl.logistic_regression_rows(test='wald',y=eval(f'mt.pheno.{pheno}'),
                  x=mt.GT.n_alt_alleles(),
                  covariates=[1.0,
                  mt.pheno.AGE,
                  mt.pheno.SEX,
                  mt.pca.PC1,
                  mt.pca.PC2,
                  mt.pca.PC3,
                  mt.pca.PC4 ] + [eval(f'mt.pheno.{cov}')], pass_through=['RSID',mt.variant_qc.AF])
                  gwas = gwas.filter(hl.is_nan(gwas.beta), keep = False)
        if pd.isnull(cov) and int(mod) == 2:
                  
                  f = mt.aggregate_cols(hl.agg.counter(eval(f'mt.pheno.{pheno}')))
                  control = f[0]
                  case = f[1]

                  gwas = hl.logistic_regression_rows(test='wald',y=eval(f'mt.pheno.{pheno}'),
                  x=mt.GT.n_alt_alleles(),
                  covariates=[1.0,
                  mt.pheno.AGE,
                  mt.pheno.SEX,
                  mt.pca.PC1,
                  mt.pca.PC2,
                  mt.pca.PC3,
                  mt.pca.PC4], pass_through=['RSID',mt.variant_qc.AF])
                  gwas = gwas.filter(hl.is_nan(gwas.beta), keep = False)
        if int(mod) == 2:
                  gwas = gwas.annotate(n=int(Nsamp))
                  gwas = gwas.annotate(control=int(control))
                  gwas = gwas.annotate(case=int(case))
                  gwas = gwas.select('RSID','AF','n','control','case','beta','standard_error','z_stat','p_value')
                  Nsamp = str(Nsamp) + '_logistic'        
        
        gwas.write(f'/humgen/atgu1/methods/dusoltsev/biobank/HRC/gwas_results_all/{pheno}_{Nsamp}.ht',overwrite=True)
        print(f'GWAS {pheno} hail written')
        gwas = hl.read_table(f'/humgen/atgu1/methods/dusoltsev/biobank/HRC/gwas_results_all/{pheno}_{Nsamp}.ht')
        gwas.export(f'/humgen/atgu1/methods/dusoltsev/biobank/HRC/gwas_results_all/{pheno}_{Nsamp}.txt.bgz')
        print(f'GWAS {pheno} txt written')
        p = hl.plot.manhattan(gwas.p_value)
        save(p, filename=f'/humgen/atgu1/methods/dusoltsev/biobank/HRC/gwas_results_all/graphs/{pheno}_{Nsamp}_manhattan.html')
        print(f'GWAS {pheno} manhattan written')
        qq = hl.plot.qq(gwas.p_value)
        save(qq, filename=f'/humgen/atgu1/methods/dusoltsev/biobank/HRC/gwas_results_all/graphs/{pheno}_{Nsamp}_qq.html')
        print(f'GWAS {pheno} qq written')

if __name__=='__main__':

        parser = argparse.ArgumentParser(
                formatter_class=argparse.ArgumentDefaultsHelpFormatter)

        parser.add_argument('-ph', type=int, default = 1, help='GWAS phenotype id')
        parser.add_argument('-list', type=str, default = 'Pheno_for_gwas', help='path to txt file with list of pheno')
        #parser.add_argument('-mod', type=int, default = 1, help='1-linear,2=logistic')
        args = parser.parse_args()
        main(args)
