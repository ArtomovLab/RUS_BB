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

        mt = hl.read_matrix_table('/humgen/atgu1/methods/dusoltsev/biobank/HRC/esse.concat.hrc_qc_duprem_all.ht')
        related_samples_to_remove = hl.import_table('/humgen/atgu1/methods/dusoltsev/biobank/HRC/plink2.king.cutoff.out.id', impute=True, no_header=True)
        related_samples_to_remove =  related_samples_to_remove.key_by(related_samples_to_remove.f0)
        mt = mt.filter_cols(hl.is_defined(related_samples_to_remove[mt.col_key]), keep=False)
        for i in ['SPB','ORENBURG','SAMARA']:
              cl1 = hl.import_table(f'/humgen/atgu1/methods/dusoltsev/biobank/HRC/fst/samples_WGS_{i}.txt', impute=True, no_header=True)
              cl1 =  cl1.key_by(cl1.f0)
              mt_1 = mt.filter_cols(hl.is_defined(cl1[mt.col_key]), keep=True)
              mt_1 = hl.variant_qc(mt_1)
              mt_1.variant_qc.AF.export(f'/humgen/atgu1/methods/dusoltsev/biobank/HRC/AF_RUS_{i}.txt')

if __name__=='__main__':

        parser = argparse.ArgumentParser(
                formatter_class=argparse.ArgumentDefaultsHelpFormatter)

        parser.add_argument('-ph', type=int, default = 1, help='GWAS phenotype id')
        args = parser.parse_args()
        main(args)
