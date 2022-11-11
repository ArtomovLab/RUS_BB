import hail as hl
from bokeh.io import output_notebook, show
import pandas as pd
from bokeh.layouts import row, column
from bokeh.io import output_file
from bokeh.plotting import figure
from bokeh.io import save
import os
import math

from bokeh.io import export_png
from bokeh.plotting import figure
from bokeh.io import export_svgs
import svglib.svglib as svglib
from reportlab.graphics import renderPDF
from pathlib import Path

hl.init(tmp_dir='/broad/hptmp')

#FOR treemix
mt = hl.read_matrix_table('/humgen/atgu1/methods/dusoltsev/biobank/HRC/esse.concat.hrc_qc_duprem_all_1000G_WGS.ht/')
related_samples_to_remove = hl.import_table('/humgen/atgu1/methods/dusoltsev/biobank/HRC/all_relatives.king.cutoff.out.id', impute=True, no_header=True)
related_samples_to_remove =  related_samples_to_remove.key_by(related_samples_to_remove.f0)
mt = mt.filter_cols(hl.is_defined(related_samples_to_remove[mt.col_key]), keep=False)

sa_1000G = hl.import_table('/humgen/atgu1/methods/dusoltsev/biobank/1kg_annotations.txt', impute=True, key='Sample')
mt = mt.annotate_cols(pheno1 = sa_1000G[mt.s])
mt = mt.annotate_cols(SUPERPOP = hl.if_else((mt.s.startswith('HG') | mt.s.startswith('NA')), mt.pheno1.SuperPopulation, 'NaN'))
mt = mt.annotate_cols(POP = hl.if_else((mt.s.startswith('HG') | mt.s.startswith('NA')), mt.pheno1.Population, 'NaN'))
mt = mt.annotate_cols(SUPERPOP = hl.if_else((mt.s.startswith('36') | mt.s.startswith('40') | mt.s.startswith('53')), 'RUS', mt.SUPERPOP))
mt = mt.annotate_cols(POP = hl.if_else((mt.s.startswith('36')),'SAMARA', mt.POP))
mt = mt.annotate_cols(POP = hl.if_else((mt.s.startswith('40')) | (mt.s.startswith('KBL')),'SPB', mt.POP))
mt = mt.annotate_cols(POP = hl.if_else((mt.s.startswith('53')),'ORENBURG', mt.POP))

for i in [1,2,3,4,5,6]:
    cl1 = hl.import_table(f'/humgen/atgu1/methods/dusoltsev/biobank/HRC/fst/samples_WGS_cl{i}.txt', impute=True, key='x')
    cl1 = cl1.annotate(test = 0)
    mt_f = mt.annotate_cols(test = cl1[mt.s])
    mt_f.test.show()
    mt_f = mt_f.filter_cols(hl.is_nan(mt_f.test.test),keep=False)
    mt_f = hl.variant_qc(mt_f)
    mt_f.variant_qc.AC.export(f'/humgen/atgu1/methods/dusoltsev/biobank/HRC/treemix/AC/cl{i}_all.txt')

for i in set(mt.POP.collect()):
    mt_f = mt.filter_cols(mt.POP == i)
    mt_f = hl.variant_qc(mt_f)
    mt_f.variant_qc.AC.export(f'/humgen/atgu1/methods/dusoltsev/biobank/HRC/treemix/AC/{i}_all.txt')
