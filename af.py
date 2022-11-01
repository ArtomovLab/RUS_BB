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

        mt = hl.read_matrix_table('/humgen/atgu1/methods/dusoltsev/biobank/HRC/esse.concat.hrc_qc_duprem.ht')
        af = hl.read_table('/humgen/atgu1/methods/dusoltsev/biobank/HRC/HRC.all.newid.frq.ht')
        mt = mt.annotate_rows(AF_HRC = af[mt.locus,mt.alleles].AF)
        mt = hl.sample_qc(mt)
        mt.write('/humgen/atgu1/methods/dusoltsev/biobank/HRC/esse.concat.hrc_qc_duprem_all.ht')        

        rsid = hl.read_table('/humgen/atgu1/methods/dusoltsev/biobank/ANNOTATION_RSIDS.ht')
        mt = mt.annotate_rows(RSID = rsid[mt.locus,mt.alleles].rsid)
        af = hl.read_table('/humgen/atgu1/methods/dusoltsev/biobank/ANNOTATION_FREQ.ht')
        mt = mt.annotate_rows(AF = af[mt.locus,mt.alleles])

        mt = hl.variant_qc(mt)
        mt = mt.annotate_rows(MAF_RUS = hl.min(mt.variant_qc.AF))
        mt = mt.annotate_rows(MAF_FIN = hl.if_else(mt.variant_qc.AF[1] < 0.5, mt.AF.FIN,1-mt.AF.FIN))
        mt = mt.annotate_rows(MAF_NEFIN = hl.if_else(mt.variant_qc.AF[1] < 0.5, mt.AF.NEFIN,1-mt.AF.NEFIN))
        m = mt.select_rows(mt.MAF_RUS,mt.RSID,mt.MAF_FIN,mt.MAF_NEFIN,mt.variant_qc.p_value_hwe)
     
        m.export('/humgen/atgu1/methods/dusoltsev/biobank/HRC/AF.txt')
      
        mt = mt.filter_rows(mt.info.IMP == False)
        mt = mt.filter_rows(mt.variant_qc.p_value_hwe > 0.0001)
        print(mt.count())
        p1 = hl.plot.scatter(mt.MAF_RUS, mt.MAF_NEFIN) 
        p2 = hl.plot.scatter(mt.MAF_RUS, mt.MAF_FIN)
        P3 = column([row([p1,p2])])
        output_file(f'/humgen/atgu1/methods/dusoltsev/biobank/HRC/AF_FIN_NEFIN_genotyped.html')
        save(P3)
        
        mt = hl.read_matrix_table('/humgen/atgu1/methods/dusoltsev/biobank/HRC/esse.concat.hrc_qc_duprem_all.ht')
        variants = pd.read_table('/humgen/atgu1/methods/dusoltsev/biobank/HRC/variants_to_exclude.txt')
        variants = list(variants['rsid'])
        variants = hl.literal(variants)
        mt = mt.filter_rows(variants.contains(mt.rsid),keep=False)
        mt.rows().select('MAF_RUS').export('/humgen/atgu1/methods/dusoltsev/biobank/HRC/MAF_RUS.txt')
        ##################################################################################################33
        af = hl.read_table('/humgen/atgu1/methods/dusoltsev/biobank/ANNOTATION_GNOMAD_filters.ht')
        af = af.filter(af.call_rate >= 0.97)
        af = af.filter(hl.len(af.filters) == 0)

        mt = hl.read_matrix_table('/humgen/atgu1/methods/dusoltsev/biobank/HRC/esse.concat.hrc_qc_duprem_all.ht')
        variants = pd.read_table('/humgen/atgu1/methods/dusoltsev/biobank/HRC/variants_to_exclude.txt')
        variants = list(variants['rsid'])
        variants = hl.literal(variants)
        
        mt = mt.filter_rows(variants.contains(mt.rsid),keep=False)
        mt = mt.annotate_rows(GNOMAD_FILTERS = af[mt.locus,mt.alleles])
        mt = mt.filter_rows(hl.is_defined(mt.GNOMAD_FILTERS),keep=True)

        p1 = hl.plot.scatter(mt.variant_qc.AF[1], mt.AF_HRC, xlabel='AF_RUS', ylabel='AF_HRC')
        p2 = hl.plot.scatter(mt.MAF_RUS, mt.MAF_NEFIN, xlabel='MAF_RUS', ylabel='MAF_NFE')
        p3 = hl.plot.scatter(mt.MAF_RUS, mt.MAF_FIN, xlabel='MAF_RUS', ylabel='MAF_FIN')
        P3 = column([row([p1,p2,p3])])
        output_file(f'/humgen/atgu1/methods/dusoltsev/biobank/HRC/MAF_FIN_NEFIN_filtered_RF.html')
        save(P3)

        mt = mt.filter_rows(mt.info.IMP == False)
        p1 = hl.plot.scatter(mt.variant_qc.AF[1], mt.AF_HRC, xlabel='AF_RUS', ylabel='AF_HRC')
        p2 = hl.plot.scatter(mt.MAF_RUS, mt.MAF_NEFIN, xlabel='MAF_RUS', ylabel='MAF_NFE')
        p3 = hl.plot.scatter(mt.MAF_RUS, mt.MAF_FIN, xlabel='MAF_RUS', ylabel='MAF_FIN')
        P3 = column([row([p1,p2,p3])])
        output_file(f'/humgen/atgu1/methods/dusoltsev/biobank/HRC/MAF_FIN_NEFIN_filtered_RF_genotyped.html')
        save(P3)
        ##af = hl.read_table('/broad/hptmp/dusoltsev/gnomad.genomes.r2.1.1.sites.ht')
        ##af.freq.export('/humgen/atgu1/methods/dusoltsev/biobank/HRC/test.txt')

if __name__=='__main__':

        parser = argparse.ArgumentParser(
                formatter_class=argparse.ArgumentDefaultsHelpFormatter)

        parser.add_argument('-ph', type=int, default = 1, help='GWAS phenotype id')
        args = parser.parse_args()
        main(args)
