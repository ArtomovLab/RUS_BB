import hail as hl
from bokeh.layouts import row, column
from bokeh.io import output_file
from bokeh.io import save

def main():

        esse = 'hail mt file of Russian genotypes'
        HRC_AF = 'hail table file of HRC AFs'
        RSIDS = 'hail table file of RSIDS'
        GNOMAD_AF = 'hail table file of gnomAD AFs FIN,NEFIN'
        variants_to_exclude = 'excluded variants txt file'
        MAF_RUS_FIN_NEFIN = 'txt file of Russian MAFs for all variants in the Russian dataset'
        ANNOTATION_GNOMAD_filters = 'txt file with all variants from gnomAD that do not pass filters'

        hl.init(tmp_dir='tmp')

        mt = hl.read_matrix_table(f'{esse}')
        af = hl.read_table(f'{HRC_AF}')
        mt = mt.annotate_rows(AF_HRC = af[mt.locus,mt.alleles].AF)
        mt = hl.sample_qc(mt)
        mt.write('{esse}_with_HRC_AFs')        

        rsid = hl.read_table(f'{RSIDS}')
        mt = mt.annotate_rows(RSID = rsid[mt.locus,mt.alleles].rsid)
        af = hl.read_table(f'{GNOMAD_AF}')
        mt = mt.annotate_rows(AF = af[mt.locus,mt.alleles])

        mt = hl.variant_qc(mt)
        mt = mt.annotate_rows(MAF_RUS = hl.min(mt.variant_qc.AF))
        mt = mt.annotate_rows(MAF_FIN = hl.if_else(mt.variant_qc.AF[1] < 0.5, mt.AF.FIN,1-mt.AF.FIN))
        mt = mt.annotate_rows(MAF_NEFIN = hl.if_else(mt.variant_qc.AF[1] < 0.5, mt.AF.NEFIN,1-mt.AF.NEFIN))
        m = mt.select_rows(mt.MAF_RUS,mt.RSID,mt.MAF_FIN,mt.MAF_NEFIN,mt.variant_qc.p_value_hwe)
     
        m.export(f'{GNOMAD_AF}_RUS')
      
        mt = mt.filter_rows(mt.info.IMP == False)
        mt = mt.filter_rows(mt.variant_qc.p_value_hwe > 0.0001)
        print(mt.count())
        p1 = hl.plot.scatter(mt.MAF_RUS, mt.MAF_NEFIN) 
        p2 = hl.plot.scatter(mt.MAF_RUS, mt.MAF_FIN)
        P3 = column([row([p1,p2])])
        output_file(f'AF_FIN_NEFIN.html')
        save(P3)
        
        mt = hl.read_matrix_table(f'{esse}')
        variants = pd.read_table(f'{variants_to_exclude}')
        variants = list(variants['rsid'])
        variants = hl.literal(variants)
        mt = mt.filter_rows(variants.contains(mt.rsid),keep=False)
        mt.rows().select('MAF_RUS').export(f'{MAF_RUS_FIN_NEFIN}')
        ##################################################################################################33
        af = hl.read_table(f'{ANNOTATION_GNOMAD_filters}')
        af = af.filter(af.call_rate >= 0.97)
        af = af.filter(hl.len(af.filters) == 0)

        mt = hl.read_matrix_table(f'{esse}')
        variants = pd.read_table(f'{variants_to_exclude}')
        variants = list(variants['rsid'])
        variants = hl.literal(variants)
        
        mt = mt.filter_rows(variants.contains(mt.rsid),keep=False)
        mt = mt.annotate_rows(GNOMAD_FILTERS = af[mt.locus,mt.alleles])
        mt = mt.filter_rows(hl.is_defined(mt.GNOMAD_FILTERS),keep=True)

        p1 = hl.plot.scatter(mt.variant_qc.AF[1], mt.AF_HRC, xlabel='AF_RUS', ylabel='AF_HRC')
        p2 = hl.plot.scatter(mt.MAF_RUS, mt.MAF_NEFIN, xlabel='MAF_RUS', ylabel='MAF_NFE')
        p3 = hl.plot.scatter(mt.MAF_RUS, mt.MAF_FIN, xlabel='MAF_RUS', ylabel='MAF_FIN')
        P3 = column([row([p1,p2,p3])])
        output_file(f'MAF_FIN_NEFIN_filtered_RF.html')
        save(P3)

        mt = mt.filter_rows(mt.info.IMP == False)
        p1 = hl.plot.scatter(mt.variant_qc.AF[1], mt.AF_HRC, xlabel='AF_RUS', ylabel='AF_HRC')
        p2 = hl.plot.scatter(mt.MAF_RUS, mt.MAF_NEFIN, xlabel='MAF_RUS', ylabel='MAF_NFE')
        p3 = hl.plot.scatter(mt.MAF_RUS, mt.MAF_FIN, xlabel='MAF_RUS', ylabel='MAF_FIN')
        P3 = column([row([p1,p2,p3])])
        output_file(f'MAF_FIN_NEFIN_filtered_RF_genotyped.html')
        save(P3)

if __name__=='__main__':

        main()
