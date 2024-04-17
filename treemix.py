import hail as hl

def main():

         G_annotations = 'annotations of 1000G populations txt file'
         esse_G1K = 'merged esse and 1000G hail mt file'
         relatives = 'IDs of relatives txt file'
         cluster = 'IDs of esse samples from cluster N txt file'
         AC <- 'path with ACs for each population'

         hl.init(tmp_dir='tmp')

         #FOR treemix
         mt = hl.read_matrix_table(esse_G1K)
         related_samples_to_remove = hl.import_table(relatives, impute=True, no_header=True)
         related_samples_to_remove =  related_samples_to_remove.key_by(related_samples_to_remove.f0)
         mt = mt.filter_cols(hl.is_defined(related_samples_to_remove[mt.col_key]), keep=False)

         sa_1000G = hl.import_table(G_annotations, impute=True, key='Sample')
         mt = mt.annotate_cols(pheno1 = sa_1000G[mt.s])
         mt = mt.annotate_cols(SUPERPOP = hl.if_else((mt.s.startswith('HG') | mt.s.startswith('NA')), mt.pheno1.SuperPopulation, 'NaN'))
         mt = mt.annotate_cols(POP = hl.if_else((mt.s.startswith('HG') | mt.s.startswith('NA')), mt.pheno1.Population, 'NaN'))
         mt = mt.annotate_cols(SUPERPOP = hl.if_else(definition_of_regions)
         mt = mt.annotate_cols(POP = hl.if_else(definition_of_regions)
         mt = mt.annotate_cols(POP = hl.if_else(definition_of_regions)
         mt = mt.annotate_cols(POP = hl.if_else(definition_of_regions)

         for i in [1,2,3,4,5,6]:
             cl1 = hl.import_table(f'{AC}/{cluster}_cl{i}.txt', impute=True, key='x')
             cl1 = cl1.annotate(test = 0)
             mt_f = mt.annotate_cols(test = cl1[mt.s])
             mt_f.test.show()
             mt_f = mt_f.filter_cols(hl.is_nan(mt_f.test.test),keep=False)
             mt_f = hl.variant_qc(mt_f)
             mt_f.variant_qc.AC.export(f'{AC}/cl{i}_all.txt')

         for i in set(mt.POP.collect()):
             mt_f = mt.filter_cols(mt.POP == i)
             mt_f = hl.variant_qc(mt_f)
             mt_f.variant_qc.AC.export(f'{AC}/{i}_all.txt')

if __name__=='__main__':

         main()
