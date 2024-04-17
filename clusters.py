import hail as hl

def main():

        esse = 'hail mt file of Russian genotypes'
        relatives = 'IDs of esse relatives txt file'
        cluster = 'IDs of esse samples from cluster N txt file'

        hl.init(tmp_dir='tmp')

        mt = hl.read_matrix_table(f'{esse}')
        related_samples_to_remove = hl.import_table(f'{relatives}', impute=True, no_header=True)
        related_samples_to_remove =  related_samples_to_remove.key_by(related_samples_to_remove.f0)
        mt = mt.filter_cols(hl.is_defined(related_samples_to_remove[mt.col_key]), keep=False)
        for i in ['cl1','cl2','cl3','cl4','cl5','cl6']:
              cl1 = hl.import_table(f'{cluster}_{i}.txt', impute=True, no_header=True)
              cl1 =  cl1.key_by(cl1.f0)
              mt_1 = mt.filter_cols(hl.is_defined(cl1[mt.col_key]), keep=True)
              mt_1 = hl.variant_qc(mt_1)
              mt_1.variant_qc.AF.export(f'AF_RUS_{i}.txt')

if __name__=='__main__':

        main()
