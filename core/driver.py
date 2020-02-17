import os
import sys

from patientBKFProcessor import PatientProcessor
from reasoner import Reasoner

class Driver:
    def __init__(self,
                 tumor_rna_expression_file=None,
                 normal_rna_expression_file=None,
                 clinical_data_file=None,
                 drug_data_file=None,
                 radiation_data_file=None,
                 mutation_data_file=None,
                 fused_bkb=None):

        self.tumor_rna_expression_file = tumor_rna_expression_file
        self.normal_rna_expression_file = normal_rna_expression_file
        self.clinical_data_file = clinical_data_file
        self.drug_data_file = drug_data_file
        self.radiation_data_file = radiation_data_file
        self.mutation_data_file = mutation_data_file

        #-- Process Data if there is any.
        PP = PatientProcessor()
        PP.processPatientGeneData(mutation_data_file, tumor_rna_expression_file, 'None')
        PP.processClinicalData(clinical_data_file)

        self.pp = PP
        self.patients = PP.patients

        if fused_bkb is not None:
            self.reasoner = Reasoner(fused_bkb, self.patients)
        else:
            raise NotImplementedError


if __name__ == '__main__':
    driver = Driver(
        tumor_rna_expression_file='/home/public/data/ncats/data_drop_02-11-2020/rnaseq_fpkm_uq_primary_tumor.csv',
        clinical_data_file='/home/public/data/ncats/data_drop_02-11-2020/clinical.csv',
        mutation_data_file='/home/public/data/ncats/data_drop_02-11-2020/wxs.csv')

    print(driver.reasoner.patients)
    print('Complete')
