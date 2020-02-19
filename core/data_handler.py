import pandas as pd
import os
import pickle
import tqdm
import argparse
import csv

from pybkb.core.cpp_base.fusion import fuse
from patientBKFProcessor import PatientProcessor
from pathwayBKFProcessor import PathwayProcessor
from reactomePathwayProcessor import ReactomePathwayProcessor

class DataDriver:
    def __init__(self, config_file):
        #-- Read in configuration from config file.
        self.config = dict()
        with open(config_file, 'r') as csv_file:
            reader = csv.reader(csv_file)
            for row in reader:
                if row[1] == '':
                    self.config[row[0]] = None
                else:
                    self.config[row[0]] = row[1]

        self.dataHandler = DataHandler(tumor_rna_expression_file=self.config['tumor_rna_expression_file'],
                                       normal_rna_expression_file=self.config['normal_rna_expression_file'],
                                       clinical_data_file=self.config['clinical_data_file'],
                                       drug_data_file=self.config['drug_data_file'],
                                       radiation_data_file=self.config['radiation_data_file'],
                                       mutation_data_file=self.config['mutation_data_file'],
                                       pathway_file=self.config['pathway_file'],
                                       reactome_g2r_file=self.config['reactome_g2r_file'],
                                       reactome_p2p_file=self.config['reactome_p2p_file'],
                                       working_dir=self.config['working_dir'])

        self.dataHandler.readData()

    def run(self, i=0):
        if i == 0:
            print('=' * 50)
            print("Welcome to the BKB Data Driver!\n\tLet's fuse some fragments...")
            input('Press any key to continue...')
        print('Choose the patients you wish to fuse:')
        for j, patient in enumerate(self.dataHandler.patientProcessor.patients):
            print('{})\t{}'.format(j, patient.patientID))
        print('{})\tAll patients.'.format(j+1))
        choice = input('Your choices (seperate by spaces): ')
        try:
            if int(choice) == j+1:
                patients = [pat for pat in self.dataHandler.patientProcessor.patients]
        except ValueError:
            choice_split = choice.split()
            patient_indices = [int(pat) for pat in choice_split]
            patients = [self.dataHandler.patientProcessor.patients[idx] for idx in patient_indices]

        while True:
            print('You choose to fuse over the following patients:')
            if len(patients) == len(self.dataHandler.patientProcessor.patients):
                print('\tAll patients.')
            else:
                for pat in patients:
                    print('\t{}'.format(pat.patientID))
            correct = input('Is that correct? ([y], n): ') or 'y'
            if correct == 'y':
                break
            if correct == 'n':
                self.run(i=i+1)
            else:
                print('Response not recognized.')

        while True:
            print('Do you want to fuse with Pathway data?')
            choice = input('Your choice ([y],n): ') or 'y'
            if choice == 'y':
                with_pathways = True
                break
            if choice == 'n':
                with_pathways = False
                break
            else:
                print('Response not recognized.')

        while True:
            print('Ready to fuse.')
            choice = input('Do you want to continue? ([y],n): ') or 'y'
            if choice == 'y':
                self.dataHandler.fuse(patients, with_pathways)
                break
            if choice == 'n':
                while True:
                    print('Would you like to restart or quit?')
                    choice = input('Your choice ([q],r): ')
                    if choice == 'q':
                        return
                    if choice == 'r':
                        self.run(i=i+1)
                    else:
                        print('Response not recognized.')
            else:
                print('Response not recognized.')
        return

class DataHandler:
    def __init__(self,
                 tumor_rna_expression_file=None,
                 normal_rna_expression_file=None,
                 clinical_data_file=None,
                 drug_data_file=None,
                 radiation_data_file=None,
                 mutation_data_file=None,
                 pathway_file=None,
                 reactome_p2p_file=None,
                 reactome_g2r_file=None,
                 working_dir='/tmp'):

        self.tumor_rna_expression_file = tumor_rna_expression_file
        self.normal_rna_expression_file = normal_rna_expression_file
        self.clinical_data_file = clinical_data_file
        self.drug_data_file = drug_data_file
        self.radiation_data_file = radiation_data_file
        self.mutation_data_file = mutation_data_file
        self.pathway_file = pathway_file
        self.reactome_p2p_file = reactome_p2p_file
        self.reactome_g2r_file = reactome_g2r_file
        self.working_dir = working_dir

    def readTcgaData(self):
        self.patientProcessor = PatientProcessor()

        if self.tumor_rna_expression_file is not None and self.mutation_data_file is not None:
            self.patientProcessor.processPatientGeneData(self.mutation_data_file,
                                                         self.tumor_rna_expression_file,
                                                         None)
        if self.clinical_data_file is not None:
            self.patientProcessor.processClinicalData(self.clinical_data_file)

        self.patientProcessor.processPatientBKF()

    def readPathwayData(self):
        self.pathwayProcessor = PathwayProcessor()
        if self.pathway_file is not None:
            self.pathwayProcessor.processPathways(self.pathway_file)
            self.pathwayProcessor.processPathwayBKF()
        else:
            print('Warning: No pathway directory given.')

        self.reactomePathwayProcessor = ReactomePathwayProcessor()
        if self.reactome_g2r_file is not None and self.reactome_p2p_file is not None:
            self.reactomePathwayProcessor.processGenePathways(self.reactome_g2r_file)
            self.reactomePathwayProcessor.processHierarchyPathways(self.reactome_p2p_file)
            self.reactomePathwayProcessor.processPathwayBKF()
        else:
            print('Warning: No reactome files given.')

    def readData(self):
        self.readTcgaData()
        self.readPathwayData()

    def fuse(self, patients, with_pathways=True):
        patient_indices = list()
        for patient_name in patients:
            idx = self.patientProcessor.patients.index(patient_name)
            patient_indices.append(idx)
        patient_bkf_files, patient_source_names, dump_loc = self.patientProcessor.SubsetBKFsToFile(self.working_dir, patient_indices)

        if with_pathways:
            pathway_bkf_files, pathway_source_names = self.pathwayProcessor.BKFsToFile(self.working_dir)
            reactome_bkf_files, reactome_source_names = self.reactomePathwayProcessor.BKFsToFile(self.working_dir)
        else:
            pathway_bkf_files, pathway_source_names = (list(), list())
            reactome_bkf_files, reactome_source_names = (list(), list())

        bkf_files = patient_bkf_files + pathway_bkf_files + reactome_bkf_files
        source_names = patient_source_names + pathway_source_names + reactome_source_names

        fuse(bkf_files,
             reliabilities=[1 for _ in range(len(bkf_files))],
             source_names=source_names)
        print('Patient data located at: {}'.format(dump_loc))
        return True


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--config_file', default=os.path.join(os.getcwd(), 'data_driver.config'), type=str)

    args = parser.parse_args()

    dataDriver = DataDriver(args.config_file)
    dataDriver.run()
