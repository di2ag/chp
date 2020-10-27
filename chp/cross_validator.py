'''
Source code developed by DI2AG.
Thayer School of Engineering at Dartmouth College
Authors:    Dr. Eugene Santos, Jr
            Mr. Chase Yakaboski,
            Mr. Gregory Hyde,
            Dr. Keum Joo Kim
'''


import os
import random
import sys
import itertools
from operator import ge, le, eq
import copy
import pandas as pd
import tqdm
import csv
import pickle
import time
import argparse
import ast
import logging

from pybkb.common.bayesianKnowledgeBase import bayesianKnowledgeBase as BKB
from pybkb.common.bayesianKnowledgeBase import BKB_I_node, BKB_component, BKB_S_node
from pybkb.python_base.fusion import fuse
from chp_data.bkb_handler import BkbDataHandler

from chp.reasoner import Reasoner
from chp.query import Query
from chp.util import process_operator
from chp.patientBKFProcessor import PatientProcessor

#-- Point to Folder with Patient BKFs
#DEFAULT_BKF_FOLDER = '/home/public/data/ncats/BabelBKBs/smallProblem'
DEFAULT_BKF_FOLDER = '/home/public/data/ncats/BabelBKBs/collapsedAll'

random.seed(123)

#-- Setup logging
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

class CrossValidator:
    def __init__(self,
                 test_ratio=None,
                 fold=None,
                 bkf_folder=None,
                 compare_with=['random_forest'],
                 result_folder=None,
                 hosts_filename=None,
                 num_procosses_per_host=0,
                 venv=None):
        self.test_ratio = test_ratio
        self.bkf_folder = bkf_folder
        self.fold = fold
        self.compare_with = compare_with
        self.target_strategy = 'explicit'
        self.interpolation = 'standard'
        self.interpolation_model = 'bigram'
        self.interpolation_selection = 'frequency_based'
        self.frequency_threshold = 5
        self.hosts_filename = hosts_filename
        self.num_processes_per_host = num_procosses_per_host
        self.venv = venv

        if bkf_folder is None:
            #bkb_handler = BkbDataHandler(bkb_version='special', dataset_version='babel-small-problem')
            bkb_handler = BkbDataHandler()
            with open(bkb_handler.patient_data_pk_path, 'rb') as f_:
                self.patient_data = pickle.load(f_)
            self.patient_processor = PatientProcessor()
            self.patient_processor.loadFromPatientData(bkb_handler.patient_data_pk_path, patient_fraction=1)
        else:
            #- Try to load in Patient Data
            try:
                with open(os.path.join(self.bkf_folder, 'patient_data.pk'), 'rb') as f_:
                    self.patient_data = pickle.load(f_)
                #self.patient_processor = PatientProcessor()
                #self.patient_processor.loadFromPatientData(os.path.join(self.bkf_folder, 'patient_data.pk'))
            except:
                raise IOError("Could not find patient_data.pk in either bkf folder or chp_data")

        #-- Make a results folder
        if result_folder is None:
            self.result_folder = os.path.join(os.getcwd(), 'results')
        else:
            self.result_folder = result_folder
        try:
            os.makedirs(self.result_folder)
        except:
            pass

        #-- Instaniate a logger file.
        #logger.basicConfig(filename=os.path.join(self.result_folder, 'cross_valid.log'), level='DEBUG')

        #-- Run BKB Setup
        self.setup()

    def setup(self):
        t1 = time.time()
        if self.test_ratio is None and self.fold is None:
            raise('Both the test ratio and cross fold can not be None.')
        elif self.test_ratio is not None and self.fold is not None:
            raise('Either test ratio and cross fold must be  None.')
        elif self.test_ratio is not None:
            if self.bkf_folder is None:
                num_tests = round(len(self.patient_processor.patients) * self.test_ratio)
                test_patients = [pat.patientHash for pat in random.sample(self.patient_processor.patients, num_tests)]
                train_patients = list(set([pat.patientHash for pat in self.patient_processor.patients]) - set(test_patients))
            else:
                num_test = round(len(self.patient_data) * self.test_ratio)
                test_patients = random.sample(list(self.patient_data.keys()), num_test)
                train_patients = list(set(self.patient_data.keys()) - set(test_patients))
            #-- Delay make BKB if entropy-based.
            if self.interpolation_selection != 'entropy_based':
                fused_bkb, bad_genes = self.make_bkb(train_patients)
                self.fused_bkbs = [fused_bkb]
                self.bad_genes = [bad_genes]
            self.test_patients = [test_patients]
            self.train_patients = [train_patients]
            logger.info('Setup Ok.')
            logger.info('Elapsed Time: {} sec'.format(time.time() - t1))
            return
        else:
            num_fold = len(self.patient_data) // self.fold
            fold_patients = list()
            patient_hashes = list(self.patient_data.keys())
            for i in range(0, step=num_fold):
                try:
                    fold_patients.append(patient_hashes[i:i+num_fold])
                except:
                    fold_patients.append(patient_hashes[i:])
            self.train_patients = list()
            self.test_patients = list()
            self.fused_bkbs = list()
            self.bad_genes = list()
            for i in range(self.fold):
                test_patients = fold_patients[i]
                train_patients = set(patient_hashes) - set(test_patients)
                #-- Delay make BKB if entropy-based.
                if self.interpolation_selection != 'entropy_based':
                    fused_bkb, bad_genes = self.make_bkb(train_patients)
                    self.fused_bkbs.append(fused_bkb)
                    self.bad_genes.append(bad_genes)
                self.train_patients.append(list(train_patients))
                self.test_patients.append(test_patients)
            logger.info('Setup Ok.')
            logger.info('Elapsed Time: {} sec'.format(time.time() - t1))
            return

    def make_bkb(self, train_patients, target=None):
            t1 = time.time()
            #-- Make Fused BKB with training paitents
            #-- Collect all patient BKFs
            bkfs = list()
            hashes = list()
            bad_genes = None
            if self.bkf_folder is not None:
                for patient_hash in tqdm.tqdm(train_patients, leave=False, desc='Preparing BKFs for Fusion'):
                    bkf_file = os.path.join(self.bkf_folder, '{}.bkf'.format(patient_hash))
                    _bkf = BKB()
                    try:
                        _bkf.load(bkf_file)
                        bkfs.append(_bkf)
                        hashes.append(str(patient_hash))
                    except:
                        continue
                #-- Try to add PatientX
                try:
                    bkf_file = os.path.join(self.bkf_folder, 'PatientX.bkf')
                    _bkf = BKB()
                    _bkf.load(bkf_file)
                    bkfs.append(_bkf)
                    hashes.append('PatientX')
                except:
                    pass
            else:
                #-- Process BKFs
                bad_genes = self.patient_processor.processPatientBKF_v2(interpolation_model=self.interpolation_model,
                                                                        interpolation_selection=self.interpolation_selection,
                                                                        frequency_threshold=self.frequency_threshold,
                                                                        entropy_criteria='min',
                                                                        phenotypic_evidence=target,
                                                                        patient_hashes_to_process=train_patients,
                                                                        gene_subset_top_k=None)
                #-- Extract patient bkfs based on hashes
                for patient, bkf in zip(self.patient_processor.requested_patients, self.patient_processor.bkfs):
                    bkfs.append(bkf)
                    hashes.append(str(patient.patientHash))
                #-- Add interpolator
                bkfs.append(self.patient_processor.interpolator)
                hashes.append('interpolator')
            fused_bkb = fuse(bkfs,
                             [1 for _ in range(len(bkfs))],
                             hashes)
            fused_bkb.save('fused_bkb_sample.bkb', use_pickle=True)
            logger.info('Make BKB Ok.')
            logger.info('Elapsed Time: {} sec'.format(time.time() - t1))
            return fused_bkb, bad_genes

    def run(self, target):
        if self.interpolation_selection == 'entropy_based':
            self.fused_bkbs = []
            self.bad_genes = []
            for train_patients_ in self.train_patients:
                fused_bkb, bad_genes = self.make_bkb(train_patients_, target)
                self.fused_bkbs.append(fused_bkb)
                self.bad_genes.append(bad_genes)
        logger.info('Runing Cross Validation...')
        target_op = process_operator(target[0][1])
        results = []
        for i, (fused_bkb, bad_genes, train_patients, test_patients) in enumerate(zip(self.fused_bkbs, self.bad_genes, self.train_patients, self.test_patients)):
            res = dict()
            with open(os.path.join(self.result_folder, 'cross_valid_data_{}.pk'.format(i)), 'wb') as f_:
                pickle.dump((fused_bkb, train_patients, test_patients), f_)
            supported_compnames = fused_bkb.getAllComponentNames()
            for patient in tqdm.tqdm(test_patients, desc='Running Inferences', leave=False):
                t1 = time.time()
                logger.info('Running patient {}'.format(patient))
                #-- Set Up Reasoner
                reasoner = Reasoner(fused_bkb=fused_bkb, patient_data=self.patient_data, hosts_filename=self.hosts_filename, num_processes_per_host=self.num_processes_per_host, venv=self.venv)
                query = self.make_query(patient, target, supported_compnames, bad_genes)
                if query is None:
                    logger.warning('Query evidence is None, likely cause no genes were supported.')
                else:
                    logger.info('Number of gene evidence = {}'.format(len(query.evidence)))
                if query is not None:
                    #-- Run query.
                    query = reasoner.analyze_query(query, target_strategy=self.target_strategy, interpolation=self.interpolation)
                    #for thing in query.result.result:
                    #    print(thing)
                    #-- Save out the query.
                    #query.save('results')
                    with open(os.path.join(self.result_folder, 'query_{}_data.pk'.format(patient)), 'wb') as f_:
                        pickle.dump(query, f_)
                    logger.info('Elapsed Time: {} sec'.format(time.time() - t1))
                    res[patient] = {'query': copy.deepcopy(query),
                                    'updates': [],
                                    'truth': []}
                    for update, prob in query.result.updates.items():
                        comp_idx, state_idx = update
                        comp_name = query.bkb.getComponentName(comp_idx)
                        state_name = query.bkb.getComponentINodeName(comp_idx, state_idx)
                        if target_op(float(self.patient_data[patient][target[0][0]]), target[0][2]) == ast.literal_eval(state_name):
                            truth = 1
                        else:
                            truth = 0
                        res[patient]['updates'].append((comp_name, state_name, prob))
                        res[patient]['truth'].append((comp_name, state_name, truth))
                else:
                    logger.warning('Patient {} had no gene mutations supported by Fused BKB.'.format(patient))
            with open(os.path.join(self.result_folder, 'cross_valid_res_{}.pk'.format(i)), 'wb') as f_:
                pickle.dump(res, f_)
            results.append(res)
        return results

    def make_query(self, patient_hash, target, supported_compnames, bad_genes, evidence_type='mutation'):
        evidence = dict()
        if evidence_type == 'mutation':
            patientMutationEvidence = dict()
            for mut in self.patient_data[patient_hash]["Patient_Genes"]:
                if mut not in bad_genes:
                    compName = '_mut_'+ mut
                    if compName in supported_compnames:
                        patientMutationEvidence[compName] = 'True'
            if len(patientMutationEvidence) == 0:
                return None
            logger.debug('{}'.format(patientMutationEvidence))
            query = Query(evidence=patientMutationEvidence,
                          targets=[],
                          meta_targets=target)
            return query


if __name__ == '__main__':
    target = [('Survival_Time', '<=', 943)]

    parser = argparse.ArgumentParser()
    parser.add_argument('--f', default=None, type=str, help='BKF Folder Path')
    parser.add_argument('--r', default=None, type=str, help='Results folder path')
    parser.add_argument('--t', default=.2, type=float, help='Test ratio')
    parser.add_argument('--c', default=None, type=int, help='K-fold Cross Validation')
    parser.add_argument('--hosts', default=None, type=str, help='Hosts file for distributed reasoning.')
    parser.add_argument('--numproc', default=0, type=int, help='Number of processes to be run on each worker during distributed reasoning')
    parser.add_argument('--venv', default=None, type=str, help='Virtual environment path to use for distributed reasoning.')
    args = parser.parse_args()

    cross_validator = CrossValidator(test_ratio=args.t, fold=args.c, bkf_folder=args.f, result_folder=args.r, hosts_filename=args.hosts, num_procosses_per_host=args.numproc, venv=args.venv)
    #cross_validator = CrossValidator(test_ratio=args.t, fold=args.c, bkf_folder=DEFAULT_BKF_FOLDER, result_folder=args.r)
    results = cross_validator.run(target)
    for patient, res_dict in results[0].items():
        print(patient)
        print(res_dict['updates'])
        print(res_dict['truth'])
