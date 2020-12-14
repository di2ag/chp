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

from chp.reasoner import Reasoner
from chp.query import Query
from chp.util import process_operator

#-- Point to Folder with Patient BKFs
#DEFAULT_BKF_FOLDER = '/home/public/data/ncats/BabelBKBs/smallProblem'
DEFAULT_BKF_FOLDER = '/home/public/data/ncats/BabelBKBs/collapsed100'

random.seed(123)

class CrossValidator:
    def __init__(self, test_ratio=None, fold=None, bkf_folder=None, compare_with=['random_forest'], result_folder=None):
        self.test_ratio = test_ratio
        self.bkf_folder = bkf_folder
        self.fold = fold
        self.compare_with = compare_with
        self.target_strategy = 'explicit'
        self.interpolation = 'standard'

        if bkf_folder is None:
            self.bkf_folder = DEFAULT_BKF_FOLDER

        #- Try to load in Patient Data
        try:
            with open(os.path.join(self.bkf_folder, 'patient_data.pk'), 'rb') as f_:
                self.patient_data = pickle.load(f_)
        except:
            raise IOError("Could not find patient_data.pk in {}".format(self.bkf_folder))

        #-- Make a results folder
        if result_folder is None:
            self.result_folder = os.path.join(os.getcwd(), 'results')
        else:
            self.result_folder = result_folder
        try:
            os.makedirs(self.result_folder)
        except:
            pass

        #-- Instaniate a logging file.
        logging.basicConfig(filename=os.path.join(self.result_folder, 'cross_valid.log'), level='DEBUG')

        #-- Run BKB Setup
        self.setup()

    def setup(self):
        t1 = time.time()
        if self.test_ratio is None and self.fold is None:
            raise('Both the test ratio and cross fold can not be None.')
        elif self.test_ratio is not None and self.fold is not None:
            raise('Either test ratio and cross fold must be  None.')
        elif self.test_ratio is not None:
            num_test = round(len(self.patient_data) * self.test_ratio)
            test_patients = random.sample(list(self.patient_data.keys()), num_test)
            train_patients = list(set(self.patient_data.keys()) - set(test_patients))
            fused_bkb = self.make_bkb(train_patients)

            self.fused_bkbs = [fused_bkb]
            self.test_patients = [test_patients]
            self.train_patients = [train_patients]
            logging.info('Setup Ok.')
            logging.info('Elapsed Time: {} sec'.format(time.time() - t1))
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
            for i in range(self.fold):
                test_patients = fold_patients[i]
                train_patients = set(patient_hashes) - set(test_patients)
                fused_bkb = self.make_bkb(train_patients)

                self.train_patients.append(list(train_patients))
                self.test_patients.append(test_patients)
                self.fused_bkbs.append(fused_bkb)
            logging.info('Setup Ok.')
            logging.info('Elapsed Time: {} sec'.format(time.time() - t1))
            return

    def make_bkb(self, train_patients):
            t1 = time.time()
            #-- Make Fused BKB with training paitents
            #-- Collect all patient BKFs
            bkfs = list()
            hashes = list()
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

            fused_bkb = fuse(bkfs,
                             [1 for _ in range(len(bkfs))],
                             hashes)
            logging.info('Make BKB Ok.')
            logging.info('Elapsed Time: {} sec'.format(time.time() - t1))
            return fused_bkb

    def run(self, target):
        logging.info('Runing Cross Validation...')
        target_op = process_operator(target[1])
        results = []
        for i, (fused_bkb, train_patients, test_patients) in enumerate(zip(self.fused_bkbs, self.train_patients, self.test_patients)):
            res = dict()
            with open(os.path.join(self.result_folder, 'cross_valid_data_{}.pk'.format(i)), 'wb') as f_:
                pickle.dump((fused_bkb, train_patients, test_patients), f_)
            supported_compnames = fused_bkb.getAllComponentNames()
            for patient in tqdm.tqdm(test_patients, desc='Running Inferences', leave=False):
                t1 = time.time()
                logging.info('Running patient {}'.format(patient))
                #-- Set Up Reasoner
                reasoner = Reasoner(fused_bkb, patient_data=self.patient_data)
                query = self.make_query(patient, target, supported_compnames)
                logging.info('Number of gene evidence = {}'.format(len(query.evidence)))
                if query is not None:
                    #-- Save out the query.
                    with open(os.path.join(self.result_folder, 'query_{}_data.pk'.format(patient)), 'wb') as f_:
                        pickle.dump(query, f_)
                    #-- Run query.
                    query = reasoner.analyze_query(query, target_strategy=self.target_strategy, interpolation=self.interpolation)
                    logging.info('Elapsed Time: {} sec'.format(time.time() - t1))
                    res[patient] = {'query': copy.deepcopy(query),
                                    'updates': [],
                                    'truth': []}
                    for update, prob in query.result.updates.items():
                        comp_idx, state_idx = update
                        comp_name = query.bkb.getComponentName(comp_idx)
                        state_name = query.bkb.getComponentINodeName(comp_idx, state_idx)
                        if target_op(float(self.patient_data[patient][target[0]]), target[2]) == ast.literal_eval(state_name):
                            truth = 1
                        else:
                            truth = 0
                        res[patient]['updates'].append((comp_name, state_name, prob))
                        res[patient]['truth'].append((comp_name, state_name, truth))
                else:
                    logging.warning('Patient {} had no gene mutations supported by Fused BKB.'.format(patient))
            with open(os.path.join(self.result_folder, 'cross_valid_res_{}.pk'.format(i)), 'wb') as f_:
                pickle.dump(res, f_)
            results.append(res)
        return results

    def make_query(self, patient_hash, target, supported_compnames, evidence_type='mutation'):
        evidence = dict()
        if evidence_type == 'mutation':
            patientMutationEvidence = dict()
            for mut in self.patient_data[patient_hash]["Patient_Genes"]:
                compName = '_mut_'+ mut
                if compName in supported_compnames:
                    patientMutationEvidence[compName] = 'True'
            if len(patientMutationEvidence) == 0:
                return None
            #print(patientMutationEvidence)
            query = Query(evidence=patientMutationEvidence,
                          targets=[],
                          meta_targets=[target])
            return query


if __name__ == '__main__':
    target = ('Survival_Time', '<=', 943)

    parser = argparse.ArgumentParser()
    parser.add_argument('--f', default=None, type=str, help='BKF Folder Path')
    parser.add_argument('--r', default=None, type=str, help='Results folder path')
    parser.add_argument('--t', default=None, type=float, help='Test ratio')
    parser.add_argument('--c', default=None, type=int, help='K-fold Cross Validation')

    args = parser.parse_args()

    cross_validator = CrossValidator(test_ratio=args.t, fold=args.c, bkf_folder=args.f, result_folder=args.r)
    results = cross_validator.run(target)
    for patient, res_dict in results[0].items():
        print(patient)
        print(res_dict['updates'])
        print(res_dict['truth'])
