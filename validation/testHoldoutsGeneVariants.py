import os
import sys
import itertools
from operator import ge, le, eq
import copy
import pandas as pd
import tqdm
import csv
import pickle

from pybkb import bayesianKnowledgeBase as BKB
from pybkb.core.common.bayesianKnowledgeBase import BKB_I_node, BKB_component, BKB_S_node
from pybkb.core.python_base.fusion import fuse

sys.path.append('/home/ghyde/bkb-pathway-provider/core')

from reasoner import Reasoner, _constructSNodesByHead
from query import Query

def read_withheldPatientFile(file_name):
    hashes = list()
    with open(file_name, 'r') as csv_file:
        reader = csv.reader(csv_file)
        for row in reader:
            for patient_hash in row:
                hashes.append(int(patient_hash))
    return hashes

def _process_operator(op):
    if op == '>=':
        return ge
    elif op == '<=':
        return le
    elif op == '==':
        return eq
    else:
        raise ValueError('Unknown Operator')

class CrossValidator:
    def __init__(self, bkb, test_patient_hashes, patient_data_file):
        self.bkb = bkb
        self.queryCopy = None
        self.test_patient_hashes = test_patient_hashes
        self.patient_data_file = patient_data_file
        self.first = True

    def run_demo_suite(self, target, patientIDs, patientsDynamicEvidence):
        #queryResults = list()
        #probs = list()
        #summaries = list()
        for idx, patID in tqdm.tqdm(enumerate(patientIDs)):
            queryResults = list()
            probs = list()
            summaries = list()
            for patientDynamicEvidence in tqdm.tqdm(patientsDynamicEvidence[idx]):
                prob, summary = self.run_demo_only(target, patID, patientDynamicEvidence)
                queryResults.append(patID)
                probs.append(prob)
                summaries.append(summary)
            for i in range(0, len(queryResults)):
                print(probs[i])
                print(summaries[i])

    def run_demo_only(self, target, patID, evidence):
        #-- Set Up Reasoner
        self.reasoner = Reasoner(self.bkb, None)
        self.reasoner.set_src_metadata(self.patient_data_file)
        self.reasoner.cpp_reasoning = False

        print(evidence)
        print([target])
        #-- Make query and analyze
        query = Query(evidence=evidence,
                      targets=[],
                      meta_evidence=[],
                      meta_targets=[target],
                      type='updating')
        probs = list()
        summary = None
        if self.first:
            query, queryCopy= self.reasoner.analyze_query(copy.deepcopy(query), save_dir=None)
            self.queryCopy = queryCopy
            self.first = False
            summary = query.result.summary()
            for update, prob in query.result.updates.items():
                comp_idx, state_idx = update
                comp_name = query.bkb.getComponentName(comp_idx)
                state_name = query.bkb.getComponentINodeName(comp_idx, state_idx)
                probs.append((comp_name, state_name, prob))
        else:
            query = self.queryCopy
            query.meta_targets = [target]
            query.evidence = evidence
            query = self.reasoner.solve_query(copy.deepcopy(query))
            summary = query.result.summary()
            for update, prob in query.result.updates.items():
                comp_idx, state_idx = update
                comp_name = query.bkb.getComponentName(comp_idx)
                state_name = query.bkb.getComponentINodeName(comp_idx, state_idx)
                probs.append((comp_name, state_name, prob))

        return probs, summary


if __name__ == '__main__':
    fused_bkb = BKB()

    #fused_bkb.load('/home/public/data/ncats/663Pats6Holdouts/fusion.bkb')
    #patient_data_file = '/home/public/data/ncats/663Pats6Holdouts/patient_data.pk'
    #withheld_patients_file = '/home/public/data/ncats/663Pats6Holdouts/withheldPatients.csv'

    fused_bkb.load('/home/public/data/ncats/90PERCENTValidation50Patients10Validation/set1/fusion.bkb')
    patient_data_file = '/home/public/data/ncats/90PERCENTValidation50Patients10Validation/set1/patient_data.pk'
    withheld_patients_file = '/home/public/data/ncats/90PERCENTValidation50Patients10Validation/set1/withheldPatients.csv'


    compNames = fused_bkb.getAllComponentNames()
    f = open(withheld_patients_file, 'r')
    withheldPatientHashes = f.read().split(',')

    patientDict = pickle.load(open(patient_data_file, 'rb'))
    patientsDynamicEvidence = []
    patientIDs = []
    count = 0
    for withheldPatientHash in withheldPatientHashes:
        withheldPatientDict = patientDict[int(withheldPatientHash)]
        if count < 1:
            patientIDs.append(withheldPatientDict["Patient_ID"])
            pgv = withheldPatientDict["Patient_Genes"]
            patientDynamicEvidence = []
            for mut in pgv[0:5]:
                compName = 'mut_'+mut+'='
                if compName in compNames:
                    dict = {compName:'True'}
                    patientDynamicEvidence.append(dict) #          [('Patient_Gene_Variants', '==', tuple([mut]))])
            patientsDynamicEvidence.append(patientDynamicEvidence)
        count += 1
    target = ('Survival_Time', '<=', 943)

    cross_validator = CrossValidator(fused_bkb,withheldPatientHashes, patient_data_file)
    df = cross_validator.run_demo_suite(target,
                                        patientIDs,
                                        patientsDynamicEvidence)
    print(df)
