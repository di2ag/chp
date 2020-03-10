import os
import sys
import itertools
from operator import ge, le, eq
import copy
import pandas as pd
import tqdm
import csv
import pickle
import time

from pybkb import bayesianKnowledgeBase as BKB
from pybkb.core.common.bayesianKnowledgeBase import BKB_I_node, BKB_component, BKB_S_node
from pybkb.core.python_base.fusion import fuse

sys.path.append('/home/ghyde/bkb-pathway-provider/core')
#sys.path.append('/home/cyakaboski/src/python/projects/bkb-pathway-provider/core')

from reasoner import Reasoner
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
    def __init__(self, bkb, test_patient_hashes, patient_data_file, patient_dict):
        self.bkb = bkb
        self.queryCopy = None
        self.test_patient_hashes = test_patient_hashes
        self.patient_data_file = patient_data_file
        self.patient_dict = patient_dict
        self.drug = None

        self.target_strategy = 'topological'
        self.interpolation = 'independence'

    def run_demo_suite(self, target, patientIDs, patientsMutationEvidence):
        #-- Set Up Reasoner
        self.reasoner = Reasoner(self.bkb, None)
        self.reasoner.set_src_metadata(self.patient_data_file)
        self.reasoner.cpp_reasoning = False
        holdoutResults = list()

        for idx, patID in tqdm.tqdm(enumerate(patientIDs)):
            probs = list()
            for mutationDrugEvidence in patientsMutationEvidence[idx]:
                print(patientsMutationEvidence[idx])
                prob = self.run_demo_only(target, patID, mutationDrugEvidence)
                probs.append(prob)
            probOne = 1.0
            probTwo = 1.0
            for prob in probs:
                sumProbs = prob[0][2] + prob[1][2]
                probOne *= (prob[0][2]/sumProbs)
                probTwo *= (prob[1][2]/sumProbs)
                sumProbs = probOne + probTwo
                probOne /= sumProbs
                probTwo /= sumProbs
            prob = list()
            prob.append((probs[0][0], probs[0][1], probOne))
            prob.append((probs[1][0], probs[1][1], probTwo))
            holdoutResults.append((prob[0],prob[1]))

        self.getResults(target, holdoutResults, patientIDs)

    def run_demo_only(self, target, patID, evidence):
        #-- Make query and analyze
        query = Query(evidence=evidence[1],
                      targets=[],
                      meta_evidence=evidence[0],
                      meta_targets=[target],
                      type='updating')
        probs = list()
        if self.drug != evidence[0][0][2]:
            query = self.reasoner.analyze_query(copy.deepcopy(query), save_dir=None,
                                               target_strategy=self.target_strategy, interpolation=self.interpolation)
            self.processed_bkb = copy.deepcopy(query.bkb)
            self.drug == evidence[0][0][2]
            query.getReport()
            for comp_name, state_dict in query.independ_result.items():
                for state_name, prob in state_dict.items():
                    probs.append((comp_name, state_name, prob))
        else:
            query = self.reasoner.analyze_query(copy.deepcopy(query), preprocessed_bkb=self.processed_bkb,
                                               target_strategy=self.target_strategy, interpolation=self.interpolation)
            query.getReport()
            for comp_name, state_dict in query.independ_result.items():
                for state_name, prob in state_dict.items():
                    probs.append((comp_name, state_name, prob))
        print(probs)
        return (probs[0],probs[1])

    def getResults(self, target, holdoutResults, patientIDs):
        assert len(holdoutResults) == len(self.test_patient_hashes), "results mismatch with test patients"
        print("P("+target[0] + " " + target[1] + " " + str(target[2]) + " | geneVariants)")
        numCorrect = 0
        numWrong = 0
        for idx, tph in enumerate(self.test_patient_hashes):
            targetVal = int(self.patient_dict[int(tph)][target[0]])
            patID = self.patient_dict[int(tph)]['Patient_ID']
            assert patID == patientIDs[idx], "holdout patient mismatch"
            op = _process_operator(target[1])

            indexTrue = None
            indexFalse = None
            if holdoutResults[idx][0][1] == 'True':
                indexTrue = 0
                indexFalse = 1
            else:
                indexTrue = 1
                indexFalse = 0

            probOne = holdoutResults[idx][indexFalse][2] #false
            probTwo = holdoutResults[idx][indexTrue][2] #true
            sum = probOne + probTwo
            probOne /= sum
            probTwo /= sum
            if op(targetVal, target[2]):
                print(patID)
                print(holdoutResults[idx][indexFalse],holdoutResults[idx][indexTrue])
                print(probOne, probTwo)
                if holdoutResults[idx][indexTrue][2] > holdoutResults[idx][indexFalse][2]:
                    print("\ttrue -",targetVal, "CORRECT")
                    numCorrect += 1
                else:
                    print("\ttrue -",targetVal, "WRONG")
                    numWrong += 1
            else:
                print(patID)
                print(holdoutResults[idx][indexFalse],holdoutResults[idx][indexTrue])
                print(probOne, probTwo)
                if holdoutResults[idx][indexFalse][2] > holdoutResults[idx][indexTrue][2]:
                    print("\tfalse -",targetVal, "CORRECT")
                    numCorrect += 1
                else:
                    print("\tfalse -",targetVal, "WRONG")
                    numWrong += 1
        print("Number correct:", numCorrect)
        print("Number wrong:", numWrong)

if __name__ == '__main__':
    fused_bkb = BKB()

    fused_bkb.load('/home/public/data/ncats/660Pats6Holdout/fusion.bkb')
    patient_data_file = '/home/public/data/ncats/660Pats6Holdout/patient_data.pk'
    withheld_patients_file = '/home/public/data/ncats/660Pats6Holdout/withheldPatients.csv'

    #fused_bkb.load('/home/public/data/ncats/660Pats5PercentHoldout/fusion.bkb')
    #patient_data_file = '/home/public/data/ncats/660Pats5PercentHoldout/patient_data.pk'
    #withheld_patients_file = '/home/public/data/ncats/660Pats5PercentHoldout/withheldPatients.csv'

    compNames = fused_bkb.getAllComponentNames()
    f = open(withheld_patients_file, 'r')
    withheldPatientHashes = f.read().split(',')

    patientDict = pickle.load(open(patient_data_file, 'rb'))
    patientsMutationDrugEvidence = []
    patientIDs = []
    for withheldPatientHash in withheldPatientHashes:
        withheldPatientDict = patientDict[int(withheldPatientHash)]
        patientIDs.append(withheldPatientDict["Patient_ID"])
        pgv = withheldPatientDict["Patient_Genes"]
        patientMutationEvidence = dict()
        patientMutationDrugEvidence = []
        for mut in pgv:
            compName = 'mut_'+mut
            if compName in compNames:
                patientMutationEvidence[compName] = 'True'
        #patientsMutationEvidence.append(patientMutationEvidence)
        drugs = withheldPatientDict["Drug_Name(s)"]
        for drug in drugs:
            drugEvidence = [('Drug_Name(s)', '==', drug)]
            patientMutationDrugEvidence.append((drugEvidence, patientMutationEvidence))
        patientsMutationDrugEvidence.append(patientMutationDrugEvidence)
    target = ('Survival_Time', '<=', 943)

    cross_validator = CrossValidator(fused_bkb,withheldPatientHashes, patient_data_file, patientDict)
    df = cross_validator.run_demo_suite(target,
                                        patientIDs,
                                        patientsMutationDrugEvidence)
