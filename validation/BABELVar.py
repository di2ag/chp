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
        self.processed_bkb = None
        self.queryCopy = None
        self.test_patient_hashes = test_patient_hashes
        self.patient_data_file = patient_data_file
        self.patient_dict = patient_dict
        self.first = True

        self.target_strategy = 'chain'
        self.interpolation = 'standard'

    def run_demo_suite(self, target, patientIDs, patientsMutationVariantEvidence):
        #-- Set Up Reasoner
        self.reasoner = Reasoner(self.bkb, None)
        self.reasoner.set_src_metadata(self.patient_data_file)
        self.reasoner.cpp_reasoning = False
        holdoutResults = list()

        for idx, patID in tqdm.tqdm(enumerate(patientIDs)):
            prob = self.run_demo_only(target, patID, patientsMutationVariantEvidence[idx])
            holdoutResults.append((prob[0],prob[1]))

        self.getResults(target, holdoutResults, patientIDs)

    def run_demo_only(self, target, patID, evidence):
        print(evidence)
        print([target])
        #-- Make query and analyze
        query = Query(evidence=evidence,
                      targets=[],
                      meta_evidence=None,
                      meta_targets=[target],
                      type='updating')
        probs = list()
        summary = None
        if self.first:
            query = self.reasoner.analyze_query(copy.deepcopy(query),
                                                target_strategy=self.target_strategy, interpolation=self.interpolation)
            self.processed_bkb = copy.deepcopy(query.bkb)
            self.first = False
            query.result.summary()
            query.getReport()
            for update, prob in query.result.updates.items():
                comp_idx, state_idx = update
                comp_name = query.bkb.getComponentName(comp_idx)
                state_name = query.bkb.getComponentINodeName(comp_idx, state_idx)
                probs.append((comp_name, state_name, prob))

        else:
            query = self.reasoner.analyze_query(copy.deepcopy(query), preprocessed_bkb=self.processed_bkb,
                                                target_strategy=self.target_strategy, interpolation=self.interpolation)
            query.result.summary()
            query.getReport()
            for update, prob in query.result.updates.items():
                comp_idx, state_idx = update
                comp_name = query.bkb.getComponentName(comp_idx)
                state_name = query.bkb.getComponentINodeName(comp_idx, state_idx)
                probs.append((comp_name, state_name, prob))

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

    fused_bkb.load('/home/public/data/ncats/BabelBKBs/collapsedAll/fusion.bkb')
    patient_data_file = '/home/public/data/ncats/BabelBKBs/collapsedAll/patient_data.pk'
    withheld_patients_file = '/home/public/data/ncats/BabelBKBs/collapsedAll/withheldPatients.csv'

    #fused_bkb.load('/home/public/data/ncats/660Pats5PercentHoldout/fusion.bkb')
    #patient_data_file = '/home/public/data/ncats/660Pats5PercentHoldout/patient_data.pk'
    #withheld_patients_file = '/home/public/data/ncats/660Pats5PercentHoldout/withheldPatients.csv'

    compNames = fused_bkb.getAllComponentNames()
    f = open(withheld_patients_file, 'r')
    withheldPatientHashes = f.read().split(',')

    patientDict = pickle.load(open(patient_data_file, 'rb'))
    patientsMutationVariantEvidence = []
    patientIDs = []
    for withheldPatientHash in withheldPatientHashes:
        withheldPatientDict = patientDict[int(withheldPatientHash)]
        patientIDs.append(withheldPatientDict["Patient_ID"])
        genes = withheldPatientDict["Patient_Genes"]
        variants = withheldPatientDict["Patient_Variants"]
        patientMutationVariantEvidence = dict()
        for idx, mut in enumerate(genes[0:5]):
            compName = 'mut-var_'+mut
            if compName in compNames:
                compIDX = fused_bkb.getComponentIndex(compName)
                compStatesIDXs = fused_bkb.getAllComponentINodeIndices(compIDX)
                stateNames = [fused_bkb.getComponentINodeName(compIDX,csIdx) for csIdx in compStatesIDXs]
                if variants[idx] in stateNames:
                    patientMutationVariantEvidence[compName] = variants[idx]
        patientsMutationVariantEvidence.append(patientMutationVariantEvidence)
    target = ('Survival_Time', '<=', 943)

    cross_validator = CrossValidator(fused_bkb,withheldPatientHashes, patient_data_file, patientDict)
    cross_validator.run_demo_suite(target,
                                   patientIDs,
                                   patientsMutationVariantEvidence)


