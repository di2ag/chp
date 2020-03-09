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

#sys.path.append('/home/ghyde/bkb-pathway-provider/core')
sys.path.append('/home/cyakaboski/src/python/projects/bkb-pathway-provider/core')

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
    def __init__(self, bkb, test_patient_hashes, patient_data_file):
        self.bkb = bkb
        self.queryCopy = None
        self.test_patient_hashes = test_patient_hashes
        self.patient_data_file = patient_data_file
        self.first = True

        #-- You can change target_strategy to either: chain or topological
        #-- You can change interpolation to either: standard or independence
        self.target_strategy = 'topological'
        self.interpolation = 'independence'

    def run_demo_suite(self, target, patientIDs, patientsMutationEvidence):
        #-- Set Up Reasoner
        self.reasoner = Reasoner(self.bkb, None)
        self.reasoner.set_src_metadata(self.patient_data_file)
        self.reasoner.cpp_reasoning = False

        for idx, patID in tqdm.tqdm(enumerate(patientIDs)):
            queryResults = list()
            probs = list()
            summaries = list()
            evs = list()
            probOneState = ""
            probOne = 1
            probTwoState = ""
            probTwo = 1
            badGenes = list()

            prob, summary = self.run_demo_only(target, patID, patientsMutationEvidence[idx])
            evs.append(patientsMutationEvidence[idx])
            queryResults.append(patID)
            probs.append(prob)
            summaries.append(summary)
            for i in range(0, len(queryResults)):
                probOneState = probs[i][0][1]
                probTwoState = probs[i][1][1]
                one = float(probs[i][0][2])
                two = float(probs[i][1][2])
                print(target, evs[i])
                print(probOneState, one, probTwoState, two, "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$")
                print(summaries[i])
                if one == -1 and two == -1:
                    print("FUCK")
                elif one == -1 and two != -1:
                    two = 1.0
                    one = 0.0
                elif one != -1 and two == -1:
                    one = 1.0
                    two = 0
                else:
                    sum = one + two
                    one /= sum
                    two /= sum
                if one != 0 and two != 0:
                    probOne *= one
                    probTwo *= two
                    sum = probOne + probTwo
                    probOne /= sum
                    probTwo /= sum
                else:
                    badGenes.append(evs[i])
            print(patID,"@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@")
            print(probOneState,probTwoState)
            print(probOne, probTwo)
            for badGene in badGenes:
                print(badGene, end='', sep=', ')
                #print(summaries[i])

    def run_demo_only(self, target, patID, evidence):
        #print(evidence)
        #print([target])
        #-- Make query and analyze
        query = Query(evidence=evidence,
                      targets=[],
                      meta_evidence=None,
                      meta_targets=[target],
                      type='updating')
        probs = list()
        summary = None
        if self.first:
            query = self.reasoner.analyze_query(copy.deepcopy(query), save_dir=None,
                                               target_strategy=self.target_strategy, interpolation=self.interpolation)
            self.processed_bkb = copy.deepcopy(query.bkb)
            self.first = False
            query.getReport()
            #summary = query.result.summary()
            #query.result.summary()
            for comp_name, state_dict in query.independ_result.items():
                for state_name, prob in state_dict.items():
                    probs.append((comp_name, state_name, prob))
        else:
            query = self.reasoner.analyze_query(copy.deepcopy(query), preprocessed_bkb=self.processed_bkb,
                                               target_strategy=self.target_strategy, interpolation=self.interpolation)
            query.getReport()
            #summary = query.result.summary()
            #query.result.summary()
            for comp_name, state_dict in query.independ_result.items():
                for state_name, prob in state_dict.items():
                    probs.append((comp_name, state_name, prob))

        return (probs[0],probs[1]),summary


if __name__ == '__main__':
    fused_bkb = BKB()

    fused_bkb.load('/home/public/data/ncats/660Pats6Holdout/fusion.bkb')
    patient_data_file = '/home/public/data/ncats/660Pats6Holdout/patient_data.pk'
    withheld_patients_file = '/home/public/data/ncats/660Pats6Holdout/withheldPatients.csv'

    #fused_bkb.load('/home/public/data/ncats/660Pats20PercentHoldout/fusion.bkb')
    #patient_data_file = '/home/public/data/ncats/660Pats20PercentHoldout/patient_data.pk'
    #withheld_patients_file = '/home/public/data/ncats/660Pats20PercentHoldout/withheldPatients.csv'

    compNames = fused_bkb.getAllComponentNames()

    f = open(withheld_patients_file, 'r')
    withheldPatientHashes = f.read().split(',')

    patientDict = pickle.load(open(patient_data_file, 'rb'))
    patientsMutationEvidence = []
    '''
    for comp_name, inode_name in fused_bkb.getINodeNames():
        if 'mut-var_' == comp_name[:8]:
            gene = '_'.join(comp_name.split('_')[1:])
            variant = inode_name
            gene_variant = '{}-{}'.format(gene, variant)
            count = 0
            for pat_hash, data_dict in patientDict.items():
                #try:
                if gene_variant in data_dict['Patient_Gene_Variants']:
                    count += 1
                #except:
                #    continue
            if count == 0:
                print("Couldn't find: {} = {}".format(comp_name, inode_name))
                print('Rather: {} = {}'.format(gene, variant))
    input('Done')
    '''
    patientIDs = []
    count = 0
    for withheldPatientHash in withheldPatientHashes:
        withheldPatientDict = patientDict[int(withheldPatientHash)]
        if count < 102:
            patientIDs.append(withheldPatientDict["Patient_ID"])
            pgv = withheldPatientDict["Patient_Genes"]
            patientMutationEvidence = dict()
            for mut in pgv:
                compName = 'mut_'+mut
                if compName in compNames:
                    patientMutationEvidence[compName] = 'True'
                    #patientDynamicEvidence.append(dict) #          [('Patient_Gene_Variants', '==', tuple([mut]))])
            patientsMutationEvidence.append(patientMutationEvidence)
        count += 1
    target = ('Survival_Time', '<=', 943)
    #print(patientsMutationEvidence)

    cross_validator = CrossValidator(fused_bkb,withheldPatientHashes, patient_data_file)
    df = cross_validator.run_demo_suite(target,
                                        patientIDs,
                                        patientsMutationEvidence)
    #print(df)
