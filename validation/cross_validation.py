import os
import sys
import itertools
from operator import ge, le, eq
import copy
import pandas as pd
import tqdm
import csv

from pybkb import bayesianKnowledgeBase as BKB
from pybkb.core.common.bayesianKnowledgeBase import BKB_I_node, BKB_component, BKB_S_node
from pybkb.core.python_base.fusion import fuse

sys.path.append('/home/cyakaboski/src/python/projects/bkb-pathway-provider/core')

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
        self.test_patient_hashes = test_patient_hashes
        self.patient_data_file = patient_data_file

    def run_demo_suite(self, demo_target, demo_evidence_static=list(), demo_evidence_dynamic=list(), demo_evidence_values=list()):
        results = list()
        for val in tqdm.tqdm(demo_evidence_values):
            res = self.run_demo_only(demo_target, demo_evidence_static,demo_evidence_dynamic, val)
            results.append(res)

        res_collection = dict()
        for res in results:
            for key, item in res.items():
                if key not in res_collection:
                    res_collection[key] = [item]
                else:
                    res_collection[key].append(item)

        return self.results_to_dataframe(res_collection)

    def run_demo_only(self, demo_target, demo_evidence_static=list(), demo_evidence_dynamic=list(), demo_evidence_val=0):
        #-- Set Up Reasoner
        reasoner = Reasoner(self.bkb, None)
        reasoner.set_src_metadata(self.patient_data_file)
        reasoner.cpp_reasoning = False

        #-- Initialize Result dict
        result = {'X': demo_evidence_val,
                   'Compute_time': -1}

        #-- Setup demographic evidence tuple
        if len(demo_evidence_dynamic) > 0:
            demo_evidence_dynamic_run = [tuple(list(demo_evidence_dynamic) + [demo_evidence_val])]
        demo_evidence_run = demo_evidence_static + demo_evidence_dynamic_run

        #-- Make query and analyze
        query = Query(meta_evidence=demo_evidence_run,
                      meta_targets=[demo_target],
                      type='updating')
        query = reasoner.analyze_query(copy.deepcopy(query), save_dir=None)
        result['Compute_time'] = query.compute_time

        #-- Collect Updates and put into results
        for update, prob in query.result.updates.items():
            comp_idx, state_idx = update
            comp_name = query.bkb.getComponentName(comp_idx)
            state_name = query.bkb.getComponentINodeName(comp_idx, state_idx)

            result[('Updates', comp_name, state_name)] = prob
        return result

    def run_patients(self, demo_target, gene_evidence=False, demo_evidence=None, demo_value=0):
        reasoner = Reasoner(self.bkb, None)
        reasoner.set_src_metadata(self.patient_data_file)
        reasoner.cpp_reasoning = False
        #-- Comment this out for real thing
        #self.test_patient_hashes = [list(reasoner.metadata.keys())[0]]

        target_actual = '{} {} {} Actual'.format(demo_target[0], demo_target[1], demo_target[2])
        results = {'X': [demo_value for _ in range(len(self.test_patient_hashes[:10]))],
                   'Patient': list(),
                   'Compute_time': list(),
                   'Length_evid': list(),
                   target_actual: list()}

        if gene_evidence is False and demo_evidence is None:
            raise ValueError('Gene and demo evidence can not be None.')
        if demo_evidence is not None:
            title = 'Evidence = {} {} X'.format(demo_evidence[0], demo_evidence[1])
            demo_evidence_run = [tuple(list(demo_evidence) + [demo_value])]
            for patient_hash in tqdm.tqdm(self.test_patient_hashes[:10]):
                results['Patient'].append(patient_hash)
                #-- Calculate Truth
                prop_targ, op_targ, val_targ = demo_target
                op_ = _process_operator(op_targ)
                results[target_actual].append(op_(reasoner.metadata[patient_hash][prop_targ], val_targ))
                if gene_evidence:
                    gene_evidence_patient = {gene:'True' for gene in patient_data[patient_hash]['Patient_Genes']}
                    title += ' w/ GeneEvidence'
                else:
                    title += ' w/o GeneEvidence'
                    gene_evidence_patient=dict()

                #-- Make query
                query = Query(evidence=gene_evidence_patient,
                              meta_evidence=demo_evidence_run,
                              meta_targets=[demo_target],
                              type='updating')
                query = reasoner.analyze_query(copy.deepcopy(query), save_dir=None)
                results['Length_evid'].append(len(query.evidence))

                results['Compute_time'].append(query.compute_time)
                for update, prob in query.result.updates.items():
                    comp_idx, state_idx = update
                    comp_name = query.bkb.getComponentName(comp_idx)
                    state_name = query.bkb.getComponentINodeName(comp_idx, state_idx)
                    try:
                        results[('Updates', comp_name, state_name)].append(prob)
                    except:
                        results[('Updates', comp_name, state_name)] = [prob]

        else:
            title = 'No Demographic Evidence w/ Gene Evidence'
            for i, patient_hash in enumerate(self.test_patient_hashes[:10]):
                results['Patient'].append(patient_hash)
                #-- Calculate Truth
                prop_targ, op_targ, val_targ = demo_target
                op_ = _process_operator(op_targ)
                #-- Calculate Truth
                results[target_actual].append(op_(reasoner.metadata[patient_hash][prop_targ], val_targ))
                #-- Take intersection between all bkb genes and patient genes
                print('Intersecting genes')
                genes_patient = set(['mut_{}='.format(gene) for gene in reasoner.metadata[patient_hash]['Patient_Genes']])
                bkb_comps = set(self.bkb.getAllComponentNames())
                gene_intersect = genes_patient.intersection(bkb_comps)
                gene_evidence_patient = {gene: 'True' for gene in gene_intersect}
                #print(set(gene_evidence_patient.keys()) - set(self.bkb.getAllComponentNames()))
                #-- Make query
                #print(gene_evidence_patient)
                results['Length_evid'].append(len(gene_evidence_patient))
                query = Query(evidence=gene_evidence_patient,
                              meta_targets=[demo_target],
                              type='updating',
                              name='Pat50-10_600_{}'.format(i))
                print('Analyzing Query.')
                query = reasoner.analyze_query(copy.deepcopy(query), save_dir=None)
                query.getReport()
                results['Compute_time'].append(query.compute_time)
                for update, prob in query.result.updates.items():
                    comp_idx, state_idx = update
                    comp_name = query.bkb.getComponentName(comp_idx)
                    state_name = query.bkb.getComponentINodeName(comp_idx, state_idx)
                    try:
                        results[('Updates', comp_name, state_name)].append(prob)
                    except:
                        results[('Updates', comp_name, state_name)] = [prob]
                del query
                print('Complete Patient {}.'.format(i))
        return self.results_to_dataframe(results), title

    def results_to_dataframe(self, results):
        print(results)
        df = pd.DataFrame(data=results)
        df.to_csv('crossValid-results-1.csv')
        return df

if __name__ == '__main__':
    fused_bkb = BKB()
    fused_bkb.load('/home/public/data/ncats/90PERCENTValidation200Patients40Validation/set1/fusion.bkb')
    #print(fused_bkb.getAllComponentNames())
    patient_data_file = '/home/public/data/ncats/90PERCENTValidation200Patients40Validation/set1/patient_data.pk'
    test_patient_hashes = read_withheldPatientFile('/home/public/data/ncats/90PERCENTValidation200Patients40Validation/set1/withheldPatients.csv')
    #print(test_patient_hashes)
    '''
    fused_bkb = BKB()
    fused_bkb.load('10pat.bkb')
    patient_data_file = 'patient_data.pk'
    test_patient_hashes = None
    #print(test_patient_hashes)
    '''
    demo_evidence = ('Age_of_Diagnosis', '>=')
    demo_evidence_vals = [val for val in range(15000, 30000, 5000)]

    demo_target = ('Survival_Time', '>=', 300)

    cross_validator = CrossValidator(fused_bkb, test_patient_hashes, patient_data_file)
    df = cross_validator.run_demo_suite(demo_target, demo_evidence_dynamic=demo_evidence, demo_evidence_values=demo_evidence_vals)
    print(df)
