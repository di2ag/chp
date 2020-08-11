'''
Source code developed by DI2AG.
Thayer School of Engineering at Dartmouth College
	Authors:    Dr. Eugene Santos, Jr
            Mr. Chase Yakaboski,
            Mr. Gregory Hyde,
            Dr. Keum Joo Kim
'''


import json
import itertools
import tqdm
import numpy as np
import random
import copy
import uuid
import pickle
import csv
import sys

from chp.query import Query
from chp.reasoner import Reasoner

from chp_data.bkb_handler import BkbDataHandler

from pybkb.common.bayesianKnowledgeBase import bayesianKnowledgeBase as BKB
from pybkb.common.bayesianKnowledgeBase import BKB_I_node, BKB_component, BKB_S_node
from pybkb.python_base.fusion import fuse

class RankingHandler:
    def __init__(self, query, hosts_filename=None, num_processes_per_host=0):
        self.query = query
        self.qg = self.query['query_graph']
        self.kg = self.query['knowledge_graph']
        self.results = self.query['results']


        #change
        #-- Instiatate Reasoner
        self.bkb_data_handler = BkbDataHandler()
        self.reasoner = Reasoner(bkb_data_handler=self.bkb_data_handler,
                                hosts_filename=hosts_filename,
                                num_processes_per_host=num_processes_per_host)
        self.gene_curie_dict = dict()
        with open(self.bkb_data_handler.gene_curie_path, 'r') as gene_file:
            reader = csv.reader(gene_file)
            next(reader)
            for row in reader:
                self.gene_curie_dict[row[1]] = row[0]
        self.drug_curie_dict = dict()
        with open(self.bkb_data_handler.drug_curie_path, 'r') as drug_file:
            reader = csv.reader(drug_file)
            next(reader)
            for row in reader:
                self.drug_curie_dict[row[1]] = row[0]


        self.target_strategy = 'chain'
        self.interpolation = 'standard'

    def checkQuery(self):
        return True

    def buildQueries(self):
        """ Assumes Query Graph of structure: Disease -> |Intermediate Nodes, i.e. Genes, Drugs, etc.| -> UpdateTargetNode
            - Update Target Node is designated by a 't' label in the QG, i.e. t1, t2, etc.
            - The beginning disease node tells us which datasets to use.
            - The intermediate nodes dictate the evidence to use during updating.
        """
        queries = []

        # get evidence
        evidence = dict()
        meta_evidence = list()
        for node in self.qg['nodes']:
            if node['type'] == 'Gene':
                gene_curie = node['curie']
                try:
                    gene = self.gene_curie_dict[gene_curie]
                except:
                    sys.exit('Invalid ENSEMBL Identifier. Must be in form ENSEMBL:<ID>.')
                evidence["mut_" + gene] = 'True'
            if node['type'] == 'Drug':
                drug_curie = node['curie']
                try:
                    drug = self.drug_curie_dict[drug_curie]
                except:
                    sys.exit('Invalid CHEMBL Identifier. Must be in form CHEMBL:<ID>')
                meta_evidence.append(('Drug_Name(s)', '==', drug))

        # get target
        targets = list()
        for node in self.qg['nodes']:
            if node['type'] == 'Death':
                operator = node['operator']
                val = node['value']
                targets.append(('Survival_Time', operator, int(val)))

        if len(list(evidence.keys())) + len(meta_evidence) > 1:
            sys.exit('More than 1 piece of evidence')

        if len(list(evidence.keys())) > 0:
            query = Query(evidence=evidence,
                          targets=[],
                          meta_evidence=None,
                          meta_targets=targets,
                          type='updating')
        else:
            query = Query(evidence=None,
                          targets=[],
                          meta_evidence=meta_evidence,
                          meta_targets=targets,
                          type='updating')
        queries.append(query)
        self.chp_query = query
        return queries

    def runQueries(self, bogus=True):
        query = self.reasoner.analyze_query(copy.deepcopy(self.chp_query),
                                            save_dir=None,
                                            target_strategy=self.target_strategy,
                                            interpolation=self.interpolation
                                            )
        self.target_info = []
        for update, prob in query.result.updates.items():
            comp_idx, state_idx = update
            comp_name = query.bkb.getComponentName(comp_idx)
            state_name = query.bkb.getComponentINodeName(comp_idx, state_idx)
            #print(comp_name, state_name, prob)
            self.target_info.append([comp_name, state_name, prob])
        if self.target_info[0][2] != -1 and self.target_info[1][2] != -1:
            prob_sum = self.target_info[0][2] + self.target_info[1][2]
            self.target_info[0][2] /= prob_sum
            self.target_info[1][2] /= prob_sum
        elif self.target_info[0][2] == -1 and self.target_info[1][2] != -1:
            self.target_info[0][2] = 0
            prob_sum = self.target_info[0][2] + self.target_info[1][2]
            self.target_info[0][2] /= prob_sum
            self.target_info[1][2] /= prob_sum
        elif self.target_info[0][2] != -1 and self.target_info[1][2] == -1:
            self.target_info[1][2] == 0
            prob_sum = self.target_info[0][2] + self.target_info[1][2]
            self.target_info[0][2] /= prob_sum
            self.target_info[1][2] /= prob_sum

        #report = query.jsonExplanations()
        #self.patient_report = report['Patient Analysis']

    def constructDecoratedKG(self):
        self.kg = copy.deepcopy(self.qg)
        results = {'edge_bindings':[], 'node_bindings':[]}
        Q_graph_node_target = None
        K_graph_new_node_target = None
        # update target node info and form new KGraph id
        for node in self.kg['nodes']:
            if node['curie'] == 'UBERON:0000071':
                node['has_confidence_level'] = self.target_info[0][2]
                #they may want descriptions later
                #node['Description'] = self.patient_report
                Q_graph_node_target = node['id']
                K_graph_new_node_target = str(uuid.uuid4())
                node['id'] = K_graph_new_node_target

        Q_graph_edge_target = None
        K_graph_new_edge_target = None
        # update KGraph edge connection to target with new KGraph id
        for edge in self.kg['edges']:
            if edge['target_id'] == Q_graph_node_target:
                edge['target_id'] == K_graph_new_node_target
                Q_graph_edge_target = edge['id']
                K_graph_new_edge_target = str(uuid.uuid4())
                edge['id'] = K_graph_new_edge_target

        results['node_bindings'].append({'qg_id':Q_graph_node_target, 'kg_id':K_graph_new_node_target})
        results['edge_bindings'].append({'qg_id':Q_graph_edge_target, 'kg_id':K_graph_new_edge_target})

        reasoner_std = {'query_graph': self.qg,
                        'knowledge_graph': self.kg,
                        'results': results}
        return reasoner_std
