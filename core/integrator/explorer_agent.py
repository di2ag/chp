import json
import itertools
import tqdm
import numpy as np
import random
import copy
import uuid
import os

from chp.core.query import Query
from chp.core.reasoner import Reasoner
from pybkb.core.common.bayesianKnowledgeBase import bayesianKnowledgeBase as BKB
from pybkb.core.common.bayesianKnowledgeBase import BKB_I_node, BKB_component, BKB_S_node
from pybkb.core.python_base.fusion import fuse

ABSOLUTE_PATH_TO_BKB = '/home/public/data/ncats/BabelBKBs_7-7-2020/All'

class ExplorerHandler:
    def __init__(self, query):
        self.query = query
        self.qg = self.query['query_graph']
        self.kg = self.query['knowledge_graph']
        self.results = self.query['results']
        self.probability_targets = self.query['probability_targets']

        #-- Collect Node Bindings
        self.node_bindings_map = {}
        for binding in self.results['node_bindings']:
            for qnode_id, knode_id in binding.items():
                if knode_id not in self.node_bindings_map:
                    self.node_bindings_map[knode_id] = set([qnode_id])
                else:
                    self.node_bindings_map[knode_id].add(qnode_id)

        #-- Hash edges and nodes by ID
        self.kedges = {edge['id']: edge for edge in self.kg['edges']}
        self.qedges = {edge['id']: edge for edge in self.qg['edges']}
        self.knodes = {node['id']: node for node in self.kg['nodes']}
        self.qnodes = {node['id']: node for node in self.qg['nodes']}

        #-- Instiatate Reasoner
        self.fused_bkb = BKB()
        self.fused_bkb.load('/home/public/data/ncats/BabelBKBs_7-7-2020/All/fusion.bkb')
        self.patient_data_file = '/home/public/data/ncats/BabelBKBs_7-7-2020/All/patient_data.pk'
        self.reasoner = Reasoner(self.fused_bkb, None)
        self.reasoner.set_src_metadata(self.patient_data_file)
        self.reasoner.cpp_reasoning = False

        self.target_strategy = 'explicit'
        self.interpolation = 'standard'

    def checkQuery(self):
        return True

    def getOptions(self, reverse=False):
        """ Will return a dictionary with keys corresponding to QG edges mapped to a
        KG edges based on results."""
        options = dict()
        for edge in self.qg['edges']:
            options[edge['id']] = []
        for result in self.results['edge_bindings']:
            for qg_id, kg_ids in result.items():
                options[qg_id].append(kg_ids[0])
        if reverse:
            reverse_options = {}
            for key, option_list in options.items():
                for kedge in option_list:
                    reverse_options[kedge] = key
            return options, reverse_options
        return options

    def getAnswers(self, options=None):
        """ Differs from results in that it takes the results and builds consistent
        complete subgraphs of the QG, i.e. an answer."""
        #-- Get answers defined as a set of edges such that all nodes are consistent, i.e. complete subgraph that answers query.
        answers = []
        if options is None:
            options = self.getOptions()
        #-- Only use non-target edges for options
        options = {edge: values for edge, values in options.items() if edge[:2] != 'et'}
        total = np.prod([len(option) for option in options.values()])
        for answer in tqdm.tqdm(iterable=itertools.product(*list(options.values())), total=total):
            for i, kg_id in enumerate(answer):
                source_id = self.kedges[kg_id]['source_id']
                target_id = self.kedges[kg_id]['target_id']
                if i == 0:
                    prev_target_id = target_id
                    continue
                elif prev_target_id != source_id:
                    break
                elif i == len(answer) - 1:
                    answers.append(answer)
        return answers

    def extractUpdateTargetEdges(self, options=None):
        if options is None:
            options = self.getOptions()

        #-- Find target nodes
        target_edges = []
        unique_targets = set()
        for edge_id in options:
            if edge_id[:2] == 'et':
                unique_targets.add(self.qedges[edge_id]['target_id'])
                target_edges.append(edge_id)
        return target_edges, unique_targets

    def buildQueries(self):
        """ Assumes Query Graph of structure: Disease -> |Intermediate Nodes, i.e. Genes, Drugs, etc.| -> UpdateTargetNode
            - Update Target Node is designated by a 't' label in the QG, i.e. t1, t2, etc.
            - The beginning disease node tells us which datasets to use.
            - The intermediate nodes dictate the evidence to use during updating.
        """
        options, reverse_options = self.getOptions(reverse=True)
        answers = self.getAnswers(options=options)
        update_target_edges, _ = self.extractUpdateTargetEdges(options=options)

        queries = []
        for answer in answers:
            #-- Get update target
            targets = []
            '''
            meta_targets = []
            target_source_id = None
            for edge in update_target_edges:
                target_id = self.qedges[edge]['target_id']
                target_qnode = self.qnodes[target_id]
                ont, feature = target_qnode['curie'].split(':')
                if target_qnode['type'] == 'PhenotypicFeature':
                    if ont == 'CHPDART':
                        if feature == 'SURVIVAL':
                            meta_targets.append(('Survival_Time',
                                                 target_qnode['operator'],
                                                 target_qnode['value']))
                        else:
                            raise NotImplementedError
                    else:
                        raise NotImplementedError
                else:
                    raise NotImplementedError
            #-- Collect Evidence starting from target edges
            answer_edges_by_target_id = {self.kedges[edge]['target_id']: self.kedges[edge] for edge in answer}
            '''
            evidence = {}
            meta_evidence = []
            disease = None
            processed_node_ids = set()
            for edge in answer:
                target_node = self.knodes[self.kedges[edge]['target_id']]
                source_node = self.knodes[self.kedges[edge]['source_id']]
                for node in [target_node, source_node]:
                    if node['id'] in processed_node_ids:
                        continue
                    #-- Add evidence
                    if node['type'] == 'ChemicalSubstance':
                        meta_evidence.append(('Patient_Drug(s)', '==', node['id']))
                    elif node['type'] == 'Gene':
                        evidence['_mut_{}'.format(node['id'])] = 'True'
                    elif node['type'] == 'Disease':
                        disease = node['name']
                        if disease != 'BREAST CANCER':
                            raise NotImplementedError
                    else:
                        raise NotImplementedError
                    processed_node_ids.add(node['id'])
            queries.append(Query(evidence=evidence,
                                 meta_evidence=meta_evidence,
                                 targets=targets,
                                 meta_targets=self.probability_targets))
        self.answers = answers
        self.chp_queries = queries
        return queries

    def runQueries(self, bogus=True):
        queries = copy.deepcopy(self.chp_queries)

        updated_queries = []
        if bogus:
            self.bogus_updates = True
            for query in queries:
                updated_queries.append(query.make_bogus_updates())
        else:
            self.bogus_updates = False
            for query in queries:
                updated_queries.append(self.reasoner.analyze_query(copy.deepcopy(query),
                                            save_dir=None,
                                            target_strategy=self.target_strategy,
                                            interpolation=self.interpolation))
        self.updated_queries = updated_queries
        return queries

    def constructDecoratedKG(self):
        options, reverse_options = self.getOptions(reverse=True)
        target_edges, unique_targets = self.extractUpdateTargetEdges(options=options)

        #-- Going to place nodes in the KG for the probability targets but not in the query graph.

        inserted_updates = set()
        results = {'edge_bindings': [], 'node_bindings': [], 'probability_target_edges': []}
        for answer, query in zip(self.answers, self.updated_queries):
            #-- Construct update target nodes
            node_bindings = {}
            edge_bindings = {}
            probability_target_edges = []
            if self.bogus_updates:
                processed_updates = query
            else:
                processed_updates = query.result.process_updates(normalize=True)
            for comp_name, state_prob_dict in processed_updates.items():
                for state_name, prob in state_prob_dict.items():
                    comp_state_id = str(uuid.uuid4())
                    self.kg['nodes'].append({'id': comp_state_id,
                                                          'name': '{} = {}'.format(comp_name, state_name),
                                                          'type': 'probability',
                                                          'node_attributes': [{'type': 'float',
                                                                               'name': 'probability',
                                                                               'value': str(prob)}]
                                                         })
                    '''
                    if update_target_qnode not in node_bindings:
                        node_bindings[update_target_qnode] = [update_id]
                    else:
                        node_bindings[update_target_qnode].append(update_id)
                    '''
                    processed_knode_ids = set()
                    for kedge in answer:
                        qedge = reverse_options[kedge]
                        target_knode_id = self.kedges[kedge]['target_id']
                        source_knode_id = self.kedges[kedge]['source_id']
                        target_qnode_id = self.qedges[qedge]['target_id']
                        source_qnode_id = self.qedges[qedge]['source_id']
                        #-- create standard edge binding
                        edge_bindings[qedge] = [kedge]

                        for knode_id, qnode_id in zip([target_knode_id, source_knode_id], [target_qnode_id, source_qnode_id]):
                            if knode_id in processed_knode_ids:
                                continue
                            processed_knode_ids.add(knode_id)
                            new_kg_edge_id = str(uuid.uuid4())

                            #-- Make node bindings
                            if qnode_id not in node_bindings:
                                node_bindings[qnode_id] = [knode_id]
                            else:
                                node_bindings[qnode_id].append(knode_id)
                            #-- make edge from each node in the kg to the update comp_state_node
                            self.kg['edges'].append({'id': new_kg_edge_id,
                                                     'type': 'affects',
                                                     'relation': 'correlates_with',
                                                     'source_id': knode_id,
                                                     'target_id': comp_state_id})
                            probability_target_edges.append(new_kg_edge_id)
                            ''''
                            #-- Find associated edge in QG for new edge in KG.
                            found_qg_edge = None
                            for _qedge in self.qg['edges']:
                                if _qedge['source_id'] == qnode_id and _qedge['target_id'] == update_target_qnode:
                                    found_qg_edge = _qedge
                                    break
                            #-- Make node and edge bindings
                            if qnode_id not in node_bindings:
                                node_bindings[qnode_id] = [knode_id]
                            else:
                                node_bindings[qnode_id].append(knode_id)
                            if found_qg_edge['id'] not in edge_bindings:
                                edge_bindings[found_qg_edge['id']] = [new_kg_edge_id]
                            else:
                                edge_bindings[found_qg_edge['id']].append(new_kg_edge_id)
                            '''
                        #- Add other edge bindings
                        #edge_bindings[reverse_options[kedge]] = [kedge]

            #-- Add to new results
            results['edge_bindings'].append(edge_bindings)
            results['node_bindings'].append(node_bindings)
            results['probability_target_edges'].append(probability_target_edges)
        reasoner_std = {'query_graph': self.qg,
                        'knowledge_graph': self.kg,
                        'results': results}
        return reasoner_std
