import json
import itertools
import tqdm
import numpy as np
import random
import copy
import uuid

from chp.core.query import Query
from chp.core.reasoner import Reasoner
from pybkb.core.common.bayesianKnowledgeBase import bayesianKnowledgeBase as BKB
from pybkb.core.common.bayesianKnowledgeBase import BKB_I_node, BKB_component, BKB_S_node
from pybkb.core.python_base.fusion import fuse

class UnsecretHandler:
    def __init__(self, query):
        self.query = query
        self.qg = self.query['query_graph']
        self.kg = self.query['knowledge_graph']
        self.results = self.query['results']
        self.probability_targets = self.query['probability_targets']

        self.fused_bkb = BKB()
        #self.fused_bkb.load('/home/ghyde/bkb-pathway-provider/core/testOldStuffSmall/fusion.bkb')
        #self.patient_data_file = '/home/ghyde/bkb-pathway-provider/core/testOldStuffSmall/patient_data.pk'
        self.fused_bkb.load('/home/public/data/ncats/BabelBKBs/collapsedAll/fusion.bkb')
        self.patient_data_file = '/home/public/data/ncats/BabelBKBs/collapsedAll/patient_data.pk'
        self.reasoner = Reasoner(self.fused_bkb, None)
        self.reasoner.set_src_metadata(self.patient_data_file)
        self.reasoner.cpp_reasoning = False

        self.target_strategy = 'chain'
        self.interpolation = 'standard'

#        #-- Collect Node Bindings
#        self.node_bindings_map = {}
#        for binding in self.results['node_bindings']:
#            for qnode_id, knode_id in binding.items():
#                if knode_id not in self.node_bindings_map:
#                    self.node_bindings_map[knode_id] = set([qnode_id])
#                else:
#                    self.node_bindings_map[knode_id].add(qnode_id)
#
#        #-- Hash edges and nodes by ID
#        self.kedges = {edge['id']: edge for edge in self.kg['edges']}
#        self.qedges = {edge['id']: edge for edge in self.qg['edges']}
#        self.knodes = {node['id']: node for node in self.kg['nodes']}
#        self.qnodes = {node['id']: node for node in self.qg['nodes']}

    def checkQuery(self):
        return True

    def buildQueries(self):
        """ Assumes Query Graph of structure: Disease -> |Intermediate Nodes, i.e. Genes, Drugs, etc.| -> UpdateTargetNode
            - Update Target Node is designated by a 't' label in the QG, i.e. t1, t2, etc.
            - The beginning disease node tells us which datasets to use.
            - The intermediate nodes dictate the evidence to use during updating.
        """
        queries = []
        evidence = dict()
        for node in self.qg['nodes']:
            if node['type'] == 'Gene':
                evidence["mut_" + node['name']] = 'True'
        query = Query(evidence=evidence,
                      targets=[],
                      meta_evidence=None,
                      meta_targets=self.probability_targets,
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
            if prob == -1:
                prob = 0
            self.target_info.append([comp_name, state_name, prob])
        prob_sum = self.target_info[0][2] + self.target_info[1][2]
        self.target_info[0][2] /= prob_sum
        self.target_info[1][2] /= prob_sum

        report = query.jsonExplanations()
        self.patient_report = report['Patient Analysis']

    def constructDecoratedKG(self):
        self.kg = copy.deepcopy(self.qg)
        results = {'edge_bindings':[], 'node_bindings':[]}
        Q_graph_node_target = None
        K_graph_new_node_target = None
        # update target node info and form new KGraph id
        for node in self.kg['nodes']:
            if 'curie' in node.keys() and node['curie'] == 'CHPDART:SURVIVAL':
                node['has_confidence_level'] = self.target_info[0][2]
                node['Description'] = self.patient_report
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
