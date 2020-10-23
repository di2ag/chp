'''
Source code developed by DI2AG.
Thayer School of Engineering at Dartmouth College
	Authors:    Dr. Eugene Santos, Jr
            Mr. Chase Yakaboski,
            Mr. Gregory Hyde,
            Dr. Keum Joo Kim
'''

import copy
import csv
import sys
import uuid

from chp.query import Query
from chp.reasoner import Reasoner

from chp_data.bkb_handler import BkbDataHandler

#from pybkb.common.bayesianKnowledgeBase import bayesianKnowledgeBase as BKB
#from pybkb.common.bayesianKnowledgeBase import BKB_I_node, BKB_component, BKB_S_node
#from pybkb.python_base.fusion import fuse

class RankingHandler:
    def __init__(self, query, hosts_filename=None, num_processes_per_host=0):
        # query graph components
        self.query = query
        self.qg = self.query['query_graph']
        if 'knowledge_graph' not in list(self.query.keys()):
            self.kg = { "edges": dict(),
                        "nodes": dict()
                      }
        else:
            self.kg = self.query['knowledge_graph']
        if 'results' not in list(self.query.keys()):
            self.results = [{ "node_bindings": dict(),
                              "edge_bindings": dict()
                           }]
        else:
            self.results = self.query['results']

        # Instiatate Reasoner
        self.bkb_data_handler = BkbDataHandler(dataset_version='1.1')
        self.reasoner = Reasoner(bkb_data_handler=self.bkb_data_handler,
                                hosts_filename=hosts_filename,
                                num_processes_per_host=num_processes_per_host)
        # prepare curie gene dict
        self.gene_curie_dict = dict()
        with open(self.bkb_data_handler.gene_curie_path, 'r') as gene_file:
            reader = csv.reader(gene_file)
            next(reader)
            for row in reader:
                self.gene_curie_dict[row[1]] = row[0]
        # prepare curie drug dict
        self.drug_curie_dict = dict()
        with open(self.bkb_data_handler.drug_curie_path, 'r') as drug_file:
            reader = csv.reader(drug_file)
            next(reader)
            for row in reader:
                self.drug_curie_dict[row[1]] = row[0]

        # default query specification
        self.target_strategy = 'explicit'
        self.interpolation = 'standard'

    ##########################################################
    # checkQuery
    # Input:
    # Output: boolean True/False
    #--------------------------------------------------------
    # Description: NOT IMPLEMENTED. Checks the query graph
    # to ensure it meets our specifications

    def checkQuery(self):
        return True

    ##########################################################
    # buildQueries
    # Input:
    # Output: list of queries
    #--------------------------------------------------------
    # Description: Parses over sent query graph to form a BKB
    # query. The BKB query is returned. Only allows 1 gene curie
    # or 1 drug curie. Raises error if both or non. There is an
    # optional edge property 'onset_qualifier'. In the event this
    # is not provided we default to 970 days. This needs to be
    # documented in openAPI documentation?

    def buildQueries(self):
        # get evidence
        evidence = dict()
        meta_evidence = list()
        for node_id in self.qg['nodes'].keys():
            if self.qg['nodes'][node_id]['type'] == 'Gene':
                gene_curie = self.qg['nodes'][node_id]['curie']
                try:
                    gene = self.gene_curie_dict[gene_curie]
                except:
                    sys.exit('Invalid ENSEMBL Identifier. Must be in form ENSEMBL:<ID>.')
                evidence["_mut_" + gene] = 'True'
            if self.qg['nodes'][node_id]['type'] == 'Drug':
                drug_curie = self.qg['nodes'][node_id]['curie']
                try:
                    drug = self.drug_curie_dict[drug_curie]
                except:
                    sys.exit('Invalid CHEMBL Identifier. Must be in form CHEMBL:<ID>')
                meta_evidence.append(('Drug_Name(s)', '==', drug))

        # get target
        targets = list()
        for edge_id in self.qg['edges'].keys():
            if self.qg['edges'][edge_id]['type'] == 'causes':
                if 'onset_qualifier' in list(self.qg['edges'][edge_id].keys()):
                    days = self.qg['edges'][edge_id]['onset_qualifier']
                # default value - needs to be documented in openAPI?
                else:
                    days = 970
                targets.append(('Survival_Time', '>=', days))
        if len(list(evidence.keys())) + len(meta_evidence) > 1:
            sys.exit('More than 1 piece of evidence')

        if len(list(evidence.keys())) + len(meta_evidence) == 0:
            sys.exit('No evidence provided')

        if len(list(evidence.keys())) > 0:
            query = Query(evidence=evidence,
                          targets=[],
                          meta_evidence=None,
                          meta_targets=targets,
                          type='updating')
        elif len(meta_evidence) > 0:
            query = Query(evidence=None,
                          targets=[],
                          meta_evidence=meta_evidence,
                          meta_targets=targets,
                          type='updating')
        else:
            sys.exit('Problem in constructing BKB query')

        self.chp_query = query
        return query

    ##########################################################
    # runQueries
    # Input:
    # Output:
    #--------------------------------------------------------
    # Description: Runs built BKB query to calculate probability
    # of survival. Outputs are checked for -1s. If BOTH are -1
    # they are left alone (should interpret as no inference and
    # should document in openAPI). Otherwise probabilities are
    # normalized

    def runQueries(self):
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

        if self.target_info[0][1] == 'True':
            self.truth_assignment = self.target_info[0][2]
            self.false_assignment = self.target_info[1][2]
        else:
            self.truth_assignment = self.target_info[1][2]
            self.false_assignment = self.target_info[0][2]

        if self.truth_assignment != -1 and self.false_assignment != -1:
            prob_sum = self.truth_assignment + self.false_assignment
            self.truth_assignment /= prob_sum
            self.false_assignment /= prob_sum
            sensitivities = True
        elif self.truth_assignment == -1 and self.false_assignment != -1:
            self.truth_assignment = 0
            prob_sum = self.truth_assignment + self.false_assignment
            self.truth_assignment /= prob_sum
            self.false_assignment /= prob_sum
        elif self.truth_assignment != -1 and self.false_assignment == -1:
            self.false_assignment = 0
            prob_sum = self.truth_assignment + self.false_assignment
            self.truth_assignment /= prob_sum
            self.false_assignment /= prob_sum

    ##########################################################
    # constructDecoratedKG
    # Input:
    # Output: endpoint query response
    #--------------------------------------------------------
    # Description: knowledge graph is a copy of query graph where
    # we add a edge property, 'has_confidence_level' annotated with
    # our determined value. The whole query response is formed and
    # returned.

    def constructDecoratedKG(self):
        self.kg = copy.deepcopy(self.qg)
        # update target node info and form new KGraph id
        edge_ids = list()
        for edge_id in list(self.kg['edges'].keys())[:]:
            if self.kg['edges'][edge_id]['type'] == 'causes':
                self.kg['edges'][edge_id]['has_confidence_level'] = self.truth_assignment
                #they may want descriptions later
                #node['Description'] = self.patient_report
            kg_id = str(uuid.uuid4())
            self.kg['edges'][kg_id] = self.kg['edges'].pop(edge_id)
            edge_ids.append(kg_id)
        # form node pair combos for results graph
        node_ids = list()
        for node_id in list(self.kg['nodes'].keys())[:]:
            kg_id = str(uuid.uuid4())
            self.kg['nodes'][kg_id] = self.kg['nodes'].pop(node_id)
            self.kg['nodes'][kg_id]['attributes'] = [dict()]
            for node_type in list(self.kg['nodes'][kg_id].keys())[:]:
                if node_type != 'name' and node_type != 'attributes':
                    print(self.kg['nodes'][kg_id])
                    print()
                    self.kg['nodes'][kg_id]['attributes'][0][node_type] = self.kg['nodes'][kg_id].pop(node_type)
            print(self.kg['nodes'][kg_id])
            node_ids.append(kg_id)
        # set up node/edge bindings
        edge_count = 0
        for edge_id in edge_ids:
            self.results[0]['edge_bindings']['e{}'.format(edge_count)] = [{'kg_id':edge_id}]
            edge_count += 1
        node_count = 0
        for node_id in node_ids:
            self.results[0]['node_bindings']['n{}'.format(node_count)] = [{'kg_id':node_id}]
            node_count += 1

        # query response
        reasoner_std = {'query_graph': self.qg,
                        'knowledge_graph': self.kg,
                        'results': self.results}
        return {'message':reasoner_std}
