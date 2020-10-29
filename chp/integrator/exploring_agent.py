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

class ExploringHandler:
    def __init__(self, query, hosts_filename=None, num_processes_per_host=0):
        # query graph components
        self.query = query
        self.bkb_data_handler = BkbDataHandler(dataset_version='1.1')

        # Only do the rest of this if a query is passed
        if self.query is not None:
            self.qg = self.query['query_graph']
            if 'knowledge_graph' not in list(self.query.keys()):
                self.kg = { "edges": [],
                            "nodes": []
                          }
            else:
                self.kg = self.query['knowledge_graph']
            if 'results' not in list(self.query.keys()):
                self.results = { "node_bindings": [],
                                 "edge_bindings": []
                               }
            else:
                self.results = self.query['results']

            # Instiatate Reasoner
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
    # Output: built query
    #--------------------------------------------------------
    # Description: Parses over sent query graph to form a BKB
    # query. The BKB query is returned.

    def buildQueries(self):

        # ensure we are using all nodes/edges
        total_nodes = 0
        total_edges = 0

        # get phenotype node
        targets = list()
        acceptable_target_curies = ['EFO:0000714']
        for node in self.qg['nodes']:
            if node['type'] == 'PhenotypicFeature' and node['curie'] in acceptable_target_curies:
                target_id = node['id']
                total_nodes += 1
        if total_nodes == 0:
            acceptable_target_curies_print = ','.join(acceptable_target_curies)
            sys.exit('Survival Node not found. Node type muse be \'PhenotypicFeature\' and curie must be in: ' + acceptable_target_curies_print)
        elif total_nodes > 1:
            sys.exit('Too many target nodes')

        # get disease node info and ensure only 1 disease:
        acceptable_disease_curies = ['MONDO:0007254']
        for node in self.qg['nodes']:
            if node['type'] == 'disease' and node['curie'] in acceptable_disease_curies:
                disease_id = node['id']
                for edge in self.qg['edges']:
                    if edge['type'] == 'disease_to_phenotype_association' and edge['source_id'] == disease_id and edge['target_id'] == target_id and 'value' in list(edge.keys()):
                        days = edge['value']
                        if isinstance(days, str):
                            days = int(days)
                        total_edges += 1
                if total_edges == 0:
                    sys.exit('Disease and target edge not found. Edge type must be \'disease_to_phenotype_association\'')
                elif total_edges > 1:
                    sys.exit('Disease has too many outgoing edges')
                total_nodes += 1
        if total_nodes  == 1:
            acceptable_disease_curies_print = ','.join(acceptable_disease_curies)
            sys.exit('Disease node not found. Node type must be \'disease\' and curie must be in: ' + acceptable_disease_curies_print)
        elif total_nodes > 2:
            sys.exit('Too many disease nodes')
        # set BKB target
        targets.append(('Survival_Time', '>=', days))

        # get evidence
        evidence = dict()
        meta_evidence = list()
        for node in self.qg['nodes']:
            # genes
            if node['type'] == 'Gene':
                # check for appropriate gene node structure
                gene_id = node['id']
                for edge in self.qg['edges']:
                    if edge['type'] == 'gene_to_disease_association' and edge['source_id'] == gene_id and edge['target_id'] == disease_id:
                        total_edges += 1
                if total_edges == total_nodes - 1:
                    sys.exit('Gene and disease edge not found. Edge type must be \'gene_to_disease_association\'')
                elif total_edges > total_nodes:
                    sys.exit('Gene has too many outgoing edges')
                # check for appropriate gene node curie
                gene_curie = node['curie']
                try:
                    gene = self.gene_curie_dict[gene_curie]
                except:
                    sys.exit('Invalid ENSEMBL Identifier. Must be in form ENSEMBL:<ID>.')
                evidence["_mut_" + gene] = 'True'
                total_nodes += 1
            # drugs
            if node['type'] == 'Drug':
                # check for appropriate drug node structure
                drug_id = node['id']
                for edge in self.qg['edges']:
                    if edge['type'] == 'chemical_to_disease_or_phenotypic_feature_association' and edge['source_id'] == drug_id and edge['target_id'] == disease_id:
                        total_edges += 1
                if total_edges == total_nodes - 1:
                    sys.exit('Drug and disease edge not found. Edge type must be \'chemical_to_disease_or_phenotypic_feature_association\'')
                elif total_edges > total_nodes:
                    sys.exit('Drug has too many outgoing edges')
                # check for appropriate drug node curie
                drug_curie = node['curie']
                try:
                    drug = self.drug_curie_dict[drug_curie]
                except:
                    sys.exit('Invalid CHEMBL Identifier. Must be in form CHEMBL:<ID>')
                meta_evidence.append(('Drug_Name(s)', '==', drug))
                total_nodes += 1

        if total_nodes != len(self.qg['nodes']) or total_edges != len(self.qg['edges']):
            sys.exit('There are extra components in the provided QG structure')

        # produce BKB query
        if len(meta_evidence) > 0 and len(list(evidence.keys())) > 0:
            query = Query(evidence=evidence,
                          targets=[],
                          meta_evidence=meta_evidence,
                          meta_targets=targets,
                          type='updating')
        elif len(list(evidence.keys())) > 0:
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
            query = Query(evidence=None,
                          targets=[],
                          meta_evidence=None,
                          meta_targets=targets,
                          type='updating')

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
        results = {'edge_bindings':[], 'node_bindings':[]}
        # update target node info and form edge pair combos for results graph
        edge_pairs = list()
        for edge in self.kg['edges']:
            if edge['type'] == 'disease_to_phenotype_association':
                edge['has_confidence_level'] = self.truth_assignment
                #they may want descriptions later
                #node['Description'] = self.patient_report
            qg_id = edge['id']
            kg_id = str(uuid.uuid4())
            edge_pairs.append([qg_id,kg_id])
            edge['id'] = kg_id
        # form node pair combos for results graph
        node_pairs = list()
        for node in self.kg['nodes']:
            qg_id = node['id']
            kg_id = str(uuid.uuid4())
            node_pairs.append([qg_id,kg_id])

        # update edge/node bindings in results graph
        for edge_pair in edge_pairs:
            results['edge_bindings'].append({'qg_id':edge_pair[0], 'kg_id':edge_pair[1]})
        for node_pair in node_pairs:
            results['node_bindings'].append({'qg_id':node_pair[0], 'kg_id':node_pair[1]})

        # query response
        reasoner_std = {'query_graph': self.qg,
                        'knowledge_graph': self.kg,
                        'results': results}
        return reasoner_std
