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

class DefaultHandler:
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
        for node_key in self.qg['nodes'].keys():
            node = self.qg['nodes'][node_key]
            if node['category'] == 'phenotypicfeature' and node['curie'] in acceptable_target_curies:
                target_id = node_key
                total_nodes += 1
        if total_nodes == 0:
            acceptable_target_curies_print = ','.join(acceptable_target_curies)
            sys.exit('Survival Node not found. Node type muse be \'PhenotypicFeature\' and curie must be in: ' + acceptable_target_curies_print)
        elif total_nodes > 1:
            sys.exit('Too many target nodes')

        # get disease node info and ensure only 1 disease:
        acceptable_disease_curies = ['MONDO:0007254']
        for node_key in self.qg['nodes'].keys():
            node = self.qg['nodes'][node_key]
            if node['category'] == 'disease' and node['curie'] in acceptable_disease_curies:
                disease_id = node_key
                for edge_key in self.qg['edges'].keys():
                    edge = self.qg['edges'][edge_key]
                    if edge['predicate'] == 'disease_to_phenotype_association' and edge['subject'] == disease_id and edge['object'] == target_id:
                        if 'properties' in edge.keys():
                            days = edge['properties']['days']
                            qualifier = edge['properties']['qualifier']
                        else:
                            days = 970
                            qualifier = '>='
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
        targets.append(('Survival_Time', qualifier, days))

        # get evidence
        evidence = dict()
        meta_evidence = list()
        for node_key in self.qg['nodes'].keys():
            # genes
            node = self.qg['nodes'][node_key]
            if node['category'] == 'gene':
                # check for appropriate gene node structure
                gene_id = node_key
                for edge_key in self.qg['edges'].keys():
                    edge = self.qg['edges'][edge_key]
                    if edge['predicate'] == 'gene_to_disease_association' and edge['subject'] == gene_id and edge['object'] == disease_id:
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
                evidence["mut_" + gene] = 'True'
                total_nodes += 1
            # drugs
            if node['category'] == 'drug':
                # check for appropriate drug node structure
                drug_id = node_key
                for edge_key in self.qg['edges'].keys():
                    edge = self.qg['edges'][edge_key]
                    if edge['predicate'] == 'chemical_to_disease_or_phenotypic_feature_association' and edge['subject'] == drug_id and edge['object'] == disease_id:
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

        report = query.jsonExplanations(contributions_include_srcs=False,
                                        contributions_top_n_inodes=10,
                                        contributions_ignore_prefixes=['_'])
        self.report = {'Patient Analysis': report['Patient Analysis'],
                       'Contribution Analysis': report['Contributions Analysis']}

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
        # update target node info and form edge pair combos for results graph

        node_pairs = dict()
        knowledge_nodes = 0
        for node_key in list(self.kg['nodes'].keys())[:]:
            kg_id = 'kgn{}'.format(knowledge_nodes)
            knowledge_nodes += 1
            self.kg['nodes'][kg_id] = self.kg['nodes'].pop(node_key)
            node_pairs[node_key] = kg_id
            if self.kg['nodes'][kg_id]['category'] == 'gene':
                self.kg['nodes'][kg_id]['name'] = self.gene_curie_dict[self.kg['nodes'][kg_id].pop('curie')]
            elif self.kg['nodes'][kg_id]['category'] == 'drug':
                self.kg['nodes'][kg_id]['name'] = self.drug_curie_dict[self.kg['nodes'][kg_id].pop('curie')]
            else:
                self.kg['nodes'][kg_id].pop('curie')

        edge_pairs = dict()
        knowledge_edges = 0
        for edge_key in list(self.kg['edges'].keys())[:]:
            kg_id = 'kge{}'.format(knowledge_edges)
            knowledge_edges += 1
            self.kg['edges'][kg_id] = self.kg['edges'].pop(edge_key)
            self.kg['edges'][kg_id]['subject'] = node_pairs[self.kg['edges'][kg_id]['subject']]
            self.kg['edges'][kg_id]['object'] = node_pairs[self.kg['edges'][kg_id]['object']]
            edge_pairs[edge_key] = kg_id
            if self.kg['edges'][kg_id]['predicate'] == 'disease_to_phenotype_association':
                self.kg['edges'][kg_id]['has_confidence_level'] = self.truth_assignment
                if 'properties' in self.kg['edges'][kg_id].keys() and 'contributions' in self.kg['edges'][kg_id]['properties'].keys() and self.kg['edges'][kg_id]['properties']['contributions'] == True:
                    self.kg['edges'][kg_id]['Description'] = self.report

        edge_bindings = 0
        for edge_pair_key in edge_pairs:
            self.results[0]['edge_bindings']['eb{}'.format(edge_bindings)] = { 'qg_id': edge_pair_key,
                                                                               'kg_id': edge_pairs[edge_pair_key]
                                                                             }
            edge_bindings += 1
        node_bindings = 0
        for node_pair_key in node_pairs:
            self.results[0]['node_bindings']['nb{}'.format(node_bindings)] = { 'qg_id': node_pair_key,
                                                                               'kg_id': node_pairs[node_pair_key]
                                                                             }
            node_bindings += 1

        # query response
        reasoner_std = {'query_graph': self.qg,
                        'knowledge_graph': self.kg,
                        'results': self.results}
        reasoner_std = {'message' : reasoner_std}
        return reasoner_std
