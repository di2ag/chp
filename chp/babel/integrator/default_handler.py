"""
    Source code developed by DI2AG.
    Thayer School of Engineering at Dartmouth College
    Authors:    Dr. Eugene Santos, Jr
                Mr. Chase Yakaboski,
                Mr. Gregory Hyde,
                Mr. Luke Veenhuis,
                Dr. Keum Joo Kim
"""

import copy
import csv
import sys
import uuid

from chp.query import Query
from chp.reasoner import Reasoner

from chp_data.bkb_handler import BkbDataHandler

class DefaultHandler:
    """WildCardHandler is the handler for gene wildcards. That is
        query graphs (QGs) that consists of 4 nodes and 3 edges.

        :param query: the query graph sent by the ARA.
        :type query: dict
        :param hosts_filename: a filename for a stored QG. Defaults to None
        :type hosts_filename: str
        :param num_processes_per_host: Not implemented thouroughly, but would be
            used for distributed reasoning.
        :type num_processes_per_host: int
        :param max_results: specific to 1-hop queries, specifies the number of
            wildcard genes to return.
        :type max_results: int
    """

    def __init__(self, query, hosts_filename=None, num_processes_per_host=0):
        # query graph components
        self.query = query
        self.bkb_data_handler = BkbDataHandler(dataset_version='1.3')

        # Only do the rest of this if a query is passed
        if self.query is not None:
            self.qg = self.query['query_graph']
            if 'knowledge_graph' not in list(self.query.keys()):
                self.kg = { "edges": dict(),
                            "nodes": dict()
                          }
            else:
                self.kg = self.query['knowledge_graph']
            if 'results' not in list(self.query.keys()):
                self.results = []
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

    def checkQuery(self):
        """ Currently not implemented. Would check validity of query.
        """
        return True

    def buildQueries(self):
        """ Parses over the sent query graph to form a BKB query.

            :return: A  internal CHP query.
            :rtype: Query
        """

        # ensure we are using all nodes/edges
        total_nodes = 0
        total_edges = 0

        # get phenotype node
        targets = list()
        acceptable_target_curies = ['EFO:0000714']
        for node_key in self.qg['nodes'].keys():
            node = self.qg['nodes'][node_key]
            if node['category'] == 'biolink:PhenotypicFeature' and node['id'] in acceptable_target_curies:
                target_id = node_key
                total_nodes += 1
        if total_nodes == 0:
            acceptable_target_curies_print = ','.join(acceptable_target_curies)
            sys.exit('Survival Node not found. Node category must be \'biolink:PhenotypicFeature\' and id must be in: ' + acceptable_target_curies_print)
        elif total_nodes > 1:
            sys.exit('Too many target nodes')

        # get disease node info and ensure only 1 disease:
        acceptable_disease_curies = ['MONDO:0007254']
        for node_key in self.qg['nodes'].keys():
            node = self.qg['nodes'][node_key]
            if node['category'] == 'biolink:Disease' and node['id'] in acceptable_disease_curies:
                disease_id = node_key
                for edge_key in self.qg['edges'].keys():
                    edge = self.qg['edges'][edge_key]
                    if edge['predicate'] == 'biolink:DiseaseToPhenotypicFeatureAssociation' and edge['subject'] == disease_id and edge['object'] == target_id:
                        if 'properties' in edge.keys():
                            days = edge['properties']['days']
                            qualifier = edge['properties']['qualifier']
                        else:
                            days = 970
                            qualifier = '>='
                        total_edges += 1
                if total_edges == 0:
                    sys.exit('Disease and target edge not found. Edge type must be \'biolink:DiseaseToPhenotypicFeatureAssociation\'')
                elif total_edges > 1:
                    sys.exit('Disease has too many outgoing edges')
                total_nodes += 1
        if total_nodes  == 1:
            acceptable_disease_curies_print = ','.join(acceptable_disease_curies)
            sys.exit('Disease node not found. Node type must be \'biolink:Disease\' and curie must be in: ' + acceptable_disease_curies_print)
        elif total_nodes > 2:
            sys.exit('Too many disease nodes')
        # set BKB target
        targets.append(('survival_time', qualifier, days))

        # get evidence
        evidence = dict()
        meta_evidence = list()
        for node_key in self.qg['nodes'].keys():
            # genes
            node = self.qg['nodes'][node_key]
            if node['category'] == 'biolink:Gene':
                # check for appropriate gene node structure
                gene_id = node_key
                for edge_key in self.qg['edges'].keys():
                    edge = self.qg['edges'][edge_key]
                    if edge['predicate'] == 'biolink:GeneToDiseaseAssociation' and edge['subject'] == gene_id and edge['object'] == disease_id:
                        total_edges += 1
                if total_edges == total_nodes - 1:
                    sys.exit('Gene and disease edge not found. Edge type must be \'biolink:GeneToDiseaseAssociation\'')
                elif total_edges > total_nodes:
                    sys.exit('Gene has too many outgoing edges')
                # check for appropriate gene node curie
                gene_curie = node['id']
                if gene_curie in self.gene_curie_dict.keys():
                    gene = gene_curie
                else:
                    sys.exit('Invalid ENSEMBL Identifier. Must be in form ENSEMBL:<ID>.')
                evidence["_mut_" + gene] = 'True'
                total_nodes += 1
            # drugs
            if node['category'] == 'biolink:Drug':
                # check for appropriate drug node structure
                drug_id = node_key
                for edge_key in self.qg['edges'].keys():
                    edge = self.qg['edges'][edge_key]
                    if edge['predicate'] == 'biolink:ChemicalToDiseaseOrPhenotypicFeatureAssociation' and edge['subject'] == drug_id and edge['object'] == disease_id:
                        total_edges += 1
                if total_edges == total_nodes - 1:
                    sys.exit('Drug and disease edge not found. Edge type must be \'biolink:ChemicalToDiseaseOrPhenotypicFeatureAssociation\'')
                elif total_edges > total_nodes:
                    sys.exit('Drug has too many outgoing edges')
                # check for appropriate drug node curie
                drug_curie = node['id']
                if drug_curie in self.drug_curie_dict.keys():
                    drug = drug_curie
                else:
                    sys.exit('Invalid CHEMBL Identifier. Must be in form CHEMBL:<ID>')
                meta_evidence.append(('drug_curies', '==', drug))
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

    def runQueries(self):
        """ Runs build BKB query to calculate probability of survival.
            A probability is returned to specificy survival time w.r.t evidence.
            Traditional bkb contributions are evaluated and will include contributions
            for all pieces in the evidence.
        """
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
        self.report = {'patient_analysis': report['Patient Analysis'],
                       'contribution_analysis': report['Contributions Analysis']}

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
        """ Knowledge graph is a copy of query graph where
            we add a edge property, 'has_confidence_level' annotated with
            our determined value. The whole query response is formed and
            returned. If specified in the QG, contributions will be stored
            under the KG edge between disease and survival nodes.

            :return: reasoner_std is our API response message combining KG and results.
            :rtype: dict
        """

        self.kg = copy.deepcopy(self.qg)
        # update target node info and form edge pair combos for results graph

        node_pairs = dict()
        for node_key in list(self.kg['nodes'].keys())[:]:
            qg_node_curie = self.kg['nodes'][node_key].pop('id')
            self.kg['nodes'][qg_node_curie] = self.kg['nodes'].pop(node_key)
            node_pairs[node_key] = qg_node_curie
            if self.kg['nodes'][qg_node_curie]['category'] == 'biolink:Gene':
                self.kg['nodes'][qg_node_curie]['name'] = self.gene_curie_dict[qg_node_curie]
            elif self.kg['nodes'][qg_node_curie]['category'] == 'biolink:Drug':
                self.kg['nodes'][qg_node_curie]['name'] = self.drug_curie_dict[qg_node_curie]

        edge_pairs = dict()
        knowledge_edges = 0
        for edge_key in list(self.kg['edges'].keys())[:]:
            kg_id = 'kge{}'.format(knowledge_edges)
            knowledge_edges += 1
            self.kg['edges'][kg_id] = self.kg['edges'].pop(edge_key)
            self.kg['edges'][kg_id]['subject'] = node_pairs[self.kg['edges'][kg_id]['subject']]
            self.kg['edges'][kg_id]['object'] = node_pairs[self.kg['edges'][kg_id]['object']]
            edge_pairs[edge_key] = kg_id
            if self.kg['edges'][kg_id]['predicate'] == 'biolink:DiseaseToPhenotypicFeatureAssociation':
                self.kg['edges'][kg_id]['has_confidence_level'] = self.truth_assignment
                if 'properties' in self.kg['edges'][kg_id].keys() and 'contributions' in self.kg['edges'][kg_id]['properties'].keys() and self.kg['edges'][kg_id]['properties']['contributions'] == True:
                    self.kg['edges'][kg_id]['properties'] = {'contributions':self.report}

        self.results.append({'edge_bindings':dict(),
                             'node_bindings':dict()})
        for edge_pair_key in edge_pairs:
            self.results[0]['edge_bindings'][edge_pair_key] = [{ 'id': edge_pairs[edge_pair_key]}]
        for node_pair_key in node_pairs:
            self.results[0]['node_bindings'][node_pair_key] = [{ 'id': node_pairs[node_pair_key]}]

        # query response
        reasoner_std = {'query_graph': self.qg,
                        'knowledge_graph': self.kg,
                        'results': self.results}
        reasoner_std = {'message' : reasoner_std}
        return reasoner_std
