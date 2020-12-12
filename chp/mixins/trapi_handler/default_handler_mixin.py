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
import json
from collections import defaultdict

from chp_data.bkb_handler import BkbDataHandler

from chp.query import Query
from chp.reasoner import ChpDynamicReasoner, ChpJointReasoner


class DefaultHandlerMixin:
    def _setup_handler(self):
        # Only do the rest of this if a query is passed
        if self.init_query is not None:
            # Setup queries
            self._setup_queries()

            # Instiatate Reasoners
            if 'default' in self.query_dict:
                if self.dynamic_reasoner is None:
                    self.dynamic_reasoner = ChpDynamicReasoner(
                        bkb_handler=self.bkb_data_handler,
                        hosts_filename=self.hosts_filename,
                        num_processes_per_host=self.num_processes_per_host)
            if 'simple' in self.query_dict:
                if self.joint_reasoner is None:
                    self.joint_reasoner = ChpJointReasoner(
                        bkb_handler=self.bkb_data_handler,
                        hosts_filename=self.hosts_filename,
                        num_processes_per_host=self.num_processes_per_host)

    def _setup_queries(self):
        if type(self.init_query) == list:
            self.query_dict = defaultdict(list)
            self.query_map = []
            for query in self.init_query:
                self.query_map.append(query["query_id"])
                if self._is_simple_query(query):
                    self.query_dict['simple'].append(self._setup_single_query(query))
                else:
                    self.query_dict['default'].append(self._setup_single_query(query))
        else:
            if self._is_simple_query(self.init_query):
                self.query_dict = {"simple": [self._setup_single_query(self.init_query)]}
            else:
                self.query_dict = {"default": [self._setup_single_query(self.init_query)]}

    def _is_simple_query(self, query):
        """ Check if this is a {0 or 1} drug, {0 or 1} gene, one outcome standard query.
        """
        _found_outcome = False
        _found_disease = False
        _found_gene = False
        _found_drug = False
        for node_key, node in query["query_graph"]["nodes"].items():
            if node["category"] == 'biolink:PhenotypicFeature':
                # If we've already found the target and there's another phenotypic feature, then this isn't simple.
                if _found_outcome:
                    return False
                else:
                    _found_outcome = True
            if node['category'] == 'biolink:Disease':
                # If we've already found disease and there's another disease, then this isn't simple.
                if _found_disease:
                    return False
                else:
                    _found_disease = True
            if node["category"] == 'biolink:Gene':
                if _found_gene:
                    return False
                else:
                    _found_gene = True
            if node['category'] == 'biolink:Drug':
                if _found_drug:
                    return False
                else:
                    _found_drug = True
        return True

    def _extract_chp_query(self, query, query_type=None):
        evidence = {}
        targets = []
        dynamic_evidence = {}
        dynamic_targets = {}
        # ensure we are using all nodes/edges
        total_nodes = 0
        total_edges = 0

        # get phenotype node
        targets = list()
        acceptable_target_curies = ['EFO:0000714']
        for node_key in query["query_graph"]['nodes'].keys():
            node = query["query_graph"]['nodes'][node_key]
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
        for node_key in query["query_graph"]['nodes'].keys():
            node = query["query_graph"]['nodes'][node_key]
            if node['category'] == 'biolink:Disease' and node['id'] in acceptable_disease_curies:
                disease_id = node_key
                for edge_key in query["query_graph"]['edges'].keys():
                    edge = query["query_graph"]['edges'][edge_key]
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
        dynamic_targets[node["id"]] = {
            "op": qualifier,
            "value": days,
        }
        truth_target = (node["id"], '{} {}'.format(qualifier, days))

        # get evidence
        for node_key in query["query_graph"]['nodes'].keys():
            # genes
            node = query["query_graph"]['nodes'][node_key]
            if node['category'] == 'biolink:Gene':
                # check for appropriate gene node structure
                gene_id = node_key
                for edge_key in query["query_graph"]['edges'].keys():
                    edge = query["query_graph"]['edges'][edge_key]
                    if edge['predicate'] == 'biolink:GeneToDiseaseAssociation' and edge['subject'] == gene_id and edge['object'] == disease_id:
                        total_edges += 1
                if total_edges == total_nodes - 1:
                    sys.exit('Gene and disease edge not found. Edge type must be \'biolink:GeneToDiseaseAssociation\'')
                elif total_edges > total_nodes:
                    sys.exit('Gene has too many outgoing edges')
                # check for appropriate gene node curie
                gene_curie = node['id']
                if gene_curie in self.curies["gene"]:
                    gene = gene_curie
                else:
                    sys.exit('Invalid ENSEMBL Identifier. Must be in form ENSEMBL:<ID>.')
                evidence["_" + gene] = 'True'
                total_nodes += 1
            # drugs
            if node['category'] == 'biolink:Drug':
                # check for appropriate drug node structure
                drug_id = node_key
                for edge_key in query["query_graph"]['edges'].keys():
                    edge = query["query_graph"]['edges'][edge_key]
                    if edge['predicate'] == 'biolink:ChemicalToDiseaseOrPhenotypicFeatureAssociation' and edge['subject'] == drug_id and edge['object'] == disease_id:
                        total_edges += 1
                if total_edges == total_nodes - 1:
                    sys.exit('Drug and disease edge not found. Edge type must be \'biolink:ChemicalToDiseaseOrPhenotypicFeatureAssociation\'')
                elif total_edges > total_nodes:
                    sys.exit('Drug has too many outgoing edges')
                # check for appropriate drug node curie
                drug_curie = node['id']
                if drug_curie in self.curies["chemical_substance"]:
                    drug = drug_curie
                else:
                    sys.exit('Invalid CHEMBL Identifier: {}. Must be in form CHEMBL:<ID>'.format(drug_curie))
                evidence[node["id"]] = 'True'
                total_nodes += 1

        if total_nodes != len(query["query_graph"]['nodes']) or total_edges != len(query["query_graph"]['edges']):
            sys.exit('There are extra components in the provided QG structure')

        # produce BKB query
        chp_query = Query(
            evidence=evidence,
            targets=targets,
            dynamic_evidence=dynamic_evidence,
            dynamic_targets=dynamic_targets,
            type='updating')
        # Set some other helpful attributes
        chp_query.truth_target = truth_target
        chp_query.query_id = query["query_id"] if 'query_id' in query else None
        return chp_query

    def _run_query(self, chp_query, query_type):
        if query_type == 'simple':
            chp_query = self.joint_reasoner.run_query(chp_query)
            # If a probability was found for the target
            if len(chp_query.result) > 0:
                # If a probability was found for the truth target
                if chp_query.truth_target in chp_query.result:
                    chp_query.truth_prob = max([0, chp_query.result[(chp_query.truth_target)]])
                else:
                    chp_query.truth_prob = 0
            chp_query.report = None
        else:
            chp_query = self.dynamic_reasoner.run_query(chp_query)
            chp_res_dict = chp_query.result.process_updates()
            chp_query.truth_prob = max([0, chp_res_dict[chp_query.truth_target[0]][chp_query.truth_target[1]]])
            chp_query.report = None
        return chp_query

    def _construct_trapi_response(self, chp_query, query_type=None):
        # Get orginal query
        if len(self.init_query) == 1:
            query = self.init_query[0]
            query_id = None
        else:
            for _query in self.init_query:
                if _query["query_id"] == chp_query.query_id:
                    query = _query
                    query_id = query["query_id"]
                    break

        kg = copy.deepcopy(query["query_graph"])
        # update target node info and form edge pair combos for results graph

        node_pairs = dict()
        for node_key in list(kg['nodes'].keys())[:]:
            qg_node_curie = kg['nodes'][node_key].pop('id')
            kg['nodes'][qg_node_curie] = kg['nodes'].pop(node_key)
            node_pairs[node_key] = qg_node_curie
            if kg['nodes'][qg_node_curie]['category'] == 'biolink:Gene':
                kg['nodes'][qg_node_curie]['name'] = self._get_curie_name('gene', qg_node_curie)
            elif kg['nodes'][qg_node_curie]['category'] == 'biolink:Drug':
                kg['nodes'][qg_node_curie]['name'] = self._get_curie_name('chemical_substance', qg_node_curie)

        edge_pairs = dict()
        knowledge_edges = 0
        for edge_key in list(kg['edges'].keys())[:]:
            kg_id = 'kge{}'.format(knowledge_edges)
            knowledge_edges += 1
            kg['edges'][kg_id] = kg['edges'].pop(edge_key)
            kg['edges'][kg_id]['subject'] = node_pairs[kg['edges'][kg_id]['subject']]
            kg['edges'][kg_id]['object'] = node_pairs[kg['edges'][kg_id]['object']]
            edge_pairs[edge_key] = kg_id
            if kg['edges'][kg_id]['predicate'] == 'biolink:DiseaseToPhenotypicFeatureAssociation':
                kg['edges'][kg_id]['has_confidence_level'] = chp_query.truth_prob
                if 'properties' in kg['edges'][kg_id].keys() and 'contributions' in kg['edges'][kg_id]['properties'].keys() and kg['edges'][kg_id]['properties']['contributions'] == True:
                    kg['edges'][kg_id]['properties'] = {'contributions':chp_query.report}

        results = []
        results.append({
            'edge_bindings': {},
            'node_bindings': {},
        })
        for edge_pair_key in edge_pairs:
            results[0]['edge_bindings'][edge_pair_key] = [{ 'id': edge_pairs[edge_pair_key]}]
        for node_pair_key in node_pairs:
            results[0]['node_bindings'][node_pair_key] = [{ 'id': node_pairs[node_pair_key]}]

        # query response
        trapi_message = {'query_graph': query["query_graph"],
                        'knowledge_graph': kg,
                        'results': results}
        trapi_response = {'message' : trapi_message}
        return query_id, trapi_response
