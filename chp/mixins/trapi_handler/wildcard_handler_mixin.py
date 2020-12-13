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
import pickle
from collections import defaultdict
import json

from chp_data.bkb_handler import BkbDataHandler

from chp.query import Query
from chp.reasoner import ChpDynamicReasoner


class WildCardHandlerMixin:
    def _setup_handler(self):
        # Only do the rest of this if a query is passed
        if self.init_query is not None:
            # Setup queries
            self._setup_queries()

            # Instiatate Reasoners
            if self.dynamic_reasoner is None:
                self.dynamic_reasoner = ChpDynamicReasoner(
                    bkb_handler=self.bkb_data_handler,
                    hosts_filename=self.hosts_filename,
                    num_processes_per_host=self.num_processes_per_host)

    def _setup_queries(self):
        if type(self.init_query) == list:
            self.query_dict = defaultdict(list)
            self.query_map = []
            for query in self.init_query:
                self.query_map.append(query["query_id"])
                self.query_dict[self._get_wildcard_type(query)].append(self._setup_single_query(query))
        else:
            self.query_dict[self._get_wildcard_type(query)].append(self._setup_single_query(query))

    def _get_wildcard_type(self, query):
        wildcard_type = None
        for node_id, node in query["query_graph"]["nodes"].items():
            if 'id' not in node:
                if wildcard_type is None:
                    wildcard_type = node['category']
                else:
                    sys.exit('You can only have one contribution target. Make sure to leave only one node with a black curie.')
        if wildcard_type == 'biolink:Drug':
            return 'drug'
        elif wildcard_type == 'biolink:Gene':
            return 'gene'
        else:
            raise ValueError('Did not understand wildcard type {}.'.format(wildcard_type))

    def _extract_chp_query(self, query, query_type):
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
                if query_type != 'gene':
                    gene_curie = node['id']
                    if gene_curie in self.curies["biolink:Gene"]:
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
                if query_type != 'drug':
                    drug_curie = node['id']
                    if drug_curie in self.curies["biolink:Drug"]:
                        drug = drug_curie
                    else:
                        sys.exit('Invalid CHEMBL Identifier: {}. Must be in form CHEMBL:<ID>'.format(drug_curie))
                    evidence['_' + drug] = 'True'
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
        """ Runs build BKB query to calculate probability of survival.
            A probability is returned to specificy survival time w.r.t a drug.
            Contributions for each gene are calculuated and classified under
            their true/false target assignments.
        """
        if query_type == 'gene':
            chp_query = self.dynamic_reasoner.run_query(chp_query, bkb_type='drug')
        elif query_type == 'drug':
            chp_query = self.dynamic_reasoner.run_query(chp_query, bkb_type='gene')
        chp_res_dict = chp_query.result.process_updates()
        #chp_query.result.summary()
        chp_res_contributions = chp_query.result.process_inode_contributions()
        chp_query.truth_prob = max([0, chp_res_dict[chp_query.truth_target[0]][chp_query.truth_target[1]]])

        # Collect all source inodes and process patient hashes
        patient_contributions = defaultdict(lambda: defaultdict(int))
        for target, contrib_dict in chp_res_contributions.items():
            target_comp_name, target_state_name = target
            for inode, contrib in contrib_dict.items():
                comp_name, state_name = inode
                if '_Source_' in comp_name:
                    # Split source state name to get patient hashes
                    source_hashes_str = state_name.split('_')[-1]
                    source_hashes = [int(source_hash) for source_hash in source_hashes_str.split(',')]
                    # Process patient contributions
                    for _hash in source_hashes:
                        # Normalize to get relative contribution
                        patient_contributions[target][_hash] += contrib / chp_res_dict[target_comp_name][target_state_name]

        # Now iterate through the patient data to translate patient contributions to drug/gene contributions
        wildcard_contributions = defaultdict(lambda: defaultdict(int))
        for target, patient_contrib_dict in patient_contributions.items():
            for patient, contrib in patient_contrib_dict.items():
                if query_type == 'gene':
                    for gene_curie in self.dynamic_reasoner.raw_patient_data[patient]["gene_curies"]:
                        wildcard_contributions[target][gene_curie] += contrib
                elif query_type == 'drug':
                    for drug_curie in self.dynamic_reasoner.raw_patient_data[patient]["drug_curies"]:
                        wildcard_contributions[target][drug_curie] += contrib

        chp_query.report = None
        chp_query.wildcard_contributions = wildcard_contributions

        return chp_query

    def _construct_trapi_response(self, chp_query, query_type):
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

        # Construct first result which is the result of the standard probablistic query.
        kg = copy.deepcopy(query["query_graph"])
        # Process Nodes
        node_pairs = defaultdict(None)
        contrib_qg_id = None
        for node_key in list(kg["nodes"].keys())[:]:
            qg_node_curie = kg['nodes'][node_key].pop('id', None)
            if qg_node_curie is not None:
                kg['nodes'][qg_node_curie] = kg['nodes'].pop(node_key)
                if kg['nodes'][qg_node_curie]['category'] == 'biolink:Gene':
                    kg['nodes'][qg_node_curie]['name'] = self.curies["biolink:Gene"][qg_node_curie]
                elif kg['nodes'][qg_node_curie]['category'] == 'biolink:Drug':
                    kg['nodes'][qg_node_curie]['name'] = self.curies["biolink:Drug"][qg_node_curie]
                node_pairs[node_key] = qg_node_curie
            else:
                kg["nodes"].pop(node_key)

        # Process Edges
        edge_pairs = dict()
        knowledge_edges = 0
        for edge_key in list(kg['edges'].keys())[:]:
            if kg['edges'][edge_key]['predicate'] == 'biolink:GeneToDiseaseAssociation' and query_type == 'gene':
                kg['edges'].pop(edge_key)
            elif kg['edges'][edge_key]['predicate'] == 'biolink:ChemicalToDiseaseOrPhenotypicFeatureAssociation' and query_type == 'drug':
                kg['edges'].pop(edge_key)
            else:
                kg_id = 'kge{}'.format(knowledge_edges)
                knowledge_edges += 1
                kg['edges'][kg_id] = kg['edges'].pop(edge_key)
                kg['edges'][kg_id]['subject'] = node_pairs[kg['edges'][kg_id]['subject']]
                kg['edges'][kg_id]['object'] = node_pairs[kg['edges'][kg_id]['object']]
                edge_pairs[edge_key] = kg_id
                if kg['edges'][kg_id]['predicate'] == 'biolink:DiseaseToPhenotypicFeatureAssociation':
                    kg['edges'][kg_id]['has_confidence_level'] = chp_query.truth_prob
                    if 'properties' in kg['edges'][kg_id].keys() and 'contributions' in kg['edges'][kg_id]['properties'].keys() and kg['edges'][kg_id]['properties']['contributions'] == True:
                        kg['edges'][kg_id]['Description'] = chp_query.report

        # Put first result of standard prob query of only curie nodes (i.e. no wildcard nodes where used as evidence)
        results = []
        results.append({'edge_bindings':dict(),
                        'node_bindings':dict()})
        for edge_pair_key in edge_pairs:
            results[0]['edge_bindings'][edge_pair_key] = [{ 'id': str(edge_pairs[edge_pair_key])}]
        for node_pair_key in node_pairs:
            results[0]['node_bindings'][node_pair_key] = [{ 'id': str(node_pairs[node_pair_key])}]

        # Build relative contribution results and added associated edges into knowledge graph
        unsorted_wildcard_contributions = []
        for wildcard, contrib in chp_query.wildcard_contributions[chp_query.truth_target].items():
            unsorted_wildcard_contributions.append((contrib, wildcard))
        sorted_wildcard_contributions = sorted(unsorted_wildcard_contributions, reverse=True)
        for contrib, wildcard in sorted_wildcard_contributions[:self.max_results]:
            rg = copy.deepcopy(query["query_graph"])
            _node_pairs = {}
            _edge_pairs = {}
            # Process node pairs
            for node_id, node in rg["nodes"].items():
                if node["category"] == 'biolink:Gene' and query_type == 'gene':
                    kg["nodes"][wildcard] = copy.deepcopy(node)
                    kg["nodes"][wildcard].update({"name": self.curies["biolink:Gene"][wildcard]})
                    _node_pairs[node_id] = wildcard
                elif node["category"] == 'biolink:Drug' and query_type == 'drug':
                    kg["nodes"][wildcard] = copy.deepcopy(node)
                    kg["nodes"][wildcard].update({"name": self.curies["biolink:Drug"][wildcard]})
                    _node_pairs[node_id] = wildcard
                else:
                    _node_pairs[node_id] = node_pairs[node_id]
            # Process edge pairs
            for edge_id, edge in rg["edges"].items():
                if query_type == 'gene'  and edge["predicate"] == 'biolink:GeneToDiseaseAssociation':
                    knowledge_edges += 1
                    kg_edge_id = 'kge{}'.format(knowledge_edges)
                    kg["edges"][kg_edge_id] = copy.deepcopy(edge)
                    kg["edges"][kg_edge_id]["subject"] = _node_pairs[kg["edges"][kg_edge_id]["subject"]]
                    kg["edges"][kg_edge_id]["object"] = _node_pairs[kg["edges"][kg_edge_id]["object"]]
                    kg["edges"][kg_edge_id]["value"] = contrib
                    _edge_pairs[edge_id] = kg_edge_id
                elif query_type == 'drug'  and edge["predicate"] == 'biolink:ChemicalToDiseaseOrPhenotypicFeatureAssociation':
                    knowledge_edges += 1
                    kg_edge_id = 'kge{}'.format(knowledge_edges)
                    kg["edges"][kg_edge_id] = copy.deepcopy(edge)
                    kg["edges"][kg_edge_id]["subject"] = _node_pairs[kg["edges"][kg_edge_id]["subject"]]
                    kg["edges"][kg_edge_id]["object"] = _node_pairs[kg["edges"][kg_edge_id]["object"]]
                    kg["edges"][kg_edge_id]["value"] = contrib
                    _edge_pairs[edge_id] = kg_edge_id
                else:
                    _edge_pairs[edge_id] = edge_pairs[edge_id]
            # Process node and edge binding results
            _res = {"edge_bindings": {},
                    "node_bindings": {}}
            for edge_pair_key in _edge_pairs:
                _res["edge_bindings"][edge_pair_key] = [{ "id": str(_edge_pairs[edge_pair_key])}]
            for node_pair_key in _node_pairs:
                _res["node_bindings"][node_pair_key] = [{ "id": str(_node_pairs[node_pair_key])}]
            results.append(_res)

        # query response
        trapi_message = {'query_graph': query["query_graph"],
                        'knowledge_graph': kg,
                        'results': results}
        trapi_response = {'message' : trapi_message}
        return query_id, trapi_response
