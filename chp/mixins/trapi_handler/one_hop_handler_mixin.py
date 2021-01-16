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

from chp_data.trapi_constants import *

from chp.query import Query
from chp.reasoner import ChpDynamicReasoner
from chp_data.bkb_handler import BkbDataHandler

class OneHopHandlerMixin:
    """ OneHopeHandler is the handler for 1-hop queries. That is
        query graphs (QGs) that consists of 2 nodes and a single edge.

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
    def _setup_handler(self):
        self.default_survival_target = {
            "EFO:0000714": {
                "op": '>=',
                "value": 970
            }
        }

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
        if wildcard_type == BIOLINK_DRUG:
            return 'drug'
        elif wildcard_type == BIOLINK_GENE:
            return 'gene'
        else:
            raise ValueError('Did not understand wildcard type {}.'.format(wildcard_type))

    def check_query(self):
        """ Currently not implemented. Would check validity of query.
        """
        return True

    def _extract_chp_query(self, query, query_type=None):
        evidence = {}
        dynamic_targets = {}

        if len(query["query_graph"]['nodes']) > 2 or len(query["query_graph"]['edges']) > 1:
            sys.exit('1 hop quries can only have 2 nodes and 1 edge')

        # check edge for source and target
        edge_key = list(query["query_graph"]["edges"].keys())[0]
        edge = query["query_graph"]['edges'][edge_key]
        if 'subject' not in edge.keys() or 'object' not in edge.keys():
            sys.exit('Edge must have both a \'subject\' and and \'object\' key')
        subject = edge['subject']
        obj = edge['object']

        # Get non-wildcard node
        if query_type == 'gene':
            if query["query_graph"]['nodes'][subject]['category'] != BIOLINK_GENE:
                sys.exit('Subject node must be \'category\' {}'.format(BIOLINK_GENE))
            drug_curie = query["query_graph"]['nodes'][obj]['id']
            if drug_curie not in self.curies[BIOLINK_DRUG]:
                sys.exit('Invalid CHEMBL Identifier. Must be CHEMBL:<ID>')
            evidence['_{}'.format(drug_curie)] = 'True'
        elif query_type == 'drug':
            if query["query_graph"]['nodes'][subject]['category'] != BIOLINK_DRUG:
                sys.exit('Subject node must be \'category\' {}'.format(BIOLINK_DRUG))
            gene_curie = query["query_graph"]['nodes'][obj]['id']
            if gene_curie not in self.curies[BIOLINK_GENE]:
                sys.exit('Invalid ENSEMBL Identifier. Must be ENSEMBL:<ID>')
            evidence['_{}'.format(gene_curie)] = 'True'

        # default survival time
        dynamic_targets.update(self.default_survival_target)
        truth_target = ('EFO:0000714', '{} {}'.format(self.default_survival_target["EFO:0000714"]["op"],
                                                      self.default_survival_target["EFO:0000714"]["value"]))

        chp_query = Query(evidence=evidence,
                      targets=None,
                      dynamic_evidence=None,
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
        chp_res_norm_dict = chp_query.result.process_updates(normalize=True)
        #chp_query.result.summary()
        chp_res_contributions = chp_query.result.process_inode_contributions()
        chp_query.truth_prob = max([0, chp_res_norm_dict[chp_query.truth_target[0]][chp_query.truth_target[1]]])

        #print(chp_res_contributions)

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
                    hash_len = len(source_hashes)
                    # Process patient contributions
                    for _hash in source_hashes:
                        # Normalize to get relative contribution
                        patient_contributions[target][_hash] += contrib/hash_len #/ chp_res_dict[target_comp_name][target_state_name]

        # Now iterate through the patient data to translate patient contributions to drug/gene contributions
        wildcard_contributions = defaultdict(lambda: defaultdict(int))
        for target, patient_contrib_dict in patient_contributions.items():
            for patient, contrib in patient_contrib_dict.items():
                if query_type == 'gene':
                    for gene_curie in self.dynamic_reasoner.raw_patient_data[patient]["gene_curies"]:
                        wildcard_contributions[gene_curie][target] += contrib
                elif query_type == 'drug':
                    for drug_curie in self.dynamic_reasoner.raw_patient_data[patient]["drug_curies"]:
                        wildcard_contributions[drug_curie][target] += contrib

        # normalize gene contributions by the target and take relative difference
        for curie in wildcard_contributions.keys():
            truth_target_gene_contrib = 0
            nontruth_target_gene_contrib = 0
            for target, contrib in wildcard_contributions[curie].items():
                if target[0] == chp_query.truth_target[0] and target[1] == chp_query.truth_target[1]:
                    truth_target_gene_contrib += contrib / chp_res_dict[target[0]][target[1]]
                else:
                    nontruth_target_gene_contrib += contrib / chp_res_dict[target[0]][target[1]]
            wildcard_contributions[curie]['relative'] = truth_target_gene_contrib - nontruth_target_gene_contrib

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

        kg = copy.deepcopy(query["query_graph"])

        edge_bindings = {}
        node_bindings = {}

        # get edge subject, object, edge label and pop edge
        edge_key = list(kg['edges'].keys())[0]
        edge = kg['edges'][edge_key]
        edge_label = edge['predicate']
        subject = edge['subject']
        obj = edge['object']
        kg['edges'].pop(edge_key)

        # move curie to key
        non_wildcard_curie = kg['nodes'][obj].pop('id')
        kg['nodes'][non_wildcard_curie] = kg['nodes'].pop(obj)
        if query_type == 'gene':
            kg['nodes'][non_wildcard_curie]['name'] = self._get_curie_name(BIOLINK_DRUG, non_wildcard_curie)[0]
        elif query_type == 'drug':
            kg['nodes'][non_wildcard_curie]['name'] = self._get_curie_name(BIOLINK_GENE, non_wildcard_curie)[0]
        node_bindings[obj] = non_wildcard_curie

        # remove wildcard gene node from kg
        kg['nodes'].pop(subject)

        # Build relative contribution results and added associated edges into knowledge graph
        unsorted_wildcard_contributions = []
        for wildcard, contrib_dict in chp_query.wildcard_contributions.items():
            unsorted_wildcard_contributions.append((contrib_dict['relative'], wildcard))
        sorted_wildcard_contributions = [(contrib,wildcard) for contrib, wildcard in sorted(unsorted_wildcard_contributions, key=lambda x: abs(x[0]), reverse=True)]

        # add kg gene nodes and edges
        edge_count = 0
        node_count = 1
        results = []
        for contrib, wildcard in sorted_wildcard_contributions[:self.max_results]:
            if query_type == 'gene':
                kg['nodes'][wildcard] = {
                    "name" : self._get_curie_name(BIOLINK_GENE, wildcard)[0],
                    "category" : BIOLINK_GENE
                }
                # add edge
                kg['edges']['kge{}'.format(edge_count)] = {
                    "predicate" : BIOLINK_CHEMICAL_TO_GENE_PREDICATE,
                    "subject" : wildcard,
                    "object" : non_wildcard_curie,
                    "attributes":[{'name':'Contribution',
                                   'type':BIOLINK_CONTRIBUTION,
                                   'value':contrib}]
                }
            elif query_type == 'drug':
                kg['nodes'][wildcard] = {
                    "name" : self._get_curie_name(BIOLINK_DRUG, wildcard)[0],
                    "category" : BIOLINK_DRUG
                }
                # add edge
                kg['edges']['kge{}'.format(edge_count)] = {
                    "predicate" : BIOLINK_CHEMICAL_TO_GENE_PREDICATE,
                    "subject" : wildcard,
                    "object" : non_wildcard_curie,
                    "attributes":[{'name':'Contribution',
                                   'type':BIOLINK_CONTRIBUTION,
                                   'value':contrib}]
                }
            # add to results
            node_binding = {obj : [{'id': non_wildcard_curie}],
                            subject : [{'id':wildcard}]}
            edge_binding = {edge_key : [{'id':'kge{}'.format(edge_count)}]}
            results.append({'node_bindings': node_binding,
                            'edge_bindings': edge_binding})

            edge_count += 1
            node_count += 1

        # query response
        trapi_message = {'query_graph': query["query_graph"],
                        'knowledge_graph': kg,
                        'results': results}
        trapi_response = {'message' : trapi_message}
        return query_id, trapi_response
