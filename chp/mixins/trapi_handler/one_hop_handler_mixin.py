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

# Default functions
def get_default_predicate_proxy():
    return 'EFO:0000714'

def get_default_operator(predicate_proxy):
    if predicate_proxy == 'EFO:0000714':
        return '>='
    else:
        raise ValueError('Unknown predicate proxy: {}'.format(predicate_proxy))

def get_default_value(predicate_proxy):
    if predicate_proxy == 'EFO:0000714':
        return 978
    else:
        raise ValueError('Unknown predicate proxy: {}'.format(predicate_proxy))

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
        if self.messages is not None:
            # Setup queries
            self._setup_messages()

            # Instiatate Reasoners
            if self.dynamic_reasoner is None:
                self.dynamic_reasoner = ChpDynamicReasoner(
                    bkb_handler=self.bkb_data_handler,
                    hosts_filename=self.hosts_filename,
                    num_processes_per_host=self.num_processes_per_host)

    def _setup_messages(self):
        self.messages_dict = defaultdict(list)
        for message in self.messages:
            self.query_dict[self._get_onehop_type(message)].append(self._setup_single_message(message))

    def _get_onehop_type(self, message):
        wildcard_type = None
        for node_id, node in message.query_graph.nodes.items():
            if node.ids[0] is None:
                if wildcard_type is None:
                    wildcard_type = node.categories[0]
        # If standard onehop query
        if wildcard_type is None:
            return 'standard'
        elif wildcard_type == BiolinkEntity(BIOLINK_DRUG):
            return 'drug'
        elif wildcard_type == BiolinkEntity(BIOLINK_GENE):
            return 'gene'
        else:
            raise ValueError('Did not understand wildcard type {}.'.format(wildcard_type))

    def check_query(self):
        """ Currently not implemented. Would check validity of query.
        """
        return True

    @staticmethod
    def _process_predicate_proxy(qedge):
        dynamic_targets = {}
        predicate_proxy_constraint = qedge.find_constraint('CHP:PredicateProxy')
        if predicate_proxy_constraint is None:
            predicate_proxy = get_default_predicate_proxy()
            proxy_constraint = qedge.find_constraint(predicate_proxy)
        else:
            predicate_proxy = predicate_proxy_constraint.value
            proxy_constraint = None
        if proxy_constraint is None:
            proxy_operator = get_default_operator(predicate_proxy)
            proxy_value = get_default_value(predicate_proxy)
        else:
            proxy_operator = proxy_constraint.operator
            proxy_value = proxy_constraint.value
        # Setup dynamic target
        dynamic_targets[predicate_proxy] = {
                "op": proxy_operator,
                "value": proxy_value,
                }
        return dynamic_targets

    @staticmethod
    def _process_predicate_context(qedge):
        evidence = {}
        predicate_context_constraint = qedge.find_constraint('CHP:PredicateContext')
        if predicate_context_constraint is not None:
            for context in predicate_context_constraint.value:
                context_curie = BiolinkEntity(context)
                context_constraint = qedge.find_constraint(context)
                if context_constraint is None:
                    raise ValueError('Provided no context details for {}'.format(context))
                if context_curie == BiolinkEntity(BIOLINK_GENE) or context_curie == BiolinkEntity(BIOLINK_DRUG):
                    if type(context_constraint.value) is list:
                        for _curie in context_constraint.value:
                            evidence['_{}'.format(context_constraint.value)] = 'true'
                    else:
                        evidence['_{}'.format(context_constraint.value)] = 'true'
                else:
                    raise ValueError('Unsupported context type: {}'.format(context_curie))
        return evidence

    def _extract_chp_query(self, message, message_type):
        evidence = {}
        dynamic_targets = {}

        if message_type == 'standard':
            # Setup gene and drug evidence
            for qnode_id, qnode in message.query_graph.nodes.items():
                if qnode.categories[0] == BiolinkEntity(BIOLINK_GENE) or qnode.categories[0] == BiolinkEntity(BIOLINK_DRUG):
                    evidence['_{}'.format(qnode.ids[0])] = 'True'
        elif message_type == 'gene':
            for qnode_id, qnode in message.query_graph.nodes.items():
                if qnode.categories[0] == BiolinkEntity(BIOLINK_DRUG):
                    evidence['_{}'.format(qnode.ids[0])] = 'True'
        elif message_type == 'drug':
            for qnode_id, qnode in message.query_graph.nodes.items():
                if qnode.categories[0] == BiolinkEntity(BIOLINK_GENE):
                    evidence['_{}'.format(qnode.ids[0])] = 'True'
        # Grab edge
        for qedge_id, qedge in message.query_graph.edges.items():
            break
        # Process predicate proxy
        dynamic_targets = self._process_predicate_proxy(qedge)
        # Process predicate context
        evidence.update(self._process_predicate_context(qedge))

        target = list(dynamic_targets.keys())[0]
        truth_target = (target, '{} {}'.format(dynamic_targets[target]["op"], dynamic_targets[target]["value"]))

        chp_query = Query(evidence=evidence,
                      targets=None,
                      dynamic_evidence=None,
                      dynamic_targets=dynamic_targets,
                      type='updating')
        # Set some other helpful attributes
        chp_query.truth_target = truth_target
        return chp_query

    def _run_query(self, chp_query, query_type):
        """ Runs build BKB query to calculate probability of survival.
            A probability is returned to specificy survival time w.r.t a drug.
            Contributions for each gene are calculuated and classified under
            their true/false target assignments.
        """
        if query_type == 'gene':
            chp_query = self.dynamic_reasoner.run_query(chp_query, bkb_type='drug')
        elif query_type == 'drug' or query_type == 'standard':
            chp_query = self.dynamic_reasoner.run_query(chp_query, bkb_type='gene')
        chp_res_dict = chp_query.result.process_updates()
        chp_res_norm_dict = chp_query.result.process_updates(normalize=True)
        #chp_query.result.summary()
        chp_res_contributions = chp_query.result.process_inode_contributions()
        chp_query.truth_prob = max([0, chp_res_norm_dict[chp_query.truth_target[0]][chp_query.truth_target[1]]])

        #print(chp_res_contributions)
        if query_type != 'standard':
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

    def _construct_trapi_response(self, chp_query, message, query_type):
        
        qg = message.query_graph
        kg = message.knowledge_graph

        edge_bindings = {}
        node_bindings = {}
        
        # Process nodes
        for qnode_id, qnode in qg.nodes.items():
            if qnode.ids[0] is not None:
                if qnode.categories[0] == BiolinkEntity(BIOLINK_GENE):
                    knode_key = kg.add_node(
                            qnode.ids[0],
                            self.curies[BIOLINK_GENE][qnode.ids[0]][0],
                            qnode.categories[0].get_curie(),
                            )
                elif qnode.categories[0] == BiolinkEntity(BIOLINK_DRUG):
                    knode_key = kg.add_node(
                            qnode.ids[0],
                            self.curies[BIOLINK_DRUG][qnode.ids[0]][0],
                            qnode.categories[0].get_curie(),
                            )
                elif qnode.categories[0] == BiolinkEntity(BIOLINK_DISEASE):
                    knode_key = kg.add_node(
                            qnode.ids[0],
                            self.curies[BIOLINK_DISEASE][qnode.ids[0]][0],
                            qnode.categories[0].get_curie(),
                            )
                node_binding[qnode_id] = [kedge_key]
            else:
                wildcard_node = qnode
        if query_type == 'standard':
            for qedge_key, qedge in qg.edges.items():
                kedge_key = kg.add_edge(
                        node_bindings[qedge.subject][0],
                        node_bindings[qedge.object][0],
                        predicate=qedge.predicates[0].get_curie(),
                        relation=qedge.relation,
                        )
                edge_bindings[qedge_key] = [kedge_key]
                # Add Attribute
                kg.edges[kedge_key].add_attribute(
                        attribute_type_id='Probability of Survival',
                        value=chp_query.truth_prob,
                        value_type_id=BiolinkEntity(BIOLINK_PROBABILITY, is_slot=True).get_curie(),
                        )
            message.results.add_result(
                    node_bindings,
                    edge_bindings,
                    )
        else:
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
                _node_bindings = {}
                _edge_bindings = {}
                # Process node bindings
                for qnode_id, qnode in qg.nodes.items():
                    if qnode.categories[0] == BiolinkEntity(BIOLINK_GENE) and query_type == 'gene':
                        knode_id = kg.add_node(
                                wildcard,
                                self.curies[BIOLINK_GENE][wildcard][0],
                                qnode.categories[0].get_curie(),
                                )
                        _node_bindings[qnode_id] = [knode_id]
                    elif qnode.categories[0] == BiolinkEntity(BIOLINK_DRUG) and query_type == 'drug':
                        knode_id = kg.add_node(
                                wildcard,
                                self.curies[BIOLINK_DRUG][wildcard][0],
                                qnode.categories[0].get_curie(),
                                )
                        _node_bindings[qnode_id] = [knode_id]
                    else:
                        _node_bindings[qnode_id] = node_bindings[qnode_id]
                # Process edge bindings
                for qedge_id, qedge in qg.edges.items():
                    kedge_id = kg.add_edge(
                            _node_bindings[qedge.subject][0],
                            _node_bindings[qedge.object][0],
                            predicate=qedge.predicates[0],
                            relation=qedge.relation,
                            )
                    kg.edges[kedge_id].add_attribute(
                            attribute_type_id='Contribution',
                            value=contrib,
                            value_type_id=BiolinkEntity(BIOLINK_CONTRIBUTION, is_slot=True).get_curie(),
                            )
                    _edge_bindings[qedge_id] = [kedge_id]
                # Process node and edge binding results
                message.results.add_result(
                        _node_bindings,
                        _edge_bindings,
                        )
        return message
