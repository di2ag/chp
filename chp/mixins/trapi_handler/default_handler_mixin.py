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
import logging
import sys
import uuid
import json
from collections import defaultdict

from trapi_model.biolink.constants import *
from chp_data.bkb_handler import BkbDataHandler

from chp.query import Query as ChpQuery
from chp.reasoner import ChpDynamicReasoner, ChpJointReasoner

# Setup logging
logger = logging.getLogger(__name__)

class DefaultHandlerMixin:
    def _setup_handler(self):
        # Only do the rest of this if a message is passed
        if self.messages is not None:
            # Setup messages
            self._setup_messages()

            # Instiatate Reasoners
            if 'default' in self.message_dict:
                if self.dynamic_reasoner is None:
                    self.dynamic_reasoner = ChpDynamicReasoner(
                        bkb_handler=self.bkb_data_handler,
                        hosts_filename=self.hosts_filename,
                        num_processes_per_host=self.num_processes_per_host)
            if 'simple' in self.message_dict:
                if self.joint_reasoner is None:
                    self.joint_reasoner = ChpJointReasoner(
                        bkb_handler=self.bkb_data_handler,
                        hosts_filename=self.hosts_filename,
                        num_processes_per_host=self.num_processes_per_host)

    def _setup_messages(self):
        self.message_dict = defaultdict(list)
        for message in self.messages:
            if self._is_simple_message(message):
                self.message_dict['simple'].append(message)
            else:
                self.message_dict['default'].append(message)

    def _is_simple_message(self, message):
        """ Check if this is a {0 or 1} drug, {0 or 1} gene, one outcome standard message.
        """
        _found_outcome = False
        _found_disease = False
        _found_gene = False
        _found_drug = False
        query_graph = message.query_graph
        for node_key, node in query_graph.nodes.items():
            if node.categories[0] == BIOLINK_PHENOTYPIC_FEATURE_ENTITY:
                # If we've already found the target and there's another phenotypic feature, then this isn't simple.
                if _found_outcome:
                    return False
                else:
                    _found_outcome = True
            if node.categories[0] == BIOLINK_DISEASE_ENTITY:
                # If we've already found disease and there's another disease, then this isn't simple.
                if _found_disease:
                    return False
                else:
                    _found_disease = True
            if node.categories[0] == BIOLINK_GENE_ENTITY:
                if _found_gene:
                    return False
                else:
                    _found_gene = True
            if node.categories[0] == BIOLINK_DRUG_ENTITY:
                if _found_drug:
                    return False
                else:
                    _found_drug = True
        return True

    def _extract_chp_query(self, message, message_type=None):
        # Initialize Chp Query
        chp_query = ChpQuery(reasoning_type='updating')
        # Ensure we are using all nodes/edges
        total_nodes = 0
        total_edges = 0

        query_graph = message.query_graph

        # get phenotype node
        targets = list()
        for node_key in query_graph.nodes.keys():
            node = query_graph.nodes[node_key]
            if node.categories[0] == BIOLINK_PHENOTYPIC_FEATURE_ENTITY:
                target_id = node_key
                total_nodes += 1

        survival_value = 970
        survival_operator = '>='

        # get disease node info and ensure only 1 disease:
        for node_key in query_graph.nodes.keys():
            node = query_graph.nodes[node_key]
            if node.categories[0] == BIOLINK_DISEASE_ENTITY:
                disease_id = node_key
                for edge_key in query_graph.edges.keys():
                    edge = query_graph.edges[edge_key]
                    if self.check_predicate_support(edge.predicates[0], BIOLINK_HAS_PHENOTYPE_ENTITY) and edge.subject == disease_id and edge.object == target_id:
                        survival_time_constraint = edge.find_constraint(name='survival_time')
                        if survival_time_constraint is not None:
                            survival_value = survival_time_constraint.value
                            survival_operator = survival_time_constraint.operator
                            if survival_operator == 'matches':
                                survival_operator = '=='
                        total_edges += 1
                total_nodes += 1
        # set BKB target
        chp_query.add_dynamic_target(node.ids[0], survival_operator, survival_value)
        truth_target = (node.ids[0], '{} {}'.format(survival_operator, survival_value))

        # get evidence
        for node_key in query_graph.nodes.keys():
            # genes
            node = query_graph.nodes[node_key]
            if node.categories[0] == BIOLINK_GENE_ENTITY:
                # check for appropriate gene node structure
                gene_id = node_key
                for edge_key in query_graph.edges.keys():
                    edge = query_graph.edges[edge_key]
                    if self.check_predicate_support(edge.predicates[0], BIOLINK_GENE_ASSOCIATED_WITH_CONDITION_ENTITY) and edge.subject == gene_id and edge.object == disease_id:
                        total_edges += 1
                # check for appropriate gene node curie
                gene_curie = node.ids[0]
                gene = gene_curie
                chp_query.add_meta_evidence(gene, 'True')
                total_nodes += 1
            # drugs
            if node.categories[0] == BIOLINK_DRUG_ENTITY:
                # check for appropriate drug node structure
                drug_id = node_key
                for edge_key in query_graph.edges.keys():
                    edge = query_graph.edges[edge_key]
                    if self.check_predicate_support(edge.predicates[0], BIOLINK_TREATS_ENTITY) and edge.subject == drug_id and edge.object == disease_id:
                        total_edges += 1
                # check for appropriate drug node curie
                drug_curie = node.ids[0]
                drug = drug_curie
                chp_query.add_dynamic_evidence(node.ids[0], '==', 'True')
                total_nodes += 1

        # Set some other helpful attributes
        chp_query.truth_target = truth_target
        return chp_query

    def _run_query(self, chp_query, query_type):
        if query_type == 'simple':
            chp_query = self.joint_reasoner.run_query(chp_query)
            # If a probability was found for the target
            if len(chp_query.result) > 0:
                # If a probability was found for the truth target
                if chp_query.truth_target in chp_query.result:
                    total_unnormalized_prob = 0
                    for target, contrib in chp_query.result.items():
                        prob = max(0, contrib)
                        total_unnormalized_prob += prob
                    chp_query.truth_prob = max([0, chp_query.result[(chp_query.truth_target)]])/total_unnormalized_prob
                else:
                    chp_query.truth_prob = 0
            else:
                chp_query.truth_prob = -1
            chp_query.report = None
        else:
            chp_query = self.dynamic_reasoner.run_query(chp_query)
            chp_res_dict = chp_query.result.process_updates(normalize=True)
            try:
                chp_query.truth_prob = max([0, chp_res_dict[chp_query.truth_target[0]][chp_query.truth_target[1]]])
            except KeyError:
                # May need to come back and fix this.
                chp_query.truth_prob = -1

            chp_query.report = None
        return chp_query

    def _construct_trapi_message(self, chp_query, message, query_type=None):

        # update target node info and form edge pair combos for results graph

        qg = message.query_graph
        kg = message.knowledge_graph
        node_bindings = {}
        for qnode_key, qnode in qg.nodes.items():
            if qnode.categories[0] == BIOLINK_GENE_ENTITY:
                knode_key = kg.add_node(
                        qnode.ids[0],
                        self.curies[BIOLINK_GENE_ENTITY.get_curie()][qnode.ids[0]][0],
                        qnode.categories[0].get_curie(),
                        )
            elif qnode.categories[0] == BIOLINK_DRUG_ENTITY:
                knode_key = kg.add_node(
                        qnode.ids[0],
                        self.curies[BIOLINK_DRUG_ENTITY.get_curie()][qnode.ids[0]][0],
                        qnode.categories[0].get_curie(),
                        )
            else:
                knode_key = kg.add_node(
                        qnode.ids[0],
                        qnode.ids[0],
                        qnode.categories[0].get_curie(),
                        )
            node_bindings[qnode_key] = [knode_key]

        edge_bindings = {}
        for qedge_key, qedge in qg.edges.items():
            kedge_key = kg.add_edge(
                    node_bindings[qedge.subject][0],
                    node_bindings[qedge.object][0],
                    predicate=qedge.predicates[0].get_curie(),
                    relation=qedge.relation,
                    )
            edge_bindings[qedge_key] = [kedge_key]
            # Add Attribute
            if self.check_predicate_support(qedge.predicates[0], BIOLINK_HAS_PHENOTYPE_ENTITY):
                kg.edges[kedge_key].add_attribute(
                        attribute_type_id='Probability of Survival',
                        value=chp_query.truth_prob,
                        value_type_id=BIOLINK_HAS_CONFIDENCE_LEVEL_ENTITY.get_curie(),
                        )
        # Proces results
        message.results.add_result(
                node_bindings,
                edge_bindings,
                )
        return message
