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

from trapi_model.biolink.constants import *

from chp.query import Query as ChpQuery
from chp.reasoner import ChpDynamicReasoner, ChpJointReasoner
from chp_data.bkb_handler import BkbDataHandler
from pybkb.python_base.utils import get_operator, get_opposite_operator

# Setup logging
logger = logging.getLogger(__name__)

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

def get_default_two_hop_proxy(message_type):
    if message_type == 'gene_two_hop':
        return BIOLINK_DRUG_ENTITY
    elif message_type == 'drug_two_hop':
        return BIOLINK_GENE_ENTITY
    else:
        return None

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
        if self.queries is not None:
            # Setup queries
            self._setup_messages()

            # Instiatate Reasoners
            if self.dynamic_reasoner is None:
                self.dynamic_reasoner = ChpDynamicReasoner(
                    bkb_handler=self.bkb_data_handler,
                    hosts_filename=self.hosts_filename,
                    num_processes_per_host=self.num_processes_per_host)
            if self.joint_reasoner is None:
                self.joint_reasoner = ChpJointReasoner(
                    bkb_handler=self.bkb_data_handler,
                    hosts_filename=self.hosts_filename,
                    num_processes_per_host=self.num_processes_per_host)

    def _setup_messages(self):
        self.queries_dict = defaultdict(list)
        for query in self.queries:
            self.queries_dict[self._get_onehop_type(query.message)].append(query)

    def _get_onehop_type(self, message):
        wildcard_type = None
        node_types = []
        all_node_categories = []
        for node_id, node in message.query_graph.nodes.items():
            if node.ids is None:
                if wildcard_type is None:
                    wildcard_type = node.categories[0]
                node_types.append(node.categories[0])
            all_node_categories.append(node.categories[0])

        # implicit 2-hop-queries
        if all(category == BIOLINK_GENE_ENTITY for category in all_node_categories):
            return 'gene_two_hop'
        elif all(category == BIOLINK_DRUG_ENTITY for category in all_node_categories):
            return 'drug_two_hop'

        # If standard onehop query
        if wildcard_type is None:
            return 'standard'
        elif wildcard_type == BIOLINK_DRUG_ENTITY:
            return 'drug'
        elif wildcard_type == BIOLINK_GENE_ENTITY:
            return 'gene'
        else:
            raise ValueError('Did not understand wildcard type {}.'.format(wildcard_type))

    def check_query(self):
        """ Currently not implemented. Would check validity of query.
        """
        return True

    @staticmethod
    def _process_predicate_proxy(qedge, chp_query):
        predicate_proxy_constraint = qedge.find_constraint('predicate_proxy')
        if predicate_proxy_constraint is None:
            predicate_proxy = get_default_predicate_proxy()
            proxy_constraint = qedge.find_constraint(predicate_proxy)
        else:
            predicate_proxy = predicate_proxy_constraint.value[0]
            proxy_constraint = qedge.find_constraint(predicate_proxy)
        if proxy_constraint is None:
            proxy_operator = get_default_operator(predicate_proxy)
            proxy_value = get_default_value(predicate_proxy)
        else:
            proxy_operator = proxy_constraint.operator
            proxy_value = proxy_constraint.value
        # Setup dynamic target
        chp_query.add_dynamic_target(predicate_proxy, proxy_operator, proxy_value)
        return chp_query


    @staticmethod
    def _process_predicate_context(qedge, message_type, chp_query):
        evidence = {}
        dynamic_evidence = {}
        predicate_context_constraint = qedge.find_constraint('predicate_context')

        if predicate_context_constraint is not None:
            for context in predicate_context_constraint.value:
                context_curie = get_biolink_entity(context)
                context_constraint = qedge.find_constraint(context)
                # used 2 hop structure where context curie is the proxy
                if context_constraint is None:
                     continue
                if context_curie == BIOLINK_GENE_ENTITY:
                    if message_type == 'gene' or message_type == 'drug_two_hop':
                            if type(context_constraint.value) is list:
                                for _curie in context_constraint.value:
                                    chp_query.add_dynamic_evidence(_curie, '==', 'True')
                            else:
                                chp_query.add_dynamic_evidence(context_constraint.value, '==', 'True')
                    else:
                        if type(context_constraint.value) is list:
                            for _curie in context_constraint.value:
                                chp_query.add_meta_evidence(_curie, 'True')
                        else:
                            chp.add_meta_evidence(_curie, 'True')
                elif context_curie == BIOLINK_DRUG_ENTITY:
                    if message_type == 'drug' or message_type == 'gene_two_hop':
                        if type(context_constraint.value) is list:
                            for _curie in context_constraint.value:
                                chp_query.add_dynamic_evidence(_curie, '==', 'True')
                        else:
                            chp_query.add_dynamic_evidence(context_constraint.value, '==', 'True')
                    else:
                        if type(context_constraint.value) is list:
                            for _curie in context_constraint.value:
                                chp_query.add_meta_evidence(_curie, 'True')
                        else:
                            chp_query.add_meta_evidence(_curie, 'True')
                else:
                    raise ValueError('Unsupported context type: {}'.format(context_curie))
        return chp_query

    def _extract_chp_query(self, query, message_type):
        # Extract Message
        message = query.message

        # Initialize CHP BKB Query
        chp_query = ChpQuery(reasoning_type='updating')

        # Grab edge
        for qedge_id, qedge in message.query_graph.edges.items():
            break
        # Process predicate proxy
        chp_query = self._process_predicate_proxy(qedge, chp_query)
        # Process predicate context
        chp_query = self._process_predicate_context(qedge, message_type, chp_query)
        #TODO: Probably need a more robust solution for when no context is provided in wildcard queries and you need it.
        #if len(evidence) == 0:
        #    raise ValueError('Did not supply context with a query that required context.')

        if message_type == 'standard':
            # Setup gene and drug evidence
            for qnode_id, qnode in message.query_graph.nodes.items():
                if qnode.categories[0] == BIOLINK_GENE_ENTITY or qnode.categories[0] == BIOLINK_DRUG_ENTITY:
                    chp_query.add_meta_evidence(qnode.ids[0], 'True')
        elif message_type == 'gene' or message_type == 'drug_two_hop':
            for qnode_id, qnode in message.query_graph.nodes.items():
                if qnode.categories[0] == BIOLINK_DRUG_ENTITY:
                    if qnode.ids is not None:
                        chp_query.add_meta_evidence(qnode.ids[0], 'True')
        elif message_type == 'drug' or message_type == 'gene_two_hop':
            for qnode_id, qnode in message.query_graph.nodes.items():
                if qnode.categories[0] == BIOLINK_GENE_ENTITY:
                    if qnode.ids is not None:
                        chp_query.add_meta_evidence(qnode.ids[0], 'True')

        target = list(chp_query.dynamic_targets.keys())[0]
        truth_target = (target, '{} {}'.format(chp_query.dynamic_targets[target]["op"], chp_query.dynamic_targets[target]["value"]))
        # Set some other helpful attributes
        chp_query.truth_target = truth_target
        return chp_query

    def _run_query(self, chp_query, query_type):
        """ Runs build BKB query to calculate probability of survival.
            A probability is returned to specificy survival time w.r.t a drug.
            Contributions for each gene are calculuated and classified under
            their true/false target assignments.
        """
        if query_type == 'standard':
            chp_query = self.joint_reasoner.run_query(chp_query, interpolation_type='gene')
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
            return chp_query
        else:
            # Do this if a disease node is present
            if query_type == 'gene' or query_type == 'drug_two_hop':
                chp_query = self.joint_reasoner.run_query(chp_query, interpolation_type='drug', contribution_type='gene')
            elif query_type == 'drug' or query_type == 'gene_two_hop':
                chp_query = self.joint_reasoner.run_query(chp_query, interpolation_type='gene', contribution_type='drug')

            chp_res_dict = chp_query.result
            if chp_query.truth_target in chp_res_dict:
                unnormalized_truth_prob = chp_res_dict[chp_query.truth_target]
            else:
                unnormalized_truth_prob = 0

            normalize = 0
            for target, prob in chp_res_dict.items():
                normalize += prob
            unnormalized_false_prob = normalize - unnormalized_truth_prob

            for target in chp_res_dict.keys():
                chp_res_dict[target] /= normalize

            if chp_query.truth_target in chp_res_dict:
                chp_query.truth_prob = chp_res_dict[chp_query.truth_target]
            else:
               chp_query.truth_target = 0

        # organize the contributions over curie then target
        wildcard_contributions = defaultdict(lambda: defaultdict(int))
        for target, curies in chp_query.contributions.items():
            for curie, contrib in curies.items():
                if curie[1] == 'True':
                    wildcard_contributions[curie[0]][target] = contrib

        #for curie in wildcard_contributions.keys():
        # normalize gene contributions by the target and take relative difference
        for curie in wildcard_contributions.keys():
            truth_target_gene_contrib = 0
            nontruth_target_gene_contrib = 0
            for target, contrib in wildcard_contributions[curie].items():
                try:
                    if target[0] == chp_query.truth_target[0] and target[1] == chp_query.truth_target[1]:
                        truth_target_gene_contrib += contrib / unnormalized_truth_prob
                    else:
                        nontruth_target_gene_contrib += contrib / unnormalized_false_prob
                except ZeroDivisionError:
                    continue
            wildcard_contributions[curie]['relative'] = truth_target_gene_contrib - nontruth_target_gene_contrib

        if query_type == 'drug_two_hop' or query_type == 'gene_two_hop':
            # Build relative contribution results and added associated edges into knowledge graph
            unsorted_wildcard_contributions = []
            for wildcard, contrib_dict in wildcard_contributions.items():
                unsorted_wildcard_contributions.append((contrib_dict['relative'], wildcard))
            truncated_sorted_wildcard_contributions = [(contrib,wildcard) for contrib, wildcard in sorted(unsorted_wildcard_contributions, key=lambda x: abs(x[0]), reverse=True)][:self.max_results]
            truncated_contribution_list = [curie[1] for curie in truncated_sorted_wildcard_contributions]

            chp_query.contributions = None
            wildcard_contributions = defaultdict(lambda: defaultdict(int))
            for contrib in truncated_contribution_list:
                chp_query_extended = copy.deepcopy(chp_query)
                chp_query_extended.add_meta_evidence(contrib, 'True')
                if query_type == 'drug_two_hop':
                    chp_query_extended = self.joint_reasoner.run_query(chp_query_extended, contribution_type='drug')
                else:
                    chp_query_extended = self.joint_reasoner.run_query(chp_query_extended, contribution_type='gene')

                chp_res_dict = chp_query_extended.result
                if chp_query_extended.truth_target in chp_res_dict:
                    extended_unnormalized_truth_prob = chp_res_dict[chp_query_extended.truth_target]
                else:
                    extended_unnormalized_truth_prob = 0
                normalize = 0
                for target, prob in chp_res_dict.items():
                    normalize += prob
                extended_unnormalized_false_prob = normalize - extended_unnormalized_truth_prob

                # organize the contributions over curie then target
                extended_wildcard_contributions = defaultdict(lambda: defaultdict(int))
                for target, curies in chp_query_extended.contributions.items():
                    for curie, contrib in curies.items():
                        if curie[1] == 'True':
                            extended_wildcard_contributions[curie[0]][target] = contrib

                # normalize gene contributions by the target and take relative difference
                for curie in extended_wildcard_contributions.keys():
                    truth_target_gene_contrib = 0
                    nontruth_target_gene_contrib = 0
                    for target, contrib in extended_wildcard_contributions[curie].items():
                        try:
                            if target[0] == chp_query_extended.truth_target[0] and target[1] == chp_query_extended.truth_target[1]:
                                truth_target_gene_contrib += contrib / extended_unnormalized_truth_prob * unnormalized_truth_prob
                            else:
                                nontruth_target_gene_contrib += contrib / extended_unnormalized_false_prob * unnormalized_false_prob
                        except ZeroDivisionError:
                            continue
                    wildcard_contributions[curie]['relative'] += (truth_target_gene_contrib - nontruth_target_gene_contrib)

        chp_query.report = None
        chp_query.wildcard_contributions = wildcard_contributions

        return chp_query

    def _construct_trapi_message(self, chp_query, query, query_type):

        # Helpful short cuts
        message = query.message
        qg = message.query_graph
        kg = message.knowledge_graph

        edge_bindings = {}
        node_bindings = {}

        # Process nodes
        for qnode_id, qnode in qg.nodes.items():
            if qnode.ids is not None:
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
                elif qnode.categories[0] == BIOLINK_DISEASE_ENTITY:
                    #TODO: Add diseases to curies and fix name hack below.
                    knode_key = kg.add_node(
                            qnode.ids[0],
                            qnode.ids[0], #TODO: Once curies is fixed, make this a name.
                            qnode.categories[0].get_curie(),
                            )
                node_bindings[qnode_id] = [knode_key]
            else:
                wildcard_node = qnode
        if query_type == 'standard':
            for qedge_key, qedge in qg.edges.items():
                kedge_key = kg.add_edge(
                        node_bindings[qedge.subject][0],
                        node_bindings[qedge.object][0],
                        predicate=qedge.predicates[0].get_curie(),
                        )
                edge_bindings[qedge_key] = [kedge_key]
                # Add Attribute
                kg.edges[kedge_key].add_attribute(
                        attribute_type_id='Probability of Survival',
                        value=chp_query.truth_prob,
                        value_type_id=BIOLINK_HAS_CONFIDENCE_LEVEL_ENTITY.get_curie(),
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
                bad_wildcard = False
                for qnode_id, qnode in qg.nodes.items():
                    if qnode.categories[0] == BIOLINK_GENE_ENTITY and qnode.ids is None and (query_type == 'gene' or query_type == 'gene_two_hop'):
                        try:
                            knode_id = kg.add_node(
                                    wildcard,
                                    self.curies[BIOLINK_GENE_ENTITY.get_curie()][wildcard][0],
                                    qnode.categories[0].get_curie(),
                                    )
                            _node_bindings[qnode_id] = [knode_id]
                        except KeyError:
                            logger.info("Couldn't find {} in curies[{}]".format(wildcard, BIOLINK_GENE))
                            bad_wildcard = True
                    elif qnode.categories[0] == BIOLINK_DRUG_ENTITY and qnode.ids is None and (query_type == 'drug' or query_type == 'drug_two_hop'):
                        knode_id = kg.add_node(
                                wildcard,
                                self.curies[BIOLINK_DRUG_ENTITY.get_curie()][wildcard][0],
                                qnode.categories[0].get_curie(),
                                )
                        _node_bindings[qnode_id] = [knode_id]
                    else:
                        _node_bindings[qnode_id] = node_bindings[qnode_id]
                if bad_wildcard:
                    continue
                # Process edge bindings
                for qedge_id, qedge in qg.edges.items():
                    kedge_id = kg.add_edge(
                            _node_bindings[qedge.subject][0],
                            _node_bindings[qedge.object][0],
                            predicate=qedge.predicates[0],
                            )
                    kg.edges[kedge_id].add_attribute(
                            attribute_type_id='Contribution',
                            value=contrib,
                            value_type_id=BIOLINK_HAS_EVIDENCE_ENTITY.get_curie(),
                            description='Inference is derived from cancer data. Contribution reflects sensitivity of nodes to patient survival time. Contribution values closer to 1 are correlated with longer survival time and contribution values closer to -1 are correlated with shorter survival time.'
                            )
                    _edge_bindings[qedge_id] = [kedge_id]
                # Process node and edge binding results
                message.results.add_result(
                        _node_bindings,
                        _edge_bindings,
                        ) 
        return query
