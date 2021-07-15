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

from trapi_model.biolink.constants import *
from chp_data.bkb_handler import BkbDataHandler

from chp.query import Query as ChpQuery
from chp.reasoner import ChpDynamicReasoner
from pybkb.python_base.utils import get_operator, get_opposite_operator

# Setup logging
logger = logging.getLogger(__name__)

class WildCardHandlerMixin:
    def _setup_handler(self):
        # Only do the rest of this if a query is passed
        if self.messages is not None:
            # Setup messages
            self._setup_messages()

            # Instiatate Reasoners
            if self.dynamic_reasoner is None:
                self.dynamic_reasoner = ChpDynamicReasoner(
                    bkb_handler=self.bkb_data_handler,
                    hosts_filename=self.hosts_filename,
                    num_processes_per_host=self.num_processes_per_host)

    def _setup_messages(self):
        if type(self.messages) == list:
            self.message_dict = defaultdict(list)
            for message in self.messages:
                self.message_dict[self._get_wildcard_type(message)].append(message)

    def _get_wildcard_type(self, message):
        wildcard_type = None
        for node_id, node in message.query_graph.nodes.items():
            if node.ids is None:
                if wildcard_type is None:
                    wildcard_type = node.categories[0]
        if wildcard_type == BIOLINK_DRUG_ENTITY:
            return 'drug'
        elif wildcard_type == BIOLINK_GENE_ENTITY:
            return 'gene'
        else:
            raise ValueError('Did not understand wildcard type {}.'.format(wildcard_type))

    def _extract_chp_query(self, message, message_type):
        # Initialize CHP BKB Query
        chp_query = ChpQuery(reasoning_type='updating')
        # ensure we are using all nodes/edges
        total_nodes = 0
        total_edges = 0

        query_graph = message.query_graph
        # get phenotype node
        targets = list()
        acceptable_target_curies = ['EFO:0000714']
        self.implicit_survival_node = False
        for node_key in query_graph.nodes.keys():
            node = query_graph.nodes[node_key]
            if node.categories[0] == BIOLINK_PHENOTYPIC_FEATURE_ENTITY and node.ids[0] in acceptable_target_curies:
                target_id = node_key
                total_nodes += 1
        if total_nodes == 0:
            # Use Default Survival
            self.implicit_survival_node = True
            total_nodes += 1
            #acceptable_target_curies_print = ','.join(acceptable_target_curies)
            #sys.exit("Survival Node not found. Node category must be '{}' and id must be in: {}".format(Biolink(BIOLINK_PHENOTYPIC_FEATURE),
            #                                                                                            acceptable_target_curies_print))

        survival_value = 970
        survival_operator = '>='
        # get disease node info and ensure only 1 disease:
        acceptable_disease_curies = ['MONDO:0007254']
        for node_key in query_graph.nodes.keys():
            node = query_graph.nodes[node_key]
            if node.categories[0] == BIOLINK_DISEASE_ENTITY and node.ids[0] in acceptable_disease_curies:
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

        if self.implicit_survival_node:
            days=970
            qualifier = '>='
            total_edges += 1

        # set BKB target
        chp_query.add_dynamic_target('EFO:0000714', survival_operator, survival_value)
        truth_target = ('EFO:0000714', '{} {}'.format(survival_operator, survival_value))

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
                if message_type != 'gene':
                    gene_curie = node.ids[0]
                    if gene_curie in self.curies[BIOLINK_GENE_ENTITY.get_curie()]:
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
                if message_type != 'drug':
                    drug_curie = node.ids[0]
                    if drug_curie in self.curies[BIOLINK_DRUG_ENTITY.get_curie()]:
                        drug = drug_curie
                    chp_query.add_meta_evidence(drug, 'True')
                total_nodes += 1

        # Temporary solution to no evidence linking
        if len(chp_query.evidence.keys()) == 0 and len(chp_query.dynamic_evidence.keys()) == 0:
            self.no_evidence_probability_check = True
        else:
            self.no_evidence_probability_check = False

        # Set some other helpful attributes
        chp_query.truth_target = truth_target
        return chp_query

    def _run_query(self, chp_query, query_type):
        """ Runs build BKB query to calculate probability of survival.
            A probability is returned to specificy survival time w.r.t a drug.
            Contributions for each gene are calculuated and classified under
            their true/false target assignments.
        """

        # temporary solution to no evidence linking
        if not self.no_evidence_probability_check:
            if query_type == 'gene':
                chp_query = self.dynamic_reasoner.run_query(chp_query, bkb_type='drug')
            elif query_type == 'drug':
                chp_query = self.dynamic_reasoner.run_query(chp_query, bkb_type='gene')
            chp_res_dict = chp_query.result.process_updates()
            chp_res_norm_dict = chp_query.result.process_updates(normalize=True)
            #chp_query.result.summary()
            chp_res_contributions = chp_query.result.process_inode_contributions()
            try:
                chp_query.truth_prob = max([0, chp_res_dict[chp_query.truth_target[0]][chp_query.truth_target[1]]])
            except KeyError:
                # May need to come back and fix this.
                chp_query.truth_prob = -1

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

        else:
            # probability of survival
            num_survived = 0
            num_all = len(self.dynamic_reasoner.raw_patient_data.keys())
            str_op = chp_query.dynamic_targets['EFO:0000714']['op']
            opp_op = get_opposite_operator(str_op)
            op = get_operator(str_op)
            days = chp_query.dynamic_targets['EFO:0000714']['value']
            for patient, pat_dict in self.dynamic_reasoner.raw_patient_data.items():
                if op(pat_dict['survival_time'], days):
                    num_survived += 1
            chp_query.truth_prob = num_survived/num_all

            # patient_contributions
            patient_contributions = defaultdict(lambda: defaultdict(int))
            for patient, pat_dict in self.dynamic_reasoner.raw_patient_data.items():
                if op(pat_dict['survival_time'], days):
                    if num_survived == 0:
                        patient_contributions[('EFO:0000714', '{} {}'.format(str_op, days))][patient] = 0
                    else:
                        patient_contributions[('EFO:0000714', '{} {}'.format(str_op, days))][patient] = chp_query.truth_prob/num_survived
                else:
                    if num_survived == 0:
                        patient_contributions[('EFO:0000714', '{} {}'.format(opp_op, days))][patient] = (1-chp_query.truth_prob)/num_all
                    else:
                        patient_contributions[('EFO:0000714', '{} {}'.format(opp_op, days))][patient] = (1-chp_query.truth_prob)/(num_all-num_survived)

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
                    truth_target_gene_contrib += contrib / chp_query.truth_prob
                else:
                    nontruth_target_gene_contrib += contrib / (1 - chp_query.truth_prob)
            wildcard_contributions[curie]['relative'] = truth_target_gene_contrib - nontruth_target_gene_contrib

        chp_query.report = None
        chp_query.wildcard_contributions = wildcard_contributions

        return chp_query

    def _construct_trapi_message(self, chp_query, message, query_type=None):

        # update target node info and form edge pair combos for results graph

        qg = message.query_graph
        kg = message.knowledge_graph

        # Process Standard QUery as first result.
        # Process Nodes
        node_bindings = {}
        contrib_qg_id = None
        for qnode_key, qnode in qg.nodes.items():
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
                else:
                    knode_key = kg.add_node(
                            qnode.ids[0],
                            qnode.ids[0],
                            qnode.categories[0].get_curie(),
                            )
                node_bindings[qnode_key] = [knode_key]
        if not self.implicit_survival_node:
            # Process Edges
            edge_bindings = {}
            knowledge_edges = 0
            for qedge_key, qedge in qg.edges.items():
                if not qedge.subject in node_bindings or not qedge.object in node_bindings:
                    continue
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
                '''
                subject_node = kg['edges'][edge_key]['subject']
                if kg['edges'][edge_key]['predicate'] == BIOLINK_GENE_ENTITY_TO_DISEASE_PREDICATE, is_slot=True) and query['query_graph']['nodes'][subject_node]['category'] == BIOLINK_GENE_ENTITY and query_type == 'gene':
                    kg['edges'].pop(edge_key)
                elif kg['edges'][edge_key]['predicate'] == BIOLINK_CHEMICAL_TO_DISEASE_OR_PHENOTYPIC_FEATURE_PREDICATE, is_slot=True) and query['query_graph']['nodes'][subject_node]['category'] == BIOLINK_DRUG_ENTITY and query_type == 'drug':
                    kg['edges'].pop(edge_key)
                else:
                    kg_id = 'kge{}'.format(knowledge_edges)
                    knowledge_edges += 1
                    kg['edges'][kg_id] = kg['edges'].pop(edge_key)
                    kg['edges'][kg_id]['subject'] = node_pairs[kg['edges'][kg_id]['subject']]
                    kg['edges'][kg_id]['object'] = node_pairs[kg['edges'][kg_id]['object']]
                    edge_pairs[edge_key] = kg_id
                    if kg['edges'][kg_id]['predicate'] == BIOLINK_DISEASE_ENTITY_TO_PHENOTYPIC_FEATURE_PREDICATE, is_slot=True):
                        if 'properties' in kg['edges'][kg_id].keys():
                            kg['edges'][kg_id].pop('properties')
                        kg['edges'][kg_id]['attributes'] = [{'name':'Probability of Survival',
                                                             'type':BIOLINK_PROBABILITY,
                                                             'value':chp_query.truth_prob}]
                '''
            # Proces results
            message.results.add_result(
                    node_bindings,
                    edge_bindings,
                    )
            '''
            # Put first result of standard prob query of only curie nodes (i.e. no wildcard nodes where used as evidence)
            results = []
            results.append({'edge_bindings':dict(),
                            'node_bindings':dict()})
            for edge_pair_key in edge_pairs:
                results[0]['edge_bindings'][edge_pair_key] = [{ 'id': str(edge_pairs[edge_pair_key])}]
            for node_pair_key in node_pairs:
                results[0]['node_bindings'][node_pair_key] = [{ 'id': str(node_pairs[node_pair_key])}]
            '''
        #else:
        #    knowledge_edges = 0
        #    kg['edges'] = {}
        #    results = []

        # Build relative contribution results and added associated edges into knowledge graph
        unsorted_wildcard_contributions = []
        for wildcard, contrib_dict in chp_query.wildcard_contributions.items():
            unsorted_wildcard_contributions.append((contrib_dict['relative'], wildcard))
        sorted_wildcard_contributions = [(contrib,wildcard) for contrib, wildcard in sorted(unsorted_wildcard_contributions, key=lambda x: abs(x[0]), reverse=True)]

        for contrib, wildcard in sorted_wildcard_contributions[:self.max_results]:
            #TODO: Fix this!
            if wildcard == 'missing':
                continue
            #rg = copy.deepcopy(query["query_graph"])
            _node_bindings = {}
            _edge_bindings = {}
            # Process node bindings
            bad_wildcard = False
            for qnode_id, qnode in qg.nodes.items():
                if qnode.categories[0] == BIOLINK_GENE_ENTITY and query_type == 'gene':
                    try:
                        knode_id = kg.add_node(
                                wildcard,
                                self.curies[BIOLINK_GENE_ENTITY.get_curie()][wildcard][0],
                                qnode.categories[0].get_curie(),
                                )
                        _node_bindings[qnode_id] = [knode_id]
                    except KeyError:
                        logger.info("Couldn't find {} in curies[{}]".format(wildcard, BIOLINK_GENE_ENTITY.get_curie()))
                        bad_wildcard = True
                elif qnode.categories[0] == BIOLINK_DRUG_ENTITY and query_type == 'drug':
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
                subject_node = qedge.subject
                object_node = qedge.object
                if query_type == 'gene' and self.check_predicate_support(qedge.predicates[0], BIOLINK_GENE_ASSOCIATED_WITH_CONDITION_ENTITY) and qg.nodes[subject_node].categories[0] == BIOLINK_GENE_ENTITY:
                    kedge_id = kg.add_edge(
                            _node_bindings[qedge.subject][0],
                            _node_bindings[qedge.object][0],
                            predicate=qedge.predicates[0],
                            relation=qedge.relation,
                            )
                    kg.edges[kedge_id].add_attribute(
                            attribute_type_id='Contribution',
                            value=contrib,
                            value_type_id=BIOLINK_HAS_EVIDENCE_ENTITY.get_curie(),
                            )
                    _edge_bindings[qedge_id] = [kedge_id]
                elif query_type == 'gene' and self.check_predicate_support(qedge.predicates[0], BIOLINK_CONDITION_ASSOCIATED_WITH_GENE_ENTITY) and qg.nodes[object_node].categories[0] == BIOLINK_GENE_ENTITY:
                    kedge_id = kg.add_edge(
                            _node_bindings[qedge.subject][0],
                            _node_bindings[qedge.object][0],
                            predicate=qedge.predicates[0],
                            relation=qedge.relation,
                            )
                    kg.edges[kedge_id].add_attribute(
                            attribute_type_id='Contribution',
                            value=contrib,
                            value_type_id=BIOLINK_HAS_EVIDENCE_ENTITY.get_curie(),
                            )
                    _edge_bindings[qedge_id] = [kedge_id]
                elif query_type == 'drug' and self.check_predicate_support(qedge.predicates[0], BIOLINK_TREATS_ENTITY) and qg.nodes[subject_node].categories[0] == BIOLINK_DRUG_ENTITY:
                    kedge_id = kg.add_edge(
                            _node_bindings[qedge.subject][0],
                            _node_bindings[qedge.object][0],
                            predicate=qedge.predicates[0],
                            relation=qedge.relation,
                            )
                    kg.edges[kedge_id].add_attribute(
                            attribute_type_id='Contribution',
                            value=contrib,
                            value_type_id=BIOLINK_HAS_EVIDENCE_ENTITY.get_curie(),
                            )
                    _edge_bindings[qedge_id] = [kedge_id]
                elif query_type == 'drug' and self.check_predicate_support(qedge.predicates[0], BIOLINK_TREATED_BY_ENTITY) and qg.nodes[object_node].categories[0] == BIOLINK_DRUG_ENTITY:
                    kedge_id = kg.add_edge(
                            _node_bindings[qedge.subject][0],
                            _node_bindings[qedge.object][0],
                            predicate=qedge.predicates[0],
                            relation=qedge.relation,
                            )
                    kg.edges[kedge_id].add_attribute(
                            attribute_type_id='Contribution',
                            value=contrib,
                            value_type_id=BIOLINK_HAS_EVIDENCE_ENTITY.get_curie(),
                            )
                    _edge_bindings[qedge_id] = [kedge_id]
                else:
                    _edge_bindings[qedge_id] = edge_bindings[qedge_id]
            # Process node and edge binding results
            message.results.add_result(
                    _node_bindings,
                    _edge_bindings,
                    )
        return message
