'''
Source code developed by DI2AG.
Thayer School of Engineering at Dartmouth College
Authors:    Dr. Eugene Santos, Jr
            Mr. Chase Yakaboski,
            Mr. Gregory Hyde,
            Dr. Keum Joo Kim
'''
import json
import itertools
import tqdm
import numpy as np
import logging
import csv
import uuid
from collections import defaultdict
import time

from trapi_model.biolink.constants import *
from trapi_model.logger import Logger as TrapiLogger
from trapi_model.meta_knowledge_graph import MetaKnowledgeGraph
from chp_utils.conflation import ConflationMap
from chp_utils.curie_database import CurieDatabase

from chp.trapi_handlers import BaseHandler, OneHopHandler
from chp.exceptions import *

# Setup logging
logger = logging.getLogger(__name__)


class TrapiInterface:
    def __init__(self,
                 hosts_filename=None,
                 num_processes_per_host=0,
                 bkb_handler=None,
                 joint_reasoner=None,
                 dynamic_reasoner=None,
                 trapi_version='1.1',
                ):
        self.hosts_filename = hosts_filename
        self.num_processes_per_host = num_processes_per_host
        self.bkb_handler = bkb_handler
        self.joint_reasoner = joint_reasoner
        self.dynamic_reasoner = dynamic_reasoner
        self.trapi_version = trapi_version

        # Get base handler for processing curies and meta kg requests
        self.base_handler = self._get_handler()
        self.curies_db = self._get_curies()
        self.curies = self.curies_db.curies
        self.meta_knowledge_graph = self._get_meta_knowledge_graph()
        self.conflation_map = self._get_conflation_map()

        # Initialize interface level logger
        self.logger = TrapiLogger()

    def setup_trapi_queries(self, trapi_queries):
        # Setup messages
        self.queries_dict = self._setup_messages(trapi_queries)
        # Initialize necessary handlers
        self.handlers = {}
        for message_type in self.queries_dict:
            self.handlers[message_type] = self._get_handler(message_type)
        return True 

    def _setup_messages(self, queries):
        queries_dict = defaultdict(list)
        for query in queries:
            try:
                message_type = self._determine_message_type(query.message)
                queries_dict[message_type].append(query)
            except Exception as ex:
                self.logger.debug(f'CHP core could not process derived query. {str(ex)}. Derived query graph: {query.message.query_graph.to_dict()}')
        # Check if any queries where setup/processed.
        if len(queries_dict) == 0:
            raise NoSupportedQueriesFound
        return queries_dict

    def _get_biolink_entity_if_supported_prefix(self, curie_prefix):
        for biolink_entity, meta_node in self.meta_knowledge_graph.nodes.items():
            if curie_prefix in meta_node.id_prefixes:
                return biolink_entity
        return None

    def _get_more_specific_biolink_entity(self, *biolink_entities):
        unsorted_entities = []
        for entity in biolink_entities:
            unsorted_entities.append((len(entity.get_ancestors()), entity))
        return sorted(unsorted_entities, reverse=True)[0][1]

    def _get_preferred_curie(self, normalized_node_dict):
        biolink_entities = normalized_node_dict["types"]
        normalized_curie_list = normalized_node_dict["equivalent_identifiers"]
        supported_biolink_entity = set.intersection(
                *[
                    set(biolink_entities),
                    set(self.meta_knowledge_graph.nodes.keys()),
                    ]
            )
        # Biolink category isn't supported
        if len(supported_biolink_entity) == 0:
            return None
        # If multiple supported biolink categories take the most specific.
        if len(supported_biolink_entity) > 1:
            supported_biolink_entity = self._get_more_specific_biolink_entity(*supported_biolink_type)
        # If there is only one supported entity just take that.
        else:
            supported_biolink_entity = supported_biolink_entity[0]
        preferred_prefix = self.meta_knowledge_graph.nodes[supported_biolink_entity].id_prefixes[0]
        for curie in normalized_curie_list:
            prefix = curie.split(':')[0]
            if prefix == preferred_prefix:
                return curie
        return None

    def _node_normalize_message(self, message):
        client = SriNodeNormalizerApiClient()
        # Collect node curies that need to be normalized
        node_curies = []
        for node_id, node in message.query_graph.items():
            # Grab curie prefix
            node_prefix = node.ids[0].split(':')[0]
            # Check CHP support
            curie_biolink_entity = self._get_biolink_entity_if_supported_prefix(node_prefix)
            # If curie prefix is not supported
            if curie_biolink_entity is None:
                node_curies.append(node.ids[0])
            # Check query category and supported category alignment
            elif curie_biolink_entity != node.categories[0]:
                if node.categories is None:
                    self.query.info('Interpretting node: {} as having category {}.'.format(node_id, curie_biolink_entity.get_curie()))
                else:
                    self.query.info(
                            'Mismatch in passed category for node {}. \
                                    Passed {} but deduced {}. Using deduced \
                                    category {} for reasoning.'.format(
                                        node_id,
                                        node.categories[0].get_curie(),
                                        curie_biolink_entity.get_curie(),
                                        curie_biolink_entity.get_curie(),
                                        )
                                    )
                node.categories[0] = curie_biolink_entity
        if len(node_curies) == 0:
            return message
        # Query the Node Normalizer
        normalized_nodes = client.get_normalized_nodes(node_curies)
        for origin_curie, normalized_node_dict in normalized_nodes.items():
            normalized_preferred_curie = self._get_preferred_curie(normalized_node_dict)
            if normalized_preferred_curie is None:
                self.query.info('Could not find a preferred curie in normalization of curie {}'.format(origin_curie))
                continue
            # Run Find and Replace over message
            message.find_and_replace(origin_curie, normalized_preferred_curie)
            # Save out curie mappings
            self.node_normalizer_mappings[origin_curie] = normalized_preferred_curie
            self.node_normalizer_inverse_mappings[normalized_preferred_curie] = origin_curie
        return message
            
    def _determine_message_type(self, message):
        """ checks for query message types. First checks node requirements to check for query type,
            then checks structures under the assumption of query type. Also updates error
            message for return to user

            :returns: a query type or None if there is a failure in matching query type
            :rtype: string or None
        """
        if message is None:
            raise UnidentifiedQueryType
        query_graph = message.query_graph

        # Check for standard or wildcard multihop query.
        gene_nodes = []
        disease_nodes = []
        drug_nodes = []
        phenotype_nodes = []
        wildcard_node_count = 0
        wildcard_node = None

        if message is not None:
            qg = message.query_graph
            for node_id, node in qg.nodes.items():
               if node.categories is not None:
                    if node.categories[0] == BIOLINK_GENE_ENTITY:
                        gene_nodes.append(node_id)
                        if node.ids is None:
                            wildcard_node_count += 1
                            wildcard_node = node_id
                        else:
                            for curie in node.ids:
                                if curie in self.curies[BIOLINK_GENE_ENTITY]:
                                    qg.nodes[node_id].set_ids(curie)
                                else:
                                    raise(UnidentifiedGeneCurie(node.ids))
                    elif node.categories[0] == BIOLINK_DRUG_ENTITY:
                        drug_nodes.append(node_id)
                        if node.ids is None:
                            wildcard_node_count += 1
                            wildcard_node = node_id
                        else:
                            for curie in node.ids:
                                if curie in self.curies[BIOLINK_DRUG_ENTITY]:
                                    qg.nodes[node_id].set_ids(curie)
                                else:
                                    raise(UnidentifiedDrugCurie(node.ids))
                    elif node.categories[0] == BIOLINK_DISEASE_ENTITY:
                        disease_nodes.append(node_id)
                    elif node.categories[0] == BIOLINK_PHENOTYPIC_FEATURE_ENTITY:
                        phenotype_nodes.append(node_id)
                        for curie in node.ids:
                            if curie in self.curies[BIOLINK_PHENOTYPIC_FEATURE_ENTITY]:
                                qg.nodes[node_id].set_ids(curie)
                            else:
                                raise(UnidentifiedPhenotypeCurie(node.ids))
                    else:
                        raise(UnidentifiedNode(node.categories[0]))

            num_total_nodes = len(gene_nodes) + len(disease_nodes) + len(drug_nodes) + len(phenotype_nodes)

            if wildcard_node_count > 1:
                raise(TooManyContributionNodes)
            if len(disease_nodes) > 1:
                raise(TooManyDiseaseNodes)
            if len(phenotype_nodes) > 1:
                raise(TooManyPhenotypeNodes)

            if num_total_nodes == 2 and len(phenotype_nodes) == 0:
                if self._check_one_hop_query(qg, gene_nodes, drug_nodes, disease_nodes, wildcard_node):
                    return 'onehop'
            else:
                raise(UnidentifiedQueryType)
    
    def get_conflation_map(self):
        return self.conflation_map

    def _get_conflation_map(self):
        return ConflationMap(conflation_map_filename=self.bkb_handler.conflation_map_path)

    def get_curies(self):
        return self.curies_db

    def _get_curies(self):
        """ Returns the available curies and their associated names.
        """
        curies_db = CurieDatabase(curies_filename=self.bkb_handler.curies_path)
        return curies_db

    def get_meta_knowledge_graph(self):
        return self.meta_knowledge_graph

    def _get_meta_knowledge_graph(self):
        """ Returns the meta knowledge graph.
        """
        return MetaKnowledgeGraph.load(
                self.trapi_version,
                None, 
                filename=self.base_handler.bkb_data_handler.meta_knowledge_graph_path,
                )

    def checkQuery(self):
        return True

    @staticmethod
    def check_predicate_support(predicate1, predicate2, support_inverse=True):
        if predicate1 == predicate2:
            is_inverse = False
            return True, is_inverse
        elif support_inverse is True and predicate2.get_inverse() is not None:
            if predicate1 == predicate2.get_inverse():
                is_inverse = True
                return True, is_inverse
        return False, False

    def _check_one_hop_query(self, query_graph, gene_nodes, drug_nodes, disease_nodes, wildcard_node):
        for edge_id, edge in query_graph.edges.items():

            is_valid, is_inverse = self.check_predicate_support(edge.predicates[0], BIOLINK_GENE_ASSOCIATED_WITH_CONDITION_ENTITY)
            if is_valid:
                if is_inverse:
                    if edge.subject not in disease_nodes  or edge.object not in gene_nodes or (wildcard_node is not None and edge.subject == wildcard_node):
                        raise(MalformedSubjectObjectOnGeneToDisease(edge_id))
                else:
                    if edge.subject not in gene_nodes or edge.object not in disease_nodes or (wildcard_node is not None and edge.object == wildcard_node):
                        raise(MalformedSubjectObjectOnGeneToDisease(edge_id))
                continue

            is_valid, is_inverse = self.check_predicate_support(edge.predicates[0], BIOLINK_TREATS_ENTITY)
            if is_valid:
                if is_inverse:
                    if edge.subject not in disease_nodes or edge.object not in drug_nodes or (wildcard_node is not None and edge.subject == wildcard_node):
                        raise(MalformedSubjectObjectOnDrugToDisease(edge_id))
                else:
                    if edge.subject not in drug_nodes or edge.object not in disease_nodes or (wildcard_node is not None and edge.object == wildcard_node):
                        raise(MalformedSubjectObjectOnDrugToDisease(edge_id))
                continue

            is_valid, is_inverse = self.check_predicate_support(edge.predicates[0], BIOLINK_GENETICALLY_INTERACTS_WITH_ENTITY)
            if is_valid:
                if is_inverse:
                    if edge.subject not in gene_nodes or edge.object not in gene_nodes:
                        raise(MalformedSubjectObjectOnGeneToGene(edge_id))
                else:
                    if edge.subject not in gene_nodes or edge.object not in gene_nodes:
                        raise(MalformedSubjectObjectOnGeneToGene(edge_id))

                continue

            is_valid, is_inverse = self.check_predicate_support(edge.predicates[0], BIOLINK_INTERACTS_WITH_ENTITY)
            if is_valid:
                if is_inverse:
                    if len(gene_nodes) != len(drug_nodes):
                        raise(MalformedSubjectObjectOnDrugGene(edge_id))
                else:
                    if len(gene_nodes) != len(drug_nodes):
                        raise(MalformedSubjectObjectOnDrugGene(edge_id))
                continue

            # unexpected edge
            raise(UnexpectedEdgeType(edge_id))
        return True

    def _get_handler(self, message_type=None):
        if message_type == 'onehop':
            return OneHopHandler(
                queries=self.queries_dict['onehop'],
                hosts_filename=self.hosts_filename,
                num_processes_per_host=self.num_processes_per_host,
                bkb_handler=self.bkb_handler,
                joint_reasoner=self.joint_reasoner,
                dynamic_reasoner=self.dynamic_reasoner,
            )
        elif message_type is None:
            return BaseHandler()
        else:
            raise ValueError('Unrecognized message type or unsupported message: {}'.format(message_type))

    def _order_response(self, results):
        _unordered_response = []
        for query_type, reasoner_type_results in results.items():
            for reasoner_type, query_results in reasoner_type_results.items():
                for query_id, result in query_results:
                    # If single result just return the response
                    if self.query_map is None:
                        return result
                    # Else put the results back in the appropriate order
                    _unordered_response.append((self.query_map.index(query_id), result))
        response = [result for _id, result in sorted(_unordered_response)]
        return response

    def build_chp_queries(self):
        built_chp_queries = {}
        for message_type, handler in self.handlers.items():
            logger.info('Building queries for {} type message(s).'.format(message_type))
            built_chp_queries[message_type] = handler.build_queries()
        return built_chp_queries

    def run_chp_queries(self):
        ran_chp_queries = {}
        for message_type, handler in self.handlers.items():
            logger.info('Running queries for {} type message(s).'.format(message_type))
            ran_chp_queries[message_type] = handler.run_queries()
        return ran_chp_queries

    def construct_trapi_responses(self):
        responses = []
        for message_type, handler in self.handlers.items():
            logger.info('Constructing TRAPI response(s) for {} type message(s).'.format(message_type))
            responses.extend(handler.construct_trapi_responses())
        return responses
