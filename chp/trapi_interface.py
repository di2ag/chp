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

# Integrators
from chp.trapi_handlers import DefaultHandler, WildCardHandler, OneHopHandler
from chp.errors import *
from chp_data.trapi_constants import *

# Setup logging
logger = logging.getLogger(__name__)

# Helper functions

def parse_query_graph(query_graph):
    """ Will extract the parameters that were used to build the query from the client.
    """
    try:
        parsed = defaultdict(str)
        for node_id, node in query_graph["nodes"].items():
            if node["category"] == BIOLINK_PHENOTYPIC_FEATURE:
                parsed["outcome_name"] = node_id
            elif node["category"] == BIOLINK_DRUG:
                parsed["therapeutic"] = node_id
            elif node["category"] == BIOLINK_GENE:
                if 'genes' in parsed:
                    parsed["genes"].append(node_id)
                else:
                    parsed["genes"] = [node_id]
            elif node["category"] == BIOLINK_DISEASE:
                parsed["disease"] = node_id
            else:
                raise ValueError('Unrecognized category: {}'.format(node["category"]))
        # Sort the genes
        if 'genes' in parsed:
            parsed["genes"] = sorted(parsed["genes"])
        # Find the outome op and value
        for edge_id, edge in query_graph["edges"].items():
            if edge['predicate'] == BIOLINK_DISEASE_TO_PHENOTYPIC_FEATURE_PREDICATE:
                if 'properties' in edge.keys():
                    parsed["outcome_op"] = edge["properties"]["qualifier"]
                    parsed["outcome_value"] = edge["properties"]["days"]
                # default
                else:
                    parsed["outcome_op"] = ">="
                    parsed["outcome_value"] = 970
        return parsed
    except:
        return None

class TrapiInterface:
    def __init__(self,
                 query=None,
                 client_id=None,
                 hosts_filename=None,
                 num_processes_per_host=0,
                 max_results=100,
                 bkb_handler=None,
                 joint_reasoner=None,
                 dynamic_reasoner=None,
                ):
        self.client_id = client_id
        self.hosts_filename = hosts_filename
        self.num_processes_per_host = num_processes_per_host
        self.max_results = max_results
        self.bkb_handler = bkb_handler
        self.joint_reasoner = joint_reasoner
        self.dynamic_reasoner = dynamic_reasoner

        # Get default handler for processing curies and predicates requests
        self.handler = self._get_handler(None)
        self.curies = self.handler.curies

        if query is not None:
            # Analyze queries
            self.query_dict, self.query_map = self._setup_query(query)

            # Initialize necessary handlers
            self.handlers = {}
            for query_type in self.query_dict:
                self.handlers[query_type] = self._get_handler(query_type)


    def _setup_query(self, query):
        if type(query) == list:
            logger.info('Detected batch queries,')
            query_dict, query_map = self._setup_batch_queries(query)
            return query_dict, query_map
        else:
            logger.info('Detected single query.')
            query_dict = self._setup_single_query(query)
            return query_dict, None

    def _setup_batch_queries(self, queries):
        query_dict = defaultdict(list)
        query_map = []
        for query in queries:
            # Set a query ID for later reordering to match batch sequence
            query_type, query = self._determine_query_type(query)
            _id = uuid.uuid4()
            query["query_id"] = _id
            query_map.append(_id)
            query_dict[query_type].append(query)
        return query_dict, query_map

    def _setup_single_query(self, query):
        query_type, query = self._determine_query_type(query)
        _id = uuid.uuid4()
        query["query_id"] = _id
        return {query_type: [query]}

    def _determine_query_type(self, query):
        """ checks for query types. First checks node requirements to check for query type,
            then checks structures under the assumption of query type. Also updates error
            message for return to user

            :param query: single query or list of queries
            :type query: dict
            :returns: a query type or None if there is a failure in matching query type
            :rtype: string or None
        """

        gene_nodes = []
        disease_nodes = []
        drug_nodes = []
        phenotype_nodes = []
        wildcard_node_count = 0
        wildcard_node = None

        # check node criteria
        if query is not None:
            qg = query["query_graph"]
            for node_id, node in qg["nodes"].items():
                if "category" in node:
                    if node["category"] == BIOLINK_GENE:
                        gene_nodes.append(node_id)
                        if 'id' not in node:
                            wildcard_node_count += 1
                            wildcard_node = node_id
                        else:
                            if isinstance(node['id'], list):
                                found_curie = None
                                for curie in node['id']:
                                    if curie in self.curies[BIOLINK_GENE]:
                                        found_curie = curie
                                        query["query_graph"]["nodes"][node_id]['id'] = found_curie
                                if found_curie is None:
                                    raise(UnidentifiedGeneCurie(node['id']))
                            elif node['id'] not in self.curies[BIOLINK_GENE]:
                                raise(UnidentifiedGeneCurie(node['id']))
                    elif node["category"] == BIOLINK_DRUG:
                        drug_nodes.append(node_id)
                        if 'id' not in node:
                            wildcard_node_count += 1
                            wildcard_node = node_id
                        else:
                            if isinstance(node['id'], list):
                                found_curie = None
                                for curie in node['id']:
                                    if curie in self.curies[BIOLINK_DRUG]:
                                        found_curie = curie
                                        query["query_graph"]["nodes"][node_id]['id'] = found_curie
                                if found_curie is None:
                                    raise(UnidentifiedDrugCurie(node['id']))
                            elif node['id'] not in self.curies[BIOLINK_DRUG]:
                                raise(UnidentifiedDrugCurie(node['id']))
                    elif node["category"] == BIOLINK_DISEASE:
                        disease_nodes.append(node_id)
                    elif node["category"] == BIOLINK_PHENOTYPIC_FEATURE:
                        phenotype_nodes.append(node_id)
                        if isinstance(node['id'], list):
                                found_curie = None
                                for curie in node['id']:
                                    if curie in self.curies[BIOLINK_PHENOTYPIC_FEATURE]:
                                        found_curie = curie
                                        query["query_graph"]["nodes"][node_id]['id'] = found_curie
                                if found_curie is None:
                                    raise(UnidentifiedPhenotypeCurie(node['id']))
                        if node['id'] not in self.curies[BIOLINK_PHENOTYPIC_FEATURE]:
                                raise(UnidentifiedPhenotypeCurie(node['id']))

        if wildcard_node_count > 1:
            raise(TooManyContributionNodes)
        if len(disease_nodes) > 1:
            raise(TooManyDiseaseNodes)
        if len(phenotype_nodes) > 1:
            raise(TooManyPhenotypeNodes)

        if wildcard_node_count == 0 and len(phenotype_nodes) == 1 and len(disease_nodes) == 1:
            if self._check_default_query(qg, gene_nodes, drug_nodes, disease_nodes, phenotype_nodes):
                return 'default', query
        elif len(disease_nodes) == 1 and wildcard_node_count == 1:
            if self._check_wildcard_query(qg, gene_nodes, drug_nodes, disease_nodes, phenotype_nodes):
                return 'wildcard', query
        elif len(disease_nodes) == 0 and wildcard_node_count == 1:
            if self._check_one_hop_query(qg, gene_nodes, drug_nodes, disease_nodes, phenotype_nodes, wildcard_node):
                return 'onehop', query
        else:
            raise(UnidentifiedQueryType)

    def get_curies(self):
        """ Returns the available curies and their associated names.
        """
        return self.handler.curies

    def get_predicates(self):
        """ Returns the available predicates and their associated names.
        """
        with open(self.handler.bkb_data_handler.predicates_path, 'r') as predicates_file:
            return json.load(predicates_file)

    def checkQuery(self):
        return True

    def _check_wildcard_query(self, query, gene_nodes, drug_nodes, disease_nodes, phenotype_nodes):
        for edge_id, edge in query['edges'].items():
            if edge['predicate'] == BIOLINK_GENE_TO_DISEASE_PREDICATE:
                subject = edge['subject']
                object = edge['object']
                if subject not in gene_nodes or object not in disease_nodes:
                    raise(MalformedSubjectObjectOnGeneToDisease(edge_id))
            elif edge['predicate'] == BIOLINK_CHEMICAL_TO_DISEASE_OR_PHENOTYPIC_FEATURE_PREDICATE:
                subject = edge['subject']
                object = edge['object']
                if subject not in drug_nodes or object not in disease_nodes:
                    raise(MalformedSubjectObjectOnDrugToDisease(edge_id))
            elif edge['predicate'] == BIOLINK_DISEASE_TO_PHENOTYPIC_FEATURE_PREDICATE:
                subject = edge['subject']
                object = edge['object']
                if subject not in disease_nodes and object not in phenotype_nodes:
                    raise(MalformedSubjectObjectOnDiseaseToPhenotype(edge_id))
            elif edge['predicate'] == BIOLINK_CHEMICAL_TO_GENE_PREDICATE:
                raise(IncompatibleWildcardEdge(edge_id))
            else:
                raise(UnexpectedEdgeType(edge_id))
        return True

    def _check_one_hop_query(self, query, gene_nodes, drug_nodes, disease_nodes, phenotype_nodes, wildcard_node):
        for edge_id, edge in query['edges'].items():
            if edge['predicate'] == BIOLINK_GENE_TO_DISEASE_PREDICATE:
                raise(IncompatibleDrugGeneOneHopEdge(edge_id))
            elif edge['predicate'] == BIOLINK_CHEMICAL_TO_DISEASE_OR_PHENOTYPIC_FEATURE_PREDICATE:
                raise(IncompatibleDrugGeneOneHopEdge(edge_id))
            elif edge['predicate'] == BIOLINK_DISEASE_TO_PHENOTYPIC_FEATURE_PREDICATE:
                raise(IncompatibleDrugGeneOneHopEdge(edge_id))
            elif edge['predicate'] == BIOLINK_CHEMICAL_TO_GENE_PREDICATE:
                subject = edge['subject']
                object = edge['object']
                if object == wildcard_node:
                     raise(MalformedSubjectObjectOnDrugGene(edge_id))
            else:
                raise(UnexpectedEdgeType(edge_id))
        return True

    def _check_default_query(self, query, gene_nodes, drug_nodes, disease_nodes, phenotype_nodes):
        for edge_id, edge in query['edges'].items():
            if edge['predicate'] == BIOLINK_GENE_TO_DISEASE_PREDICATE:
                subject = edge['subject']
                object = edge['object']
                if subject not in gene_nodes or object not in disease_nodes:
                    raise(MalformedSubjectObjectOnGeneToDisease(edge_id))
            elif edge['predicate'] == BIOLINK_CHEMICAL_TO_DISEASE_OR_PHENOTYPIC_FEATURE_PREDICATE:
                subject = edge['subject']
                object = edge['object']
                if subject not in drug_nodes or object not in disease_nodes:
                    raise(MalformedSubjectObjectOnDrugToDisease(edge_id))
            elif edge['predicate'] == BIOLINK_DISEASE_TO_PHENOTYPIC_FEATURE_PREDICATE:
                subject = edge['subject']
                object = edge['object']
                if subject not in disease_nodes and object not in phenotype_nodes:
                    raise(MalformedSubjectObjectOnDiseaseToPhenotype(edge_id))
            elif edge['predicate'] == BIOLINK_CHEMICAL_TO_GENE_PREDICATE:
                raise(IncompatibleDefaultEdge(edge_id))
            else:
                raise(UnexpectedEdgeType(edge_id))
        return True

    def _get_handler(self, query_type):
        if query_type == 'default':
            return DefaultHandler(
                self.query_dict['default'],
                hosts_filename=self.hosts_filename,
                num_processes_per_host=self.num_processes_per_host,
                bkb_handler=self.bkb_handler,
                joint_reasoner=self.joint_reasoner,
                dynamic_reasoner=self.dynamic_reasoner,
            )
        elif query_type == 'wildcard':
            return WildCardHandler(
                self.query_dict['wildcard'],
                hosts_filename=self.hosts_filename,
                num_processes_per_host=self.num_processes_per_host,
                max_results=self.max_results,
                bkb_handler=self.bkb_handler,
                dynamic_reasoner=self.dynamic_reasoner,
            )
        elif query_type == 'onehop':
            return OneHopHandler(
                self.query_dict['onehop'],
                hosts_filename=self.hosts_filename,
                num_processes_per_host=self.num_processes_per_host,
                max_results=self.max_results,
                dynamic_reasoner=self.dynamic_reasoner,
            )
        elif query_type is None:
            return DefaultHandler(None)
        else:
            raise('Unrecognized query type or unsupported query: {}'.format(query_type))

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
        for query_type, handler in self.handlers.items():
            logger.info('Building queries for {} type query(s).'.format(query_type))
            built_chp_queries[query_type] = handler.build_queries()
        return built_chp_queries

    def run_chp_queries(self):
        ran_chp_queries = {}
        for query_type, handler in self.handlers.items():
            logger.info('Running queries for {} type query(s).'.format(query_type))
            ran_chp_queries[query_type] = handler.run_queries()
        return ran_chp_queries

    def construct_trapi_response(self):
        results = {}
        for query_type, handler in self.handlers.items():
            logger.info('Constructing TRAPI response(s) for {} type query(s).'.format(query_type))
            results[query_type] = handler.construct_trapi_response()
        return self._order_response(results)
