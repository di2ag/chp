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

from trapi_model.biolink.constants import *

from chp_data.bkb_handler import BkbDataHandler

from chp.mixins.trapi_handler.default_handler_mixin import DefaultHandlerMixin
from chp.mixins.trapi_handler.wildcard_handler_mixin import WildCardHandlerMixin
from chp.mixins.trapi_handler.one_hop_handler_mixin import OneHopHandlerMixin

# Base TRAPI Handler

class BaseHandler:
    """ Current there is only the one-hop handler.

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

    def __init__(self,
                 hosts_filename=None,
                 num_processes_per_host=0,
                 bkb_handler=None,
                 joint_reasoner=None,
                 dynamic_reasoner=None,
                 max_results=10):
        # Instantiate handler is one was not passed
        if bkb_handler is None:
            self.bkb_data_handler = BkbDataHandler(
                bkb_major_version='darwin',
                bkb_minor_version='2.0',
                disease='tcga',
            )
        else:
            self.bkb_data_handler = bkb_handler

        # Save off other attributes
        self.hosts_filename = hosts_filename
        self.num_processes_per_host = num_processes_per_host
        self.max_results = max_results
        self.joint_reasoner = joint_reasoner
        self.dynamic_reasoner = dynamic_reasoner

        # Run specific handler setup
        self._setup_handler()

        # Read in curies
        with open(self.bkb_data_handler.curies_path, 'r') as f_:
            self.curies = json.load(f_)
    
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

    def _setup_handler(self):
        pass

    def _setup_queries(self):
        """ Should setup up both batch and single queries for the handler and set an attribute called
        self.query_dict, where the keys are query types for instance 'simple' or 'default', 'gene' or 'drug', etc.
        """
        pass

    def check_query(self):
        """ Currently not implemented. Would check validity of query.
        """
        return True

    def build_queries(self):
        """ Parses over the sent query graph(s) to form a BKB query(s).

            :return: A CHP query dictionary where the keys are the type of query in the context of the handler.
            :rtype: dict
        """
        self.chp_query_dict = defaultdict(list)
        for message_type, queries in self.queries_dict.items():
            for _query in queries:
                self.chp_query_dict[message_type].append(self._extract_chp_query(_query, message_type))
        return self.chp_query_dict

    def _extract_chp_query(self, query, query_type):
        """ This method should be overwritten by the specific handler and should extract the chp query
            from an associated query graph.
        """
        pass

    def _run_query(self, chp_query, query_type):
        """ This method should be overwritten by the specific handler and should run and individual query
            and conduct any necessary post processing.
        """
        pass

    def run_queries(self):
        """ Runs built BKB query(s) in correspondence with the handlers _run_query function.
        """
        self.results = defaultdict(list)
        for message_type, chp_queries in self.chp_query_dict.items():
            for chp_query in chp_queries:
                self.results[message_type].append(self._run_query(chp_query, message_type))

    def construct_trapi_responses(self):
        """ Constructs the trapi responses for each query in correspondance with each handlers
            _construct_trapi_response function.
        """
        responses = []
        for (query_type, chp_queries), (_, queries) in zip(self.results.items(), self.queries_dict.items()):
            for chp_query, query in zip(chp_queries, queries):
                # Construct Message
                _response = self._construct_trapi_message(chp_query, query, query_type)
                # Add provenance
                _response = self.add_provenance_attributes(_response)
                # Add to list of responses
                responses.append(_response)
        return responses

    def merge_messages(self, responses):
        new_query = self.query.get_copy()
        master_message = new_query.message
        for _, responses in responses.items():
            for i, message in enumerate(responses):
                master_message.update(message.knowledge_graph, message.results)
        return new_query

    def add_provenance_attributes(self, query):
        """ Adds CHP and TCGA provenance attributes to every edge in the knowledge graph.
        """
        for kedge_id, kedge in query.message.knowledge_graph.edges.items():
            # Add CHP InfoRes Attribute
            kedge.add_attribute(
                    attribute_type_id=BIOLINK_PRIMARY_KNOWLEDGE_SOURCE_ENTITY.get_curie(),
                    value='infores:connections-hypothesis',
                    value_type_id=BIOLINK_INFORMATION_RESOURCE_ENTITY.get_curie(),
                    value_url='http://chp.thayer.dartmouth.edu',
                    description='The Connections Hypothesis Provider from NCATS Translator.',
                    attribute_source='infores:connections-hypothesis',
                    )
            # Add TCGA Supporting Data Source Attribute
            kedge.add_attribute(
                    attribute_type_id=BIOLINK_SUPPORTING_DATA_SOURCE_ENTITY.get_curie(),
                    value='infores:tcga',
                    value_type_id=BIOLINK_INFORMATION_RESOURCE_ENTITY.get_curie(),
                    value_url='https://portal.gdc.cancer.gov/',
                    description='The Cancer Genome Atlas provided by the GDC Data Portal.',
                    attribute_source='infores:gdc',
                    )
            kedge.add_attribute(
                    attribute_type_id=BIOLINK_SUPPORTING_DATA_SOURCE_ENTITY.get_curie(),
                    value='infores:tcga',
                    value_type_id=BIOLINK_INFORMATION_RESOURCE_ENTITY.get_curie(),
                    value_url='https://di2ag.github.io/chp/chp_documentation#chp-single-hop-contribution-analysis',
                    description='Documents how Connections Hypothesis Provider builds ranked contribution values over tcga data.',
                    attribute_source='infores:connections-hypothesis',
                    )
        return query

    def _get_curie_name(self, entity_type, curie):
        return self.curies[entity_type][curie]

    def _construct_trapi_message(self, chp_query, message, query_type):
        """ Should be overwitten by specific handler and return the response for a single query
            and its associated query_id which can be a uuid or just None if working with a single
            query.

            :return: (query_id, query_response)
            :rtype: tuple
        """
        pass

# Handler Mixins

class OneHopHandler(OneHopHandlerMixin, BaseHandler):
    def __init__(
            self,
            queries=None, 
            hosts_filename=None,
            num_processes_per_host=0,
            bkb_handler=None,
            joint_reasoner=None,
            dynamic_reasoner=None,
            max_results=10,
            ):
        self.queries = queries

        super(OneHopHandler, self).__init__(
            hosts_filename=hosts_filename,
            num_processes_per_host=num_processes_per_host,
            bkb_handler=bkb_handler,
            joint_reasoner=joint_reasoner,
            dynamic_reasoner=dynamic_reasoner,
            max_results=max_results,
            )
