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

from chp.mixins.trapi_handler.default_handler_mixin import DefaultHandlerMixin
from chp.mixins.trapi_handler.wildcard_handler_mixin import WildCardHandlerMixin
from chp.mixins.trapi_handler.one_hop_handler_mixin import OneHopHandlerMixin

# Base TRAPI Handler

class BaseHandler:
    """WildCardHandler is the handler for gene wildcards. That is
        query graphs (QGs) that consists of 4 nodes and 3 edges.

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
                 query,
                 hosts_filename=None,
                 num_processes_per_host=0,
                 max_results=100,
                 bkb_handler=None,
                 joint_reasoner=None,
                 dynamic_reasoner=None):
        # Save initial passed query(s)
        self.init_query = query
        # Instantiate handler is one was not passed
        if bkb_handler is None:
            self.bkb_data_handler = BkbDataHandler(
                bkb_major_version='coulomb',
                bkb_minor_version='1.0'
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

    def _handler_setup(self):
        pass

    def _setup_queries(self):
        """ Should setup up both batch and single queries for the handler and set an attribute called
        self.query_dict, where the keys are query types for instance 'simple' or 'default', 'gene' or 'drug', etc.
        """
        pass

    def _setup_single_query(self, query):
        if 'knowledge_graph' not in query:
            query["knowledge_graph"] = {
                "edges": {},
                "nodes": {},
            }
        if 'results' not in query:
            query["results"] = []
        return query

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
        for query_type, query in self.query_dict.items():
            for _query in query:
                self.chp_query_dict[query_type].append(self._extract_chp_query(_query, query_type))
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
        for query_type, chp_queries in self.chp_query_dict.items():
            for chp_query in chp_queries:
                self.results[query_type].append(self._run_query(chp_query, query_type))

    def construct_trapi_response(self):
        """ Constructs the trapi responses for each query in correspondance with each handlers
            _construct_trapi_response function.
        """
        responses = defaultdict(list)
        for query_type, chp_queries in self.results.items():
            for chp_query in chp_queries:
                responses[query_type].append(self._construct_trapi_response(chp_query, query_type))
        return responses

    def _get_curie_name(self, entity_type, curie):
        return self.curies[entity_type][curie]

    def _construct_trapi_response(self, chp_query, query_type):
        """ Should be overwitten by specific handler and return the response for a single query
            and its associated query_id which can be a uuid or just None if working with a single
            query.

            :return: (query_id, query_response)
            :rtype: tuple
        """
        pass

# Handler Mixins

class DefaultHandler(DefaultHandlerMixin, BaseHandler):
    pass

class WildCardHandler(WildCardHandlerMixin, BaseHandler):
    pass

class OneHopHandler(OneHopHandlerMixin, BaseHandler):
    pass
