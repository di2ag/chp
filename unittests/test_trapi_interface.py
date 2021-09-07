import unittest
import logging
import pickle
import json
import copy
import sys

import trapi_model
trapi_model.set_biolink_debug_mode(False)
from trapi_model.query import Query
from chp_data.bkb_handler import BkbDataHandler

from chp.trapi_interface import TrapiInterface
from chp.reasoner import ChpJointReasoner, ChpDynamicReasoner
from chp.exceptions import *

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
stream_handler = logging.StreamHandler(sys.stdout)
logger.addHandler(stream_handler)

#logging.basicConfig(level=logging.INFO)
#logger = logging.getLogger()
#logger.setLevel(logging.INFO)
logger_root = logging.getLogger()
logger_root.setLevel(logging.INFO)


class TestBaseHandler(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        super(TestBaseHandler, cls).setUpClass()
        
        cls.bkb_handler = BkbDataHandler()
        cls.dynamic_reasoner = ChpDynamicReasoner(cls.bkb_handler)
        cls.joint_reasoner = ChpJointReasoner(cls.bkb_handler)

    def test_curies(self):
        interface = TrapiInterface(
                bkb_handler=self.bkb_handler,
                dynamic_reasoner=self.dynamic_reasoner,
                joint_reasoner=self.joint_reasoner,
                )
        curies = interface.get_curies()
        #print(curies.json())

    def test_meta_knowledge_graph(self):
        interface = TrapiInterface(
                bkb_handler=self.bkb_handler,
                dynamic_reasoner=self.dynamic_reasoner,
                joint_reasoner=self.joint_reasoner,
                )
        meta_kg = interface.get_meta_knowledge_graph()
        #print(meta_kg.json())

class TestOneHopHandler(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        super(TestOneHopHandler, cls).setUpClass()
        # load in sample query graphs
        with open('query_samples/onehop/standard_queries.json', 'r') as f_:
            cls.standard_queries = json.load(f_)
        with open('query_samples/onehop/wildcard_queries.json', 'r') as f_:
            cls.wildcard_queries = json.load(f_)
        cls.bkb_handler = BkbDataHandler()
        cls.dynamic_reasoner = ChpDynamicReasoner(cls.bkb_handler)
        cls.joint_reasoner = ChpJointReasoner(cls.bkb_handler)


    def get_responses(self, queries=None, trapi_queries=None):
        # Initialize interface
        interface = TrapiInterface(
                bkb_handler=self.bkb_handler,
                dynamic_reasoner=self.dynamic_reasoner,
                joint_reasoner=self.joint_reasoner,
                )
        # Load queries
        if trapi_queries is None:
            trapi_queries = [Query.load(query["trapi_version"], None, query=query) for query in queries]
        # Process trapi query
        interface.setup_trapi_queries(trapi_queries)
        # Build CHP queries
        interface.build_chp_queries()
        # Run CHP queries
        interface.run_chp_queries()
        # Get Responses
        responses = interface.construct_trapi_responses()
        return responses


    def test_standard_query(self):
        standard_queries = copy.deepcopy(self.standard_queries)
        descriptions = [query.pop("test_description", None) for query in self.standard_queries]
        responses = self.get_responses(queries=standard_queries)
    
    def test_inverse_query(self):
        standard_queries = copy.deepcopy(self.standard_queries)
        descriptions = [query.pop("test_description", None) for query in self.standard_queries]
        # Make inverse queries
        trapi_queries = [Query.load(query["trapi_version"], None, query=query) for query in standard_queries] 
        for query in trapi_queries:
            for edge_id, edge in query.message.query_graph.edges.items():
                predicate = edge.predicates[0]
                inverse = edge.predicates[0].get_inverse()
                if inverse is not None:
                    edge.set_predicates(inverse)
                    # Switch subject and object
                    edge_subject = copy.deepcopy(edge.subject)
                    edge_object = copy.deepcopy(edge.object)
                    edge.subject = edge_object
                    edge.object = edge_subject
        responses = self.get_responses(trapi_queries=trapi_queries)
    def test_wildcard_query(self):
        wildcard_queries = copy.deepcopy(self.wildcard_queries)
        descriptions = [query.pop("test_description", None) for query in self.wildcard_queries]
        responses = self.get_responses(queries=wildcard_queries)


if __name__ == '__main__':
    unittest.main()
