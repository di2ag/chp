import unittest
import logging
import pickle
import json

from trapi_interface import TrapiInterface

logging.basicConfig(level=logging.INFO)

class TestDefaultHandler(unittest.TestCase):

    def setUp(self):
        # load in sample query graphs
        with open('query_samples/random_queries.pk', 'rb') as f_:
            self.queries = pickle.load(f_)
        #for i, q in enumerate(self.queries):
        #    print(i)
        #    print(q)
        #input()

    def test_curies(self):
        interface = TrapiInterface()
        curies = interface.get_curies()
        self.assertIsInstance(curies, dict)

    def test_predicates(self):
        interface = TrapiInterface()
        predicates = interface.get_predicates()
        self.assertIsInstance(predicates, dict)

    def test_default_single_query(self):
        # This is a non-simple query 
        message = self.queries[3]
        query = message["message"]
        interface = TrapiInterface(query=query, client_id='default')
        interface.build_chp_queries()
        interface.run_chp_queries()
        response = interface.construct_trapi_response()
        print(json.dumps(response, indent=2))

    def test_simple_single_query(self):
        # This is a simple query 
        message = self.queries[0]
        query = message["message"]
        interface = TrapiInterface(query=query, client_id='default')
        interface.build_chp_queries()
        interface.run_chp_queries()
        response = interface.construct_trapi_response()
        print(json.dumps(response, indent=2))

    def test_default_batch_query(self):
        # These are non-simple queries
        queries = [message["message"] for message in self.queries[2:8]]
        interface = TrapiInterface(query=queries, client_id='default')
        interface.build_chp_queries()
        interface.run_chp_queries()
        response = interface.construct_trapi_response()
        print(json.dumps(response, indent=2))

    def test_simple_batch_query(self):
        # These are simple queries
        queries = [message["message"] for message in self.queries[:2]]
        interface = TrapiInterface(query=queries, client_id='default')
        interface.build_chp_queries()
        interface.run_chp_queries()
        response = interface.construct_trapi_response()
        print(json.dumps(response, indent=2))

    def test_mix_batch_query(self):
        # These are mix batch of simple and default queries
        queries = [message["message"] for message in self.queries]
        interface = TrapiInterface(query=queries, client_id='default')
        interface.build_chp_queries()
        interface.run_chp_queries()
        response = interface.construct_trapi_response()
        print(json.dumps(response, indent=2))

    def test_mix_batch_reasoner_test_queries(self):
        with open('query_samples/test_reasoner_coulomb_queries.pk', 'rb') as f_:
            _queries = pickle.load(f_)
        queries = [message["message"] for message in _queries]
        interface = TrapiInterface(query=queries, client_id='default')
        interface.build_chp_queries()
        interface.run_chp_queries()
        response = interface.construct_trapi_response()
        print(json.dumps(response, indent=2))

class TestWildCardHandler(unittest.TestCase):

    def setUp(self):
        # load in sample query graphs
        with open('query_samples/random_gene_wildcard_queries.pk', 'rb') as f_:
            self.gene_queries = pickle.load(f_)
        with open('query_samples/random_drug_wildcard_queries.pk', 'rb') as f_:
            self.drug_queries = pickle.load(f_)
        #for i, q in enumerate(self.queries):
        #    print(i)
        #    print(q)
        #input()

    def test_single_gene_wildcard_query(self):
        message = self.gene_queries[0]
        query = message["message"]
        interface = TrapiInterface(query=query, client_id='default', max_results=10)
        interface.build_chp_queries()
        interface.run_chp_queries()
        response = interface.construct_trapi_response()
        print(json.dumps(response, indent=2))

    def test_single_drug_wildcard_query(self):
        message = self.drug_queries[0]
        query = message["message"]
        interface = TrapiInterface(query=query, client_id='default', max_results=10)
        interface.build_chp_queries()
        interface.run_chp_queries()
        response = interface.construct_trapi_response()
        print(json.dumps(response, indent=2))

    def test_batch_gene_wildcard_query(self):
        queries = [message["message"] for message in self.gene_queries]
        interface = TrapiInterface(query=queries, client_id='default', max_results=10)
        interface.build_chp_queries()
        interface.run_chp_queries()
        response = interface.construct_trapi_response()
        print(json.dumps(response, indent=2))

    def test_batch_drug_wildcard_query(self):
        queries = [message["message"] for message in self.gene_queries]
        interface = TrapiInterface(query=queries, client_id='default', max_results=10)
        interface.build_chp_queries()
        interface.run_chp_queries()
        response = interface.construct_trapi_response()
        print(json.dumps(response, indent=2))
