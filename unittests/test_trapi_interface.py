import unittest
import logging
import pickle
import json

from chp.trapi_interface import TrapiInterface
from chp.errors import *

logging.basicConfig(level=logging.INFO)

class TestDefaultHandler(unittest.TestCase):

    def setUp(self):
        # load in sample query graphs
        with open('query_samples/random_queries.pk', 'rb') as f_:
            self.queries = pickle.load(f_)

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

    def test_simple_single_query(self):
        # This is a simple query 
        message = self.queries[1]
        query = message["message"]
        interface = TrapiInterface(query=query, client_id='default')
        interface.build_chp_queries()
        interface.run_chp_queries()
        response = interface.construct_trapi_response()
        #print(json.dumps(response, indent=2))

    def test_default_batch_query(self):
        # These are non-simple queries
        queries = [message["message"] for message in self.queries[2:8]]
        interface = TrapiInterface(query=queries, client_id='default')
        interface.build_chp_queries()
        interface.run_chp_queries()
        response = interface.construct_trapi_response()

    def test_simple_batch_query(self):
        # These are simple queries
        queries = [message["message"] for message in self.queries[:2]]
        interface = TrapiInterface(query=queries, client_id='default')
        interface.build_chp_queries()
        interface.run_chp_queries()
        response = interface.construct_trapi_response()

    def test_mix_batch_query(self):
        # These are mix batch of simple and default queries
        queries = [message["message"] for message in self.queries]
        interface = TrapiInterface(query=queries, client_id='default')
        interface.build_chp_queries()
        interface.run_chp_queries()
        response = interface.construct_trapi_response()

    def test_mix_batch_reasoner_test_queries(self):
        with open('query_samples/test_reasoner_coulomb_queries.pk', 'rb') as f_:
            _queries = pickle.load(f_)
        queries = [message["message"] for message in _queries]
        interface = TrapiInterface(query=queries, client_id='default')
        interface.build_chp_queries()
        interface.run_chp_queries()
        response = interface.construct_trapi_response()

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

    def test_single_drug_wildcard_query(self):
        message = self.drug_queries[0]
        query = message["message"]
        interface = TrapiInterface(query=query, client_id='default', max_results=10)
        interface.build_chp_queries()
        interface.run_chp_queries()
        response = interface.construct_trapi_response()

    def test_batch_gene_wildcard_query(self):
        queries = [message["message"] for message in self.gene_queries]
        interface = TrapiInterface(query=queries, client_id='default', max_results=10)
        interface.build_chp_queries()
        interface.run_chp_queries()
        response = interface.construct_trapi_response()

    def test_batch_drug_wildcard_query(self):
        queries = [message["message"] for message in self.gene_queries]
        interface = TrapiInterface(query=queries, client_id='default', max_results=10)
        interface.build_chp_queries()
        interface.run_chp_queries()
        response = interface.construct_trapi_response()

class TestOneHopHandler(unittest.TestCase):

    def setUp(self):
        # load in sample query graphs
        with open('query_samples/random_gene_one_hop_queries.pk', 'rb') as f_:
            self.gene_queries = pickle.load(f_)
        with open('query_samples/random_drug_one_hop_queries.pk', 'rb') as f_:
            self.drug_queries = pickle.load(f_)

    def test_single_gene_onehop_query(self):
        message = self.gene_queries[0]
        query = message["message"]
        interface = TrapiInterface(query=query, client_id='default', max_results=10)
        interface.build_chp_queries()
        interface.run_chp_queries()
        response = interface.construct_trapi_response()
        #print(json.dumps(response, indent=2))

    def test_single_drug_onehop_query(self):
        message = self.drug_queries[0]
        query = message["message"]
        interface = TrapiInterface(query=query, client_id='default', max_results=10)
        interface.build_chp_queries()
        interface.run_chp_queries()
        response = interface.construct_trapi_response()
        #print(json.dumps(response, indent=2))

    def test_batch_gene_onehop_query(self):
        queries = [message["message"] for message in self.gene_queries]
        interface = TrapiInterface(query=queries, client_id='default', max_results=10)
        interface.build_chp_queries()
        interface.run_chp_queries()
        response = interface.construct_trapi_response()

    def test_batch_drug_onehop_query(self):
        queries = [message["message"] for message in self.gene_queries]
        interface = TrapiInterface(query=queries, client_id='default', max_results=10)
        interface.build_chp_queries()
        interface.run_chp_queries()
        response = interface.construct_trapi_response()

class TestHandlerErrors(unittest.TestCase):

    def test_more_than_one_contribution(self):
        with open('query_samples/error_samples/test_more_than_one_contribution.pk', 'rb') as f_:
            query = pickle.load(f_)
        with self.assertRaises(TooManyContributionNodes) as context:
            interface = TrapiInterface(query=query)

    def test_more_than_one_disease(self):
        with open('query_samples/error_samples/test_more_than_one_disease.pk', 'rb') as f_:
            query = pickle.load(f_)
        with self.assertRaises(TooManyDiseaseNodes) as context:
            interface = TrapiInterface(query=query)

    def test_more_than_one_phenotype(self):
        with open('query_samples/error_samples/test_more_than_one_phenotype.pk', 'rb') as f_:
            query = pickle.load(f_)
        with self.assertRaises(TooManyPhenotypeNodes) as context:
            interface = TrapiInterface(query=query)

    def test_more_than_one_disease(self):
        with open('query_samples/error_samples/test_more_than_one_disease.pk', 'rb') as f_:
            query = pickle.load(f_)
        with self.assertRaises(TooManyDiseaseNodes) as context:
            interface = TrapiInterface(query=query)

    def test_no_disease(self):
        with open('query_samples/error_samples/test_no_disease.pk', 'rb') as f_:
            query = pickle.load(f_)
        with self.assertRaises(UnidentifiedQueryType) as context:
            interface = TrapiInterface(query=query)

    def test_no_target(self):
        with open('query_samples/error_samples/test_no_target.pk', 'rb') as f_:
            query = pickle.load(f_)
        with self.assertRaises(UnidentifiedQueryType) as context:
            interface = TrapiInterface(query=query)

    def test_illegal_drug_to_disease_default(self):
        with open('query_samples/error_samples/test_illegal_drug_to_disease_default.pk', 'rb') as f_:
            query = pickle.load(f_)
        with self.assertRaises(MalformedSubjectObjectOnDrugToDisease) as context:
            interface = TrapiInterface(query=query)

    def test_illegal_gene_to_disease_default(self):
        with open('query_samples/error_samples/test_illegal_gene_to_disease_default.pk', 'rb') as f_:
            query = pickle.load(f_)
        with self.assertRaises(MalformedSubjectObjectOnGeneToDisease) as context:
            interface = TrapiInterface(query=query)

    def test_illegal_disease_to_phenotype_default(self):
        with open('query_samples/error_samples/test_illegal_disease_to_phenotype_default.pk', 'rb') as f_:
            query = pickle.load(f_)
        with self.assertRaises(MalformedSubjectObjectOnDiseaseToPhenotype) as context:
            interface = TrapiInterface(query=query)

    def test_illegal_edge_default(self):
        with open('query_samples/error_samples/test_illegal_edge_default.pk', 'rb') as f_:
            query = pickle.load(f_)
        with self.assertRaises(IncompatibleDefaultEdge) as context:
            interface = TrapiInterface(query=query)

    def test_unknown_edge_default(self):
        with open('query_samples/error_samples/test_unknown_edge_default.pk', 'rb') as f_:
            query = pickle.load(f_)
        with self.assertRaises(UnexpectedEdgeType) as context:
            interface = TrapiInterface(query=query)

    def test_illegal_gene_to_disease_wildcard(self):
        with open('query_samples/error_samples/test_illegal_gene_to_disease_wildcard.pk', 'rb') as f_:
            query = pickle.load(f_)
        with self.assertRaises(MalformedSubjectObjectOnGeneToDisease) as context:
            interface = TrapiInterface(query=query)

    def test_illegal_drug_to_disease_wildcard(self):
        with open('query_samples/error_samples/test_illegal_drug_to_disease_wildcard.pk', 'rb') as f_:
            query = pickle.load(f_)
        with self.assertRaises(MalformedSubjectObjectOnDrugToDisease) as context:
            interface = TrapiInterface(query=query)

    def test_illegal_disease_to_phenotype_wildcard(self):
        with open('query_samples/error_samples/test_illegal_disease_to_phenotype_wildcard.pk', 'rb') as f_:
            query = pickle.load(f_)
        with self.assertRaises(MalformedSubjectObjectOnDiseaseToPhenotype) as context:
            interface = TrapiInterface(query=query)

    def test_illegal_edge_wildcard(self):
        with open('query_samples/error_samples/test_illegal_edge_wildcard.pk', 'rb') as f_:
            query = pickle.load(f_)
        with self.assertRaises(IncompatibleWildcardEdge) as context:
            interface = TrapiInterface(query=query)

    def test_unknown_edge_wildcard(self):
        with open('query_samples/error_samples/test_unknown_edge_wildcard.pk', 'rb') as f_:
            query = pickle.load(f_)
        with self.assertRaises(UnexpectedEdgeType) as context:
            interface = TrapiInterface(query=query)

    def test_illegal_drug_to_disease_one_hop(self):
        with open('query_samples/error_samples/test_illegal_drug_to_disease_one_hop.pk', 'rb') as f_:
            query = pickle.load(f_)
        with self.assertRaises(IncompatibleDrugGeneOneHopEdge) as context:
            interface = TrapiInterface(query=query)

    def test_illegal_gene_to_drug_one_hop(self):
        with open('query_samples/error_samples/test_illegal_gene_to_drug_one_hop.pk', 'rb') as f_:
            query = pickle.load(f_)
        with self.assertRaises(IncompatibleDrugGeneOneHopEdge) as context:
            interface = TrapiInterface(query=query)

    def test_illegal_disease_to_phenotype_one_hop(self):
        with open('query_samples/error_samples/test_illegal_disease_to_phenotype_one_hop.pk', 'rb') as f_:
            query = pickle.load(f_)
        with self.assertRaises(IncompatibleDrugGeneOneHopEdge) as context:
            interface = TrapiInterface(query=query)

    def test_backwards_contribution_node_one_hop(self):
        with open('query_samples/error_samples/test_backwards_contribution_node_one_hop.pk', 'rb') as f_:
            query = pickle.load(f_)
        with self.assertRaises(MalformedSubjectObjectOnDrugGene) as context:
            interface = TrapiInterface(query=query)

    def test_unknown_edge_one_hop(self):
        with open('query_samples/error_samples/test_unknown_edge_one_hop.pk', 'rb') as f_:
            query = pickle.load(f_)
        with self.assertRaises(UnexpectedEdgeType) as context:
            interface = TrapiInterface(query=query)



if __name__ == '__main__':
    unittest.main()
