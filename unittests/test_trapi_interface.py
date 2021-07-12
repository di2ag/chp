import unittest
import logging
import pickle
import json
import copy
import sys

import trapi_model
trapi_model.set_biolink_debug_mode(False)
from trapi_model.query import Query

from chp.trapi_interface import TrapiInterface
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

class TestDefaultHandler(unittest.TestCase):

    def setUp(self):
        # load in sample query graphs
        with open('query_samples/random_queries.pk', 'rb') as f_:
            self.queries = pickle.load(f_)
        with open('query_samples/random_batch_queries.pk', 'rb') as f_:
            self.batch_queries = pickle.load(f_)

    def test_curies(self):
        interface = TrapiInterface()
        curies = interface.get_curies()
        self.assertIsInstance(curies, dict)

    def test_predicates(self):
        interface = TrapiInterface()
        predicates = interface.get_predicates()
        self.assertIsInstance(predicates, dict)

    def test_meta_knowledge_graph(self):
        interface = TrapiInterface()
        meta_kg = interface.get_meta_knowledge_graph()
        self.assertIsInstance(meta_kg, dict)

    def test_default_single_query(self):
        # This is a non-simple query
        logger.info('Running single query test.')
        for trapi_version, queries in self.queries.items():
            query = Query.load(trapi_version, None, query=queries[3])
            interface = TrapiInterface(query=query)
            interface.build_chp_queries()
            interface.run_chp_queries()
            response = interface.construct_trapi_response()
    
    def test_inverse_query(self):
        # This is a simple query 
        logger.info('Running default inverse query test.')
        for trapi_version, queries in self.queries.items():
            query = Query.load(trapi_version, None, query=queries[1])
            for edge_id, edge in query.message.query_graph.edges.items():
                predicate = edge.predicates[0]
                inverse = edge.predicates[0].get_inverse()
                edge.set_predicates(inverse)
                subject = edge.subject
                object = edge.object
                edge.subject = object
                edge.object = subject
            interface = TrapiInterface(query=query)
            interface.build_chp_queries()
            interface.run_chp_queries()
            response = interface.construct_trapi_response()

    def test_simple_single_query(self):
        # This is a simple query 
        logger.info('Running single simple query test.')
        for trapi_version, queries in self.queries.items():
            query = Query.load(trapi_version, None, query=queries[1])
            interface = TrapiInterface(query=query)
            interface.build_chp_queries()
            interface.run_chp_queries()
            response = interface.construct_trapi_response()
        
    def test_default_batch_query(self):
        logger.info('Running batch query test.')
        # These are non-simple queries
        for trapi_version, queries in self.batch_queries.items():
            query = Query.load(trapi_version, None, query=queries[1])
            interface = TrapiInterface(query=query)
            interface.build_chp_queries()
            interface.run_chp_queries()
            response = interface.construct_trapi_response()

    def test_simple_batch_query(self):
        # These are simple queries
        logger.info('Running batch simple query test.')
        for trapi_version, queries in self.batch_queries.items():
            query = Query.load(trapi_version, None, query=queries[0])
            interface = TrapiInterface(query=query)
            interface.build_chp_queries()
            interface.run_chp_queries()
            response = interface.construct_trapi_response()

class TestWildCardHandler(unittest.TestCase):

    def setUp(self):
        # load in sample query graphs
        with open('query_samples/random_gene_wildcard_queries.pk', 'rb') as f_:
            self.gene_queries = pickle.load(f_)
        with open('query_samples/random_gene_wildcard_batch_queries.pk', 'rb') as f_:
            self.gene_batch_queries = pickle.load(f_)
        with open('query_samples/random_drug_wildcard_queries.pk', 'rb') as f_:
            self.drug_queries = pickle.load(f_)
        with open('query_samples/random_drug_wildcard_batch_queries.pk', 'rb') as f_:
            self.drug_batch_queries = pickle.load(f_)

    def test_single_gene_wildcard_query(self):
        for trapi_version, queries in self.gene_queries.items():
            query = Query.load(trapi_version, None, query=queries[0])
            interface = TrapiInterface(query=query)
            interface.build_chp_queries()
            interface.run_chp_queries()
            response = interface.construct_trapi_response()

    def test_inverse_wildcard_query(self):
        for trapi_version, queries in self.gene_queries.items():
            query = Query.load(trapi_version, None, query=queries[0])
            for edge_id, edge in query.message.query_graph.edges.items():
                predicate = edge.predicates[0]
                inverse = edge.predicates[0].get_inverse()
                edge.set_predicates(inverse)
                obj = edge.object
                sub = edge.subject
                edge.object = sub
                edge.subject = obj
            interface = TrapiInterface(query=query)
            interface.build_chp_queries()
            interface.run_chp_queries()
            response = interface.construct_trapi_response()


    def test_single_drug_wildcard_query(self):
        for trapi_version, queries in self.drug_queries.items():
            query = Query.load(trapi_version, None, query=queries[0])
            interface = TrapiInterface(query=query)
            interface.build_chp_queries()
            interface.run_chp_queries()
            response = interface.construct_trapi_response()

    def test_batch_gene_wildcard_query(self):
        for trapi_version, queries in self.gene_batch_queries.items():
            query = Query.load(trapi_version, None, query=queries[0])
            interface = TrapiInterface(query=query)
            interface.build_chp_queries()
            interface.run_chp_queries()
            response = interface.construct_trapi_response()
    
    def test_batch_drug_wildcard_query(self):
        for trapi_version, queries in self.drug_batch_queries.items():
            query = Query.load(trapi_version, None, query=queries[0])
            interface = TrapiInterface(query=query)
            interface.build_chp_queries()
            interface.run_chp_queries()
            response = interface.construct_trapi_response()

class TestOneHopHandler(unittest.TestCase):

    def setUp(self):
        # load in sample query graphs
        with open('query_samples/standard_single_onehop_queries.pk', 'rb') as f_:
            self.standard_single_queries = pickle.load(f_)
        with open('query_samples/standard_batch_onehop_queries.pk', 'rb') as f_:
            self.standard_batch_queries = pickle.load(f_)
        with open('query_samples/wildcard_single_onehop_queries.pk', 'rb') as f_:
            self.wildcard_single_queries = pickle.load(f_)
        with open('query_samples/wildcard_batch_onehop_queries.pk', 'rb') as f_:
            self.wildcard_batch_queries = pickle.load(f_)
        with open('query_samples/gene_to_gene_onehop_query.pk', 'rb') as f_:
            self.gene_to_gene_query = pickle.load(f_)

    def test_standard_single_onehop_query(self):
        for trapi_version, queries in self.standard_single_queries.items():
            for name, query_dict in queries.items():
                #if name != 'gene_to_disease_proxy_context':
                #    continue
                query = Query.load(trapi_version, None, query=query_dict)
                interface = TrapiInterface(query=query)
                interface.build_chp_queries()
                interface.run_chp_queries()
                response = interface.construct_trapi_response()
    
    def test_inverse_onehop_query(self):
        for trapi_version, queries in self.standard_single_queries.items():
            for name, query_dict in queries.items():
                #if name != 'gene_to_disease_proxy_context':
                #    continue
                query = Query.load(trapi_version, None, query=query_dict)
                for edge_id, edge in query.message.query_graph.edges.items():
                    predicate = edge.predicates[0]
                    inverse = edge.predicates[0].get_inverse()
                    if inverse is not None:
                        edge.set_predicates(inverse)
                    subject = edge.subject
                    object = edge.object
                    edge.subject = object
                    edge.object = subject
                interface = TrapiInterface(query=query)
                interface.build_chp_queries()
                interface.run_chp_queries()
                response = interface.construct_trapi_response()
            
    def test_standard_batch_onehop_query(self):
        for trapi_version, queries in self.standard_batch_queries.items():
            for name, query_dict in queries.items():
                #if name != 'gene_to_disease_proxy_context':
                #    continue
                query = Query.load(trapi_version, None, query=query_dict)
                interface = TrapiInterface(query=query)
                interface.build_chp_queries()
                interface.run_chp_queries()
                response = interface.construct_trapi_response()
    
    def test_wildcard_single_onehop_query(self):
        for trapi_version, queries in self.wildcard_single_queries.items():
            for name, query_dict in queries.items():
                #if name != 'gene_to_disease_proxy_context':
                #    continue
                query = Query.load(trapi_version, None, query=query_dict)
                interface = TrapiInterface(query=query)
                interface.build_chp_queries()
                interface.run_chp_queries()
                response = interface.construct_trapi_response()


    def test_wildcard_batch_onehop_query(self):
        for trapi_version, queries in self.wildcard_batch_queries.items():
            for name, query_dict in queries.items():
                #if name != 'gene_to_disease_proxy_context':
                #    continue
                query = Query.load(trapi_version, None, query=query_dict)
                interface = TrapiInterface(query=query)
                interface.build_chp_queries()
                interface.run_chp_queries()
                response = interface.construct_trapi_response()

    def test_gene_to_gene_onehop_query(self):
        query = Query.load("1.1", None, query=self.gene_to_gene_query)
        interface = TrapiInterface(query=query)
        interface.build_chp_queries()
        interface.run_chp_queries()
        response = interface.construct_trapi_response()

class TestHandlerErrors(unittest.TestCase):

    def test_more_than_one_contribution(self):
        with open('query_samples/error_samples/test_more_than_one_contribution.pk', 'rb') as f_:
            query =  pickle.load(f_)
        query = Query.load("1.1", None, query=query)
        with self.assertRaises(TooManyContributionNodes) as context:
            interface = TrapiInterface(query=query)

    def test_more_than_one_disease(self):
        with open('query_samples/error_samples/test_more_than_one_disease.pk', 'rb') as f_:
            query = pickle.load(f_)
        query = Query.load("1.1", None, query=query)
        with self.assertRaises(TooManyDiseaseNodes) as context:
            interface = TrapiInterface(query=query)

    def test_more_than_one_phenotype(self):
        with open('query_samples/error_samples/test_more_than_one_phenotype.pk', 'rb') as f_:
            query = pickle.load(f_)
        query = Query.load("1.1", None, query=query)
        with self.assertRaises(TooManyPhenotypeNodes) as context:
            interface = TrapiInterface(query=query)

    def test_more_than_one_disease(self):
        with open('query_samples/error_samples/test_more_than_one_disease.pk', 'rb') as f_:
            query = pickle.load(f_)
        query = Query.load("1.1", None, query=query)
        with self.assertRaises(TooManyDiseaseNodes) as context:
            interface = TrapiInterface(query=query)

    def test_no_disease(self):
        with open('query_samples/error_samples/test_no_disease.pk', 'rb') as f_:
            query = pickle.load(f_)
        query = Query.load("1.1", None, query=query)
        with self.assertRaises(UnidentifiedQueryType) as context:
            interface = TrapiInterface(query=query)

    def test_no_target(self):
        with open('query_samples/error_samples/test_no_target.pk', 'rb') as f_:
            query = pickle.load(f_)
        query = Query.load("1.1", None, query=query)
        with self.assertRaises(UnidentifiedQueryType) as context:
            interface = TrapiInterface(query=query)

    def test_illegal_drug_to_disease_default(self):
        with open('query_samples/error_samples/test_illegal_drug_to_disease_default.pk', 'rb') as f_:
            query = pickle.load(f_)
        query = Query.load("1.1", None, query=query)
        with self.assertRaises(MalformedSubjectObjectOnDrugToDisease) as context:
            interface = TrapiInterface(query=query)

    def test_illegal_gene_to_disease_default(self):
        with open('query_samples/error_samples/test_illegal_gene_to_disease_default.pk', 'rb') as f_:
            query = pickle.load(f_)
        query = Query.load("1.1", None, query=query)
        with self.assertRaises(MalformedSubjectObjectOnGeneToDisease) as context:
            interface = TrapiInterface(query=query)

    def test_illegal_disease_to_phenotype_default(self):
        with open('query_samples/error_samples/test_illegal_disease_to_phenotype_default.pk', 'rb') as f_:
            query = pickle.load(f_)
        query = Query.load("1.1", None, query=query)
        with self.assertRaises(MalformedSubjectObjectOnDiseaseToPhenotype) as context:
            interface = TrapiInterface(query=query)

    def test_illegal_edge_default(self):
        with open('query_samples/error_samples/test_illegal_edge_default.pk', 'rb') as f_:
            query = pickle.load(f_)
        query = Query.load("1.1", None, query=query)
        with self.assertRaises(IncompatibleDefaultEdge) as context:
            interface = TrapiInterface(query=query)

    def test_illegal_gene_to_disease_wildcard(self):
        with open('query_samples/error_samples/test_illegal_gene_to_disease_wildcard.pk', 'rb') as f_:
            query = pickle.load(f_)
        query = Query.load("1.1", None, query=query)
        with self.assertRaises(MalformedSubjectObjectOnGeneToDisease) as context:
            interface = TrapiInterface(query=query)

    def test_illegal_drug_to_disease_wildcard(self):
        with open('query_samples/error_samples/test_illegal_drug_to_disease_wildcard.pk', 'rb') as f_:
            query = pickle.load(f_)
        query = Query.load("1.1", None, query=query)
        with self.assertRaises(MalformedSubjectObjectOnDrugToDisease) as context:
            interface = TrapiInterface(query=query)

    def test_illegal_disease_to_phenotype_wildcard(self):
        with open('query_samples/error_samples/test_illegal_disease_to_phenotype_wildcard.pk', 'rb') as f_:
            query = pickle.load(f_)
        query = Query.load("1.1", None, query=query)
        with self.assertRaises(MalformedSubjectObjectOnDiseaseToPhenotype) as context:
            interface = TrapiInterface(query=query)

    def test_illegal_edge_wildcard(self):
        with open('query_samples/error_samples/test_illegal_edge_wildcard.pk', 'rb') as f_:
            query = pickle.load(f_)
        query = Query.load("1.1", None, query=query)
        with self.assertRaises(IncompatibleWildcardEdge) as context:
            interface = TrapiInterface(query=query)

    def test_illegal_drug_to_disease_one_hop(self):
        with open('query_samples/error_samples/test_illegal_drug_to_disease_one_hop.pk', 'rb') as f_:
            query = pickle.load(f_)
        query = Query.load("1.1", None, query=query)
        with self.assertRaises(MalformedSubjectObjectOnDrugToDisease) as context:
            interface = TrapiInterface(query=query)

    def test_illegal_gene_to_disease_one_hop(self):
        with open('query_samples/error_samples/test_illegal_gene_to_disease_one_hop.pk', 'rb') as f_:
            query = pickle.load(f_)
        query = Query.load("1.1", None, query=query)
        with self.assertRaises(MalformedSubjectObjectOnGeneToDisease) as context:
            interface = TrapiInterface(query=query)

    def test_illegal_disease_to_phenotype_one_hop(self):
        with open('query_samples/error_samples/test_illegal_disease_to_phenotype_one_hop.pk', 'rb') as f_:
            query = pickle.load(f_)
        query = Query.load("1.1", None, query=query)
        with self.assertRaises(UnexpectedEdgeType) as context:
            interface = TrapiInterface(query=query)

if __name__ == '__main__':
    unittest.main()
