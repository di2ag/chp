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
from trapi_model.processing_and_validation.metakg_validation_exceptions import *
UnsupportedNodeEdgeRelationship


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
    @classmethod
    def setUpClass(cls):
        super(TestDefaultHandler, cls).setUpClass()
        # load in sample query graphs
        with open('query_samples/random_queries.json', 'r') as f_:
            cls.queries = json.load(f_)
        with open('query_samples/random_batch_queries.json', 'r') as f_:
            cls.batch_queries = json.load(f_)

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
        self.assertIsInstance(curies, dict)

    def test_predicates(self):
        interface = TrapiInterface(
                bkb_handler=self.bkb_handler,
                dynamic_reasoner=self.dynamic_reasoner,
                joint_reasoner=self.joint_reasoner,
                )
        predicates = interface.get_predicates()
        self.assertIsInstance(predicates, dict)

    def test_meta_knowledge_graph(self):
        interface = TrapiInterface(
                bkb_handler=self.bkb_handler,
                dynamic_reasoner=self.dynamic_reasoner,
                joint_reasoner=self.joint_reasoner,
                )
        meta_kg = interface.get_meta_knowledge_graph()
        self.assertIsInstance(meta_kg, dict)

    def test_default_single_query(self):
        # This is a non-simple query
        logger.info('Running single query test.')
        for trapi_version, queries in self.queries.items():
            query_dict = copy.deepcopy(queries[0])
            query = Query.load(trapi_version, None, query=query_dict)
            interface = TrapiInterface(
                    query=query,
                    bkb_handler=self.bkb_handler,
                    dynamic_reasoner=self.dynamic_reasoner,
                    joint_reasoner=self.joint_reasoner,
                    )
            interface.build_chp_queries()
            interface.run_chp_queries()
            response = interface.construct_trapi_response()
    
    def test_inverse_query(self):
        # This is a simple query 
        logger.info('Running default inverse query test.')
        for trapi_version, queries in self.queries.items():
            query_dict = copy.deepcopy(queries[0])
            query = Query.load(trapi_version, None, query=query_dict)
            for edge_id, edge in query.message.query_graph.edges.items():
                predicate = edge.predicates[0]
                inverse = edge.predicates[0].get_inverse()
                edge.set_predicates(inverse)
                # Switch subject and object
                edge_subject = copy.deepcopy(edge.subject)
                edge_object = copy.deepcopy(edge.object)
                edge.subject = edge_object
                edge.object = edge_subject
            interface = TrapiInterface(
                    query=query,
                    bkb_handler=self.bkb_handler,
                    dynamic_reasoner=self.dynamic_reasoner,
                    joint_reasoner=self.joint_reasoner,
                    )
            interface.build_chp_queries()
            interface.run_chp_queries()
            response = interface.construct_trapi_response()

    def test_simple_single_query(self):
        # This is a simple query 
        logger.info('Running single simple query test.')
        for trapi_version, queries in self.queries.items():
            query_dict = copy.deepcopy(queries[3])
            query = Query.load(trapi_version, None, query=query_dict)
            interface = TrapiInterface(
                    query=query,
                    bkb_handler=self.bkb_handler,
                    dynamic_reasoner=self.dynamic_reasoner,
                    joint_reasoner=self.joint_reasoner,
                    )
            interface.build_chp_queries()
            interface.run_chp_queries()
            response = interface.construct_trapi_response()
        
    def test_default_batch_query(self):
        logger.info('Running batch query test.')
        # These are non-simple queries
        for trapi_version, queries in self.batch_queries.items():
            query_dict = copy.deepcopy(queries[0])
            query = Query.load(trapi_version, None, query=query_dict)
            interface = TrapiInterface(
                    query=query,
                    bkb_handler=self.bkb_handler,
                    dynamic_reasoner=self.dynamic_reasoner,
                    joint_reasoner=self.joint_reasoner,
                    )
            interface.build_chp_queries()
            interface.run_chp_queries()
            response = interface.construct_trapi_response()

    def test_simple_batch_query(self):
        # These are simple queries
        logger.info('Running batch simple query test.')
        for trapi_version, queries in self.batch_queries.items():
            query_dict = copy.deepcopy(queries[1])
            query = Query.load(trapi_version, None, query=query_dict)
            interface = TrapiInterface(
                    query=query,
                    bkb_handler=self.bkb_handler,
                    dynamic_reasoner=self.dynamic_reasoner,
                    joint_reasoner=self.joint_reasoner,
                    )
            interface.build_chp_queries()
            interface.run_chp_queries()
            response = interface.construct_trapi_response()

class TestWildCardHandler(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        super(TestWildCardHandler, cls).setUpClass()
        with open('query_samples/random_gene_wildcard_queries.json', 'r') as f_:
            cls.gene_queries = json.load(f_)
        with open('query_samples/random_gene_wildcard_batch_queries.json', 'r') as f_:
            cls.gene_batch_queries = json.load(f_)
        with open('query_samples/random_drug_wildcard_queries.json', 'r') as f_:
            cls.drug_queries = json.load(f_)
        with open('query_samples/random_drug_wildcard_batch_queries.json', 'r') as f_:
            cls.drug_batch_queries = json.load(f_)
        cls.bkb_handler = BkbDataHandler()
        cls.dynamic_reasoner = ChpDynamicReasoner(cls.bkb_handler)
        cls.joint_reasoner = ChpJointReasoner(cls.bkb_handler)

    def test_single_gene_wildcard_query(self):
        for trapi_version, queries in self.gene_queries.items():
            _queries = copy.deepcopy(queries)
            query = Query.load(trapi_version, None, query=_queries[0])
            interface = TrapiInterface(
                    query=query,
                    bkb_handler=self.bkb_handler,
                    dynamic_reasoner=self.dynamic_reasoner,
                    joint_reasoner=self.joint_reasoner,
                    )
            interface.build_chp_queries()
            interface.run_chp_queries()
            response = interface.construct_trapi_response()

    def test_inverse_wildcard_query(self):
        for trapi_version, queries in self.gene_queries.items():
            _queries = copy.deepcopy(queries)
            query = Query.load(trapi_version, None, query=_queries[0])
            for edge_id, edge in query.message.query_graph.edges.items():
                predicate = edge.predicates[0]
                inverse = edge.predicates[0].get_inverse()
                edge.set_predicates(inverse)
                # Switch subject and object
                edge_subject = copy.deepcopy(edge.subject)
                edge_object = copy.deepcopy(edge.object)
                edge.subject = edge_object
                edge.object = edge_subject
            interface = TrapiInterface(
                    query=query,
                    bkb_handler=self.bkb_handler,
                    dynamic_reasoner=self.dynamic_reasoner,
                    joint_reasoner=self.joint_reasoner,
                    )
            interface.build_chp_queries()
            interface.run_chp_queries()
            response = interface.construct_trapi_response()


    def test_single_drug_wildcard_query(self):
        for trapi_version, queries in self.drug_queries.items():
            _queries = copy.deepcopy(queries)
            query = Query.load(trapi_version, None, query=_queries[0])
            interface = TrapiInterface(
                    query=query,
                    bkb_handler=self.bkb_handler,
                    dynamic_reasoner=self.dynamic_reasoner,
                    joint_reasoner=self.joint_reasoner,
                    )
            interface.build_chp_queries()
            interface.run_chp_queries()
            response = interface.construct_trapi_response()

    def test_batch_gene_wildcard_query(self):
        for trapi_version, queries in self.gene_batch_queries.items():
            _queries = copy.deepcopy(queries)
            query = Query.load(trapi_version, None, query=_queries[0])
            interface = TrapiInterface(
                    query=query,
                    bkb_handler=self.bkb_handler,
                    dynamic_reasoner=self.dynamic_reasoner,
                    joint_reasoner=self.joint_reasoner,
                    )
            interface.build_chp_queries()
            interface.run_chp_queries()
            response = interface.construct_trapi_response()
    
    def test_batch_drug_wildcard_query(self):
        for trapi_version, queries in self.drug_batch_queries.items():
            _queries = copy.deepcopy(queries)
            query = Query.load(trapi_version, None, query=_queries[0])
            interface = TrapiInterface(
                    query=query,
                    bkb_handler=self.bkb_handler,
                    dynamic_reasoner=self.dynamic_reasoner,
                    joint_reasoner=self.joint_reasoner,
                    )
            interface.build_chp_queries()
            interface.run_chp_queries()
            response = interface.construct_trapi_response()


class TestOneHopHandler(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        super(TestOneHopHandler, cls).setUpClass()
        # load in sample query graphs
        with open('query_samples/standard_single_onehop_queries.json', 'r') as f_:
            cls.standard_single_queries = json.load(f_)
        with open('query_samples/standard_batch_onehop_queries.json', 'r') as f_:
            cls.standard_batch_queries = json.load(f_)
        with open('query_samples/wildcard_single_onehop_queries.json', 'r') as f_:
            cls.wildcard_single_queries = json.load(f_)
        with open('query_samples/wildcard_batch_onehop_queries.json', 'r') as f_:
            cls.wildcard_batch_queries = json.load(f_)
        cls.bkb_handler = BkbDataHandler()
        cls.dynamic_reasoner = ChpDynamicReasoner(cls.bkb_handler)
        cls.joint_reasoner = ChpJointReasoner(cls.bkb_handler)

    def test_standard_single_onehop_query(self):
        for trapi_version, queries in self.standard_single_queries.items():
            for name, query_dict in queries.items():
                query = Query.load(trapi_version, None, query=query_dict[0])
                interface = TrapiInterface(
                        query=query,
                        bkb_handler=self.bkb_handler,
                        dynamic_reasoner=self.dynamic_reasoner,
                        joint_reasoner=self.joint_reasoner,
                        )
                interface.build_chp_queries()
                interface.run_chp_queries()
                response = interface.construct_trapi_response()
    
    def test_inverse_onehop_query(self):
        for trapi_version, queries in self.standard_single_queries.items():
            for name, query_dict in queries.items():
                #if name != 'gene_to_disease_proxy_context':
                #    continue
                query = Query.load(trapi_version, None, query=query_dict[0])
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
                interface = TrapiInterface(
                        query=query,
                        bkb_handler=self.bkb_handler,
                        dynamic_reasoner=self.dynamic_reasoner,
                        joint_reasoner=self.joint_reasoner,
                        )
                interface.build_chp_queries()
                interface.run_chp_queries()
                response = interface.construct_trapi_response()
            
    def test_standard_batch_onehop_query(self):
        for trapi_version, queries in self.standard_batch_queries.items():
            for name, query_dict in queries.items():
                #if name != 'gene_to_disease_proxy_context':
                #    continue
                query = Query.load(trapi_version, None, query=query_dict[0])
                interface = TrapiInterface(
                        query=query,
                        bkb_handler=self.bkb_handler,
                        dynamic_reasoner=self.dynamic_reasoner,
                        joint_reasoner=self.joint_reasoner,
                        )
                interface.build_chp_queries()
                interface.run_chp_queries()
                response = interface.construct_trapi_response()
    
    def test_wildcard_single_onehop_query(self):
        for trapi_version, queries in self.wildcard_single_queries.items():
            for name, query_dict in queries.items():
                #if name != 'gene_to_disease_proxy_context':
                #    continue
                query = Query.load(trapi_version, None, query=query_dict)
                interface = TrapiInterface(
                        query=query,
                        bkb_handler=self.bkb_handler,
                        dynamic_reasoner=self.dynamic_reasoner,
                        joint_reasoner=self.joint_reasoner,
                        )
                interface.build_chp_queries()
                interface.run_chp_queries()
                response = interface.construct_trapi_response()
        
    def test_wildcard_batch_onehop_query(self):
        for trapi_version, queries in self.wildcard_batch_queries.items():
            for name, query_dict in queries.items():
                #if name != 'gene_to_disease_proxy_context':
                #    continue
                query = Query.load(trapi_version, None, query=query_dict[0])
                interface = TrapiInterface(
                        query=query,
                        bkb_handler=self.bkb_handler,
                        dynamic_reasoner=self.dynamic_reasoner,
                        joint_reasoner=self.joint_reasoner,
                        )
                interface.build_chp_queries()
                interface.run_chp_queries()
                response = interface.construct_trapi_response()

    def test_gene_to_gene_query(self):
        query = Query.load('1.1', None, query = self.gene_to_gene_query)
        interface = TrapiInterface(
                        query=query,
                        bkb_handler=self.bkb_handler,
                        dynamic_reasoner=self.dynamic_reasoner,
                        joint_reasoner=self.joint_reasoner,
                        )
        interface.build_chp_queries()
        interface.run_chp_queries()
        response = interface.construct_trapi_response()

class TestHandlerErrors(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        super(TestHandlerErrors, cls).setUpClass()
        cls.bkb_handler = BkbDataHandler()
        cls.dynamic_reasoner = ChpDynamicReasoner(cls.bkb_handler)
        cls.joint_reasoner = ChpJointReasoner(cls.bkb_handler)

    def test_more_than_one_contribution(self):
        with open('query_samples/error_samples/test_more_than_one_contribution.pk', 'rb') as f_:
            query = pickle.load(f_)
        with self.assertRaises(TooManyContributionNodes) as context:
            query = Query.load('1.1', None, query = query)
            interface = TrapiInterface(
                    query=query,
                    bkb_handler=self.bkb_handler,
                    dynamic_reasoner=self.dynamic_reasoner,
                    joint_reasoner=self.joint_reasoner,
                    )

    def test_more_than_one_disease(self):
        with open('query_samples/error_samples/test_more_than_one_disease.pk', 'rb') as f_:
            query = pickle.load(f_)
        with self.assertRaises(TooManyDiseaseNodes) as context:
            query = Query.load('1.1', None, query = query)
            interface = TrapiInterface(
                    query=query,
                    bkb_handler=self.bkb_handler,
                    dynamic_reasoner=self.dynamic_reasoner,
                    joint_reasoner=self.joint_reasoner,
                    )

    def test_more_than_one_phenotype(self):
        with open('query_samples/error_samples/test_more_than_one_phenotype.pk', 'rb') as f_:
            query = pickle.load(f_)
        with self.assertRaises(TooManyPhenotypeNodes) as context:
            query = Query.load('1.1', None, query = query)
            interface = TrapiInterface(
                    query=query,
                    bkb_handler=self.bkb_handler,
                    dynamic_reasoner=self.dynamic_reasoner,
                    joint_reasoner=self.joint_reasoner,
                    )

    def test_more_than_one_disease(self):
        with open('query_samples/error_samples/test_more_than_one_disease.pk', 'rb') as f_:
            query = pickle.load(f_)
        with self.assertRaises(TooManyDiseaseNodes) as context:
            query = Query.load('1.1', None, query = query)
            interface = TrapiInterface(
                    query=query,
                    bkb_handler=self.bkb_handler,
                    dynamic_reasoner=self.dynamic_reasoner,
                    joint_reasoner=self.joint_reasoner,
                    )

    def test_no_disease(self):
        with open('query_samples/error_samples/test_no_disease.pk', 'rb') as f_:
            query = pickle.load(f_)
        with self.assertRaises(UnidentifiedQueryType) as context:
            query = Query.load('1.1', None, query = query)
            interface = TrapiInterface(
                    query=query,
                    bkb_handler=self.bkb_handler,
                    dynamic_reasoner=self.dynamic_reasoner,
                    joint_reasoner=self.joint_reasoner,
                    )

    def test_no_target(self):
        with open('query_samples/error_samples/test_no_target.pk', 'rb') as f_:
            query = pickle.load(f_)
        with self.assertRaises(UnidentifiedQueryType) as context:
            query = Query.load('1.1', None, query = query)
            interface = TrapiInterface(
                    query=query,
                    bkb_handler=self.bkb_handler,
                    dynamic_reasoner=self.dynamic_reasoner,
                    joint_reasoner=self.joint_reasoner,
                    )

    def test_illegal_drug_to_disease_default(self):
        with open('query_samples/error_samples/test_illegal_drug_to_disease_default.pk', 'rb') as f_:
            query = pickle.load(f_)
            print(json.dumps(query, indent=2))
        with self.assertRaises(UnsupportedNodeEdgeRelationship) as context:
            query = Query.load('1.1', None, query = query)
            interface = TrapiInterface(
                    query=query,
                    bkb_handler=self.bkb_handler,
                    dynamic_reasoner=self.dynamic_reasoner,
                    joint_reasoner=self.joint_reasoner,
                    )

    def test_illegal_gene_to_disease_default(self):
        with open('query_samples/error_samples/test_illegal_gene_to_disease_default.pk', 'rb') as f_:
            query = pickle.load(f_)
        with self.assertRaises(UnsupportedNodeEdgeRelationship) as context:
            query = Query.load('1.1', None, query = query)
            interface = TrapiInterface(
                    query=query,
                    bkb_handler=self.bkb_handler,
                    dynamic_reasoner=self.dynamic_reasoner,
                    joint_reasoner=self.joint_reasoner,
                    )

    def test_illegal_disease_to_phenotype_default(self):
        with open('query_samples/error_samples/test_illegal_disease_to_phenotype_default.pk', 'rb') as f_:
            query = pickle.load(f_)
        with self.assertRaises(UnsupportedNodeEdgeRelationship) as context:
            query = Query.load('1.1', None, query = query)
            interface = TrapiInterface(
                    query=query,
                    bkb_handler=self.bkb_handler,
                    dynamic_reasoner=self.dynamic_reasoner,
                    joint_reasoner=self.joint_reasoner,
                    )

    def test_illegal_edge_default(self):
        with open('query_samples/error_samples/test_illegal_edge_default.pk', 'rb') as f_:
            query = pickle.load(f_)
        with self.assertRaises(UnsupportedNodeEdgeRelationship) as context:
            query = Query.load('1.1', None, query = query)
            interface = TrapiInterface(
                    query=query,
                    bkb_handler=self.bkb_handler,
                    dynamic_reasoner=self.dynamic_reasoner,
                    joint_reasoner=self.joint_reasoner,
                    )

    def test_illegal_gene_to_disease_wildcard(self):
        with open('query_samples/error_samples/test_illegal_gene_to_disease_wildcard.pk', 'rb') as f_:
            query = pickle.load(f_)
        with self.assertRaises(UnsupportedNodeEdgeRelationship) as context:
            query = Query.load('1.1', None, query = query)
            interface = TrapiInterface(
                    query=query,
                    bkb_handler=self.bkb_handler,
                    dynamic_reasoner=self.dynamic_reasoner,
                    joint_reasoner=self.joint_reasoner,
                    )

    def test_illegal_drug_to_disease_wildcard(self):
        with open('query_samples/error_samples/test_illegal_drug_to_disease_wildcard.pk', 'rb') as f_:
            query = pickle.load(f_)
        with self.assertRaises(UnsupportedNodeEdgeRelationship) as context:
            query = Query.load('1.1', None, query = query)
            interface = TrapiInterface(
                    query=query,
                    bkb_handler=self.bkb_handler,
                    dynamic_reasoner=self.dynamic_reasoner,
                    joint_reasoner=self.joint_reasoner,
                    )

    def test_illegal_disease_to_phenotype_wildcard(self):
        with open('query_samples/error_samples/test_illegal_disease_to_phenotype_wildcard.pk', 'rb') as f_:
            query = pickle.load(f_)
        with self.assertRaises(UnsupportedNodeEdgeRelationship) as context:
            query = Query.load('1.1', None, query = query)
            interface = TrapiInterface(
                    query=query,
                    bkb_handler=self.bkb_handler,
                    dynamic_reasoner=self.dynamic_reasoner,
                    joint_reasoner=self.joint_reasoner,
                    )

    def test_illegal_edge_wildcard(self):
        with open('query_samples/error_samples/test_illegal_edge_wildcard.pk', 'rb') as f_:
            query = pickle.load(f_)
        with self.assertRaises(IncompatibleWildcardEdge) as context:
            query = Query.load('1.1', None, query = query)
            interface = TrapiInterface(
                    query=query,
                    bkb_handler=self.bkb_handler,
                    dynamic_reasoner=self.dynamic_reasoner,
                    joint_reasoner=self.joint_reasoner,
                    )

    def test_illegal_drug_to_disease_one_hop(self):
        with open('query_samples/error_samples/test_illegal_drug_to_disease_one_hop.pk', 'rb') as f_:
            query = pickle.load(f_)
        with self.assertRaises(UnsupportedNodeEdgeRelationship) as context:
            query = Query.load('1.1', None, query = query)
            interface = TrapiInterface(
                    query=query,
                    bkb_handler=self.bkb_handler,
                    dynamic_reasoner=self.dynamic_reasoner,
                    joint_reasoner=self.joint_reasoner,
                    )

    def test_illegal_disease_to_phenotype_one_hop(self):
        with open('query_samples/error_samples/test_illegal_disease_to_phenotype_one_hop.pk', 'rb') as f_:
            query = pickle.load(f_)
        with self.assertRaises(UnsupportedNodeEdgeRelationship) as context:
            query = Query.load('1.1', None, query = query)
            interface = TrapiInterface(
                    query=query,
                    bkb_handler=self.bkb_handler,
                    dynamic_reasoner=self.dynamic_reasoner,
                    joint_reasoner=self.joint_reasoner,
                    )

if __name__ == '__main__':
    unittest.main()
