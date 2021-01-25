"""
    Source code developed by DI2AG.
    Thayer School of Engineering at Dartmouth College
    Authors:    Dr. Eugene Santos, Jr
                Mr. Chase Yakaboski,
                Mr. Gregory Hyde,
                Mr. Luke Veenhuis,
                Dr. Keum Joo Kim
"""

import requests
import json
import unittest

from chp_data.trapi_constants import *

from chp.trapi_interface import TrapiInterface
from reasoner_validator import validate_Message

class testWildCardHandler(unittest.TestCase):

    """ Unit tests for the one_hop_handler.
        1. test_drug - queries Cyclophosphamide which has no sparsity issues
        2. test_sparse_drug - queries the drug Gemzar which has a sparsity issue. 
               Should be handled be interpolation now
        3. test_gene - queries TP53 (edge should swap nodes here)
    """

    def test_drug(self):

        """ queries Cyclophosphamide which has no sparsity issue. Should return
            top 2 genes
        """

        # empty response
        reasoner_std = { "query_graph": dict()}
        # empty query graph
        reasoner_std["query_graph"] = { "edges": {},
                                        "nodes": {}
                                      }
        # add in evidence drug
        drug = ('CYCLOPHOSPHAMIDE', 'CHEMBL:CHEMBL88')
        reasoner_std['query_graph']['nodes']['n{}'.format('0')] = {
                                                      'category':BIOLINK_DRUG,
                                                      'id':'{}'.format(drug[1])
        }

        # add in gene node (to be filled by contribution analysis
        reasoner_std['query_graph']['nodes']['n{}'.format('1')] = {
                                                      'category':BIOLINK_GENE,
                                                   }

        # link genes/drugs to disease
        reasoner_std['query_graph']['edges']['e{}'.format(0)] = {
                                                       'predicate':BIOLINK_CHEMICAL_TO_GENE_PREDICATE,
                                                       'subject': 'n1',
                                                       'object': 'n0'
                                                     }

        # test input is TRAPI compliant
        validate_Message(reasoner_std) # doesn't return True/False for some reason... Will just present exception if not compliant

        handler = TrapiInterface(query=reasoner_std,
                                 max_results=2)
        queries = handler.build_chp_queries()
        queries = handler.run_chp_queries()
        reasoner_std_final = handler.construct_trapi_response()

        # test output is TRAPI compliant
        validate_Message(reasoner_std_final['message'])

        KG = reasoner_std_final["message"]['knowledge_graph']
        res = reasoner_std_final["message"]['results']

        res_pretty = json.dumps(reasoner_std_final, indent=2)
        print(res_pretty)


    def test_sparse_drug(self):

        """ queries the drug Gemzar which has a sparsity issue. Used to throw an error.
            Should be handled.
        """

        # empty response
        reasoner_std = { "query_graph": {}}
        # empty query graph
        reasoner_std["query_graph"] = { "edges": {},
                                        "nodes": {}
                                      }

        # add in evidence drug
        drug = ('GEMZAR', 'CHEMBL:CHEMBL888')
        reasoner_std['query_graph']['nodes']['n{}'.format('0')] = {
                                                      'category':BIOLINK_DRUG,
                                                      'id':'{}'.format(drug[1])
        }

        # add in gene node (to be filled by contribution analysis
        reasoner_std['query_graph']['nodes']['n{}'.format('1')] = {
                                                      'category':BIOLINK_GENE,
                                                   }

        # link genes/drugs to disease
        reasoner_std['query_graph']['edges']['e{}'.format(0)] = {
                                                       'predicate':BIOLINK_CHEMICAL_TO_GENE_PREDICATE,
                                                       'subject': 'n1',
                                                       'object': 'n0'
                                                     }

        # test input is TRAPI compliant
        validate_Message(reasoner_std) # doesn't return True/False for some reason... Will just present exception if not compliant

        handler = TrapiInterface(query=reasoner_std,
                                 max_results=2)
        queries = handler.build_chp_queries()
        queries = handler.run_chp_queries()
        reasoner_std_final = handler.construct_trapi_response()

        # test output is TRAPI compliant
        validate_Message(reasoner_std_final['message'])

        KG = reasoner_std_final["message"]['knowledge_graph']
        res = reasoner_std_final["message"]['results']

        res_pretty = json.dumps(reasoner_std_final, indent=2)
        print(res_pretty)

    def test_gene(self):

        """ queries TP53. Should return 2 drugs.
        """

        # empty response
        reasoner_std = { "query_graph": dict()}
        # empty query graph
        reasoner_std["query_graph"] = { "edges": {},
                                        "nodes": {}
                                      }
        # add in evidence drug
        reasoner_std['query_graph']['nodes']['n{}'.format('0')] = {
                                                      'category':BIOLINK_DRUG
                                                   }

        # add in gene node (to be filled by contribution analysis
        reasoner_std['query_graph']['nodes']['n{}'.format('1')] = {
                                                      'category':BIOLINK_GENE,
                                                      'id':'ENSEMBL:ENSG00000141510',
                                                   }

        # link genes/drugs to disease
        reasoner_std['query_graph']['edges']['e{}'.format(0)] = {
                                                       'predicate':BIOLINK_CHEMICAL_TO_GENE_PREDICATE,
                                                       'subject': 'n0',
                                                       'object': 'n1'
                                                     }

        # test input is TRAPI compliant
        validate_Message(reasoner_std) # doesn't return True/False for some reason... Will just present exception if not compliant

        handler = TrapiInterface(query=reasoner_std,
                                 max_results=2)
        queries = handler.build_chp_queries()
        queries = handler.run_chp_queries()
        reasoner_std_final = handler.construct_trapi_response()

        # test output is TRAPI compliant
        validate_Message(reasoner_std_final['message'])

        KG = reasoner_std_final["message"]['knowledge_graph']
        res = reasoner_std_final["message"]['results']

        res_pretty = json.dumps(reasoner_std_final, indent=2)
        print(res_pretty)

if __name__ == '__main__':
    unittest.main()
