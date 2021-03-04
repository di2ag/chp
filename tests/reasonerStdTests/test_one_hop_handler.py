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
        4. test_illegal_gene_to_chemical - backwards gene to chemical edge.
        5. test_illegal_chemical_to_disease - backwards chemical to disease edge.
        6. test_illegal_disease_to_phenotypic_feature - backwards disease to phenotype edge.
        7. test_drug_backwards - backwards subject object on curie node
        8. test_illegal_edge - unknown edge.
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

    def test_illegal_gene_to_chemical(self):
        """ test illegal gene to chemical edge. Should return appropriate error message
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
                                                       'predicate':BIOLINK_GENE_TO_DISEASE_PREDICATE,
                                                       'subject': 'n1',
                                                       'object': 'n0'
                                                     }

        # test input is TRAPI compliant
        validate_Message(reasoner_std) # doesn't return True/False for some reason... Will just present exception if not compliant

        handler = TrapiInterface(query=reasoner_std,
                                 max_results=2)
        queries = handler.build_chp_queries()
        self.assertEqual('Gene-chemical onehop detected. Received edge between gene and disease (edge id: e0). This edge is incompatible with this wildcard type.', handler.error_msg)

    def test_illegal_chemical_to_disease(self):
        """ test illegal chemical to disease edge. Should return appropriate error message
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
                                                       'predicate':BIOLINK_CHEMICAL_TO_DISEASE_OR_PHENOTYPIC_FEATURE_PREDICATE,
                                                       'subject': 'n1',
                                                       'object': 'n0'
                                                     }

        # test input is TRAPI compliant
        validate_Message(reasoner_std) # doesn't return True/False for some reason... Will just present exception if not compliant

        handler = TrapiInterface(query=reasoner_std,
                                 max_results=2)
        queries = handler.build_chp_queries()
        self.assertEqual('Gene-chemical onehop detected. Received edge between chemical and disease (edge id: e0). This edge is incompatible with this wildcard type.', handler.error_msg)

    def test_illegal_disease_to_phenotypic_feature(self):
        """ test illegal disease to phenotype. Should return appropriate error message
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
                                                       'predicate':BIOLINK_DISEASE_TO_PHENOTYPIC_FEATURE_PREDICATE,
                                                       'subject': 'n1',
                                                       'object': 'n0'
                                                     }

        # test input is TRAPI compliant
        validate_Message(reasoner_std) # doesn't return True/False for some reason... Will just present exception if not compliant

        handler = TrapiInterface(query=reasoner_std,
                                 max_results=2)
        queries = handler.build_chp_queries()
        self.assertEqual('Gene-chemical onehop detected. Received edge between disease and phenotype (edge id: e0). This edge is incompatible with this wildcard type.', handler.error_msg)

    def test_drug_backwards(self):

        """ test_drug_backwards. Should return appropriate error message
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
                                                       'subject': 'n0',
                                                       'object': 'n1'
                                                     }

        # test input is TRAPI compliant
        validate_Message(reasoner_std) # doesn't return True/False for some reason... Will just present exception if not compliant

        handler = TrapiInterface(query=reasoner_std,
                                 max_results=2)
        queries = handler.build_chp_queries()
        self.assertEqual('Gene-chemical onehop detected. Gene to chemical edge (edge id: e0) malformed.', handler.error_msg)

    def test_illegal_edge(self):

        """ test illegal edge. Should return appropriate error message
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
                                                       'predicate':'biolink:unknown_edge',
                                                       'subject': 'n0',
                                                       'object': 'n1'
                                                     }

        # test input is TRAPI compliant
        validate_Message(reasoner_std) # doesn't return True/False for some reason... Will just present exception if not compliant

        handler = TrapiInterface(query=reasoner_std,
                                 max_results=2)
        queries = handler.build_chp_queries()
        self.assertEqual('Gene-chemical onehop detected. Unknown predicate type for edge (edge id: e0).', handler.error_msg)

if __name__ == '__main__':
    unittest.main()
