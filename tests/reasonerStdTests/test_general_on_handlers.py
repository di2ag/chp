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

    """ Unit tests for the wildcard_handler
        1. test_more_than_one_contribution_node
        2. test_more_than_one_disease_node
        3. test_more_than_one_phenotype_node
    """

    def test_more_than_one_contribution_node(self):
        """ two nodes don't have a curie. Should return appropriate error message
        """

        # empty response
        reasoner_std = { "query_graph": {}
                       }
        # empty query graph
        reasoner_std["query_graph"] = { "edges": {},
                                        "nodes": {}
                                      }

        # add in evidence drug
        drug = ('CYCLOPHOSPHAMIDE', 'CHEMBL:CHEMBL88')
        reasoner_std['query_graph']['nodes']['n0'] = {
                                                      'category':BIOLINK_DRUG,
                                                      'id':'{}'.format(drug[1])
                                                     }

        # add in gene node (to be filled by contribution analysis
        reasoner_std['query_graph']['nodes']['n1'] = {
                                                      'category':BIOLINK_GENE
                                                     }

        reasoner_std['query_graph']['nodes']['n4'] = {
                                                      'category':BIOLINK_DRUG
                                                     }

        #add in disease node
        disease = ('Breast_Cancer', 'MONDO:0007254')
        reasoner_std['query_graph']['nodes']['n2'] = {
                                                      'category':BIOLINK_DISEASE,
                                                      'id':'{}'.format(disease[1])
                                                     }

        # add target survival node
        phenotype = ('Survival_Time', 'EFO:0000714')
        reasoner_std['query_graph']['nodes']['n3'] = {
                                                      'category': BIOLINK_PHENOTYPIC_FEATURE,
                                                      'id': '{}'.format(phenotype[1]),
                                                     }

        # link disease to target survival node
        reasoner_std['query_graph']['edges']['e0'] = {
                                                      'predicate':BIOLINK_DISEASE_TO_PHENOTYPIC_FEATURE_PREDICATE,
                                                      'subject':'n2',
                                                      'object':'n3',
                                                      'properties': {
                                                                     'qualifier': '>=',
                                                                     'days': 1000
                                                                    }
                                                     }

        # link genes/drugs to disease
        reasoner_std['query_graph']['edges']['e1'] = {
                                                      'predicate':BIOLINK_GENE_TO_DISEASE_PREDICATE,
                                                      'subject': 'n1',
                                                      'object': 'n2'
                                                     }
        reasoner_std['query_graph']['edges']['e2'] = {
                                                      'predicate':BIOLINK_CHEMICAL_TO_DISEASE_OR_PHENOTYPIC_FEATURE_PREDICATE,
                                                      'subject': 'n0',
                                                      'object': 'n2'
                                                     }

        reasoner_std['query_graph']['edges']['e3'] = {
                                                      'predicate':BIOLINK_CHEMICAL_TO_DISEASE_OR_PHENOTYPIC_FEATURE_PREDICATE,
                                                      'subject': 'n4',
                                                      'object': 'n2'
                                                     }

        # test input is TRAPI compliant
        validate_Message(reasoner_std) # doesn't return True/False for some reason... Will just present exception if not compliant

        handler = TrapiInterface(query=reasoner_std,
                                 max_results=2)
        query = handler.build_chp_queries()
        self.assertEqual('Can only have 1 node for contributions.', handler.error_msg)

    def test_more_than_one_disease_node(self):
        """ two nodes are for disease. Currently only one is supported. Should return appropriate error message
        """

        # empty response
        reasoner_std = { "query_graph": {}
                       }
        # empty query graph
        reasoner_std["query_graph"] = { "edges": {},
                                        "nodes": {}
                                      }

        # add in evidence drug
        drug = ('CYCLOPHOSPHAMIDE', 'CHEMBL:CHEMBL88')
        reasoner_std['query_graph']['nodes']['n0'] = {
                                                      'category':BIOLINK_DRUG,
                                                      'id':'{}'.format(drug[1])
                                                     }

        # add in gene node (to be filled by contribution analysis
        reasoner_std['query_graph']['nodes']['n1'] = {
                                                      'category':BIOLINK_GENE
                                                     }

        #add in disease node
        disease = ('Breast_Cancer', 'MONDO:0007254')
        reasoner_std['query_graph']['nodes']['n2'] = {
                                                      'category':BIOLINK_DISEASE,
                                                      'id':'{}'.format(disease[1])
                                                     }

        # add target survival node
        phenotype = ('Survival_Time', 'EFO:0000714')
        reasoner_std['query_graph']['nodes']['n3'] = {
                                                      'category': BIOLINK_PHENOTYPIC_FEATURE,
                                                      'id': '{}'.format(phenotype[1]),
                                                     }

        #add in disease node
        disease = ('Breast_Cancer', 'MONDO:0007254')
        reasoner_std['query_graph']['nodes']['n4'] = {
                                                      'category':BIOLINK_DISEASE,
                                                      'id':'{}'.format(disease[1])
                                                     }

        # link disease to target survival node
        reasoner_std['query_graph']['edges']['e0'] = {
                                                      'predicate':BIOLINK_DISEASE_TO_PHENOTYPIC_FEATURE_PREDICATE,
                                                      'subject':'n2',
                                                      'object':'n3',
                                                      'properties': {
                                                                     'qualifier': '>=',
                                                                     'days': 1000
                                                                    }
                                                     }

        # link genes/drugs to disease
        reasoner_std['query_graph']['edges']['e1'] = {
                                                      'predicate':BIOLINK_GENE_TO_DISEASE_PREDICATE,
                                                      'subject': 'n1',
                                                      'object': 'n2'
                                                     }
        reasoner_std['query_graph']['edges']['e2'] = {
                                                      'predicate':BIOLINK_CHEMICAL_TO_DISEASE_OR_PHENOTYPIC_FEATURE_PREDICATE,
                                                      'subject': 'n0',
                                                      'object': 'n2'
                                                     }

        # test input is TRAPI compliant
        validate_Message(reasoner_std) # doesn't return True/False for some reason... Will just present exception if not compliant

        handler = TrapiInterface(query=reasoner_std,
                                 max_results=2)
        query = handler.build_chp_queries()
        self.assertEqual('Can only have 1 node for disease.', handler.error_msg)

    def test_more_than_one_phenotype_node(self):
        """ two nodes are for phenotype. Only one phenotype is currently handled. Should return appropriate error message.
        """

        # empty response
        reasoner_std = { "query_graph": {}
                       }
        # empty query graph
        reasoner_std["query_graph"] = { "edges": {},
                                        "nodes": {}
                                      }

        # add in evidence drug
        drug = ('CYCLOPHOSPHAMIDE', 'CHEMBL:CHEMBL88')
        reasoner_std['query_graph']['nodes']['n0'] = {
                                                      'category':BIOLINK_DRUG,
                                                      'id':'{}'.format(drug[1])
                                                     }

        # add in gene node (to be filled by contribution analysis
        reasoner_std['query_graph']['nodes']['n1'] = {
                                                      'category':BIOLINK_GENE
                                                     }

        #add in disease node
        disease = ('Breast_Cancer', 'MONDO:0007254')
        reasoner_std['query_graph']['nodes']['n2'] = {
                                                      'category':BIOLINK_DISEASE,
                                                      'id':'{}'.format(disease[1])
                                                     }

        # add target survival node
        phenotype = ('Survival_Time', 'EFO:0000714')
        reasoner_std['query_graph']['nodes']['n3'] = {
                                                      'category': BIOLINK_PHENOTYPIC_FEATURE,
                                                      'id': '{}'.format(phenotype[1]),
                                                     }

        # add target survival node
        phenotype = ('Survival_Time', 'EFO:0000714')
        reasoner_std['query_graph']['nodes']['n4'] = {
                                                      'category': BIOLINK_PHENOTYPIC_FEATURE,
                                                      'id': '{}'.format(phenotype[1]),
                                                     }

        # link disease to target survival node
        reasoner_std['query_graph']['edges']['e0'] = {
                                                      'predicate':BIOLINK_DISEASE_TO_PHENOTYPIC_FEATURE_PREDICATE,
                                                      'subject':'n2',
                                                      'object':'n3',
                                                      'properties': {
                                                                     'qualifier': '>=',
                                                                     'days': 1000
                                                                    }
                                                     }

        # link genes/drugs to disease
        reasoner_std['query_graph']['edges']['e1'] = {
                                                      'predicate':BIOLINK_GENE_TO_DISEASE_PREDICATE,
                                                      'subject': 'n1',
                                                      'object': 'n2'
                                                     }
        reasoner_std['query_graph']['edges']['e2'] = {
                                                      'predicate':BIOLINK_CHEMICAL_TO_DISEASE_OR_PHENOTYPIC_FEATURE_PREDICATE,
                                                      'subject': 'n0',
                                                      'object': 'n2'
                                                     }
        # test input is TRAPI compliant
        validate_Message(reasoner_std) # doesn't return True/False for some reason... Will just present exception if not compliant

        handler = TrapiInterface(query=reasoner_std,
                                 max_results=2)
        query = handler.build_chp_queries()
        self.assertEqual('Can only have up to a single node for phenotype.', handler.error_msg)

if __name__ == '__main__':
    unittest.main()
