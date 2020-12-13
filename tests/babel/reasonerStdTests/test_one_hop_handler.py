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

from chp.reasoner_std import ReasonerStdHandler

class testWildCardHandler(unittest.TestCase):

    """ Unit tests for the one_hop_handler.
        1. test_drug - queries Cyclophosphamide which has no sparsity issues
        2. test_sparse_drug - queries the drug Gemzar which has a sparsity issue. 
               Used to throw an error. Should be handled.
    """

    def test_drug(self):

        """ queries Cyclophosphamide which has no sparsity issue. Should return
            top 5 genes
        """

        # empty response
        reasoner_std = { "query_graph": {},
                         "knowledge_graph": {},
                         "results": []
                       }
        # empty query graph
        reasoner_std["query_graph"] = { "edges": {},
                                        "nodes": {}
                                      }
        # empty knowledge graph
        reasoner_std["knowledge_graph"] = { "edges": {},
                                            "nodes": {}
                                          }

        # add in evidence drug
        drug = ('CYCLOPHOSPHAMIDE', 'CHEMBL:CHEMBL88')
        reasoner_std['query_graph']['nodes']['n{}'.format('0')] = {
                                                      'category':'biolink:Drug',
                                                      'id':'{}'.format(drug[1])
        }

        # add in gene node (to be filled by contribution analysis
        reasoner_std['query_graph']['nodes']['n{}'.format('1')] = {
                                                      'category':'biolink:Gene',
                                                   }

        # link genes/drugs to disease
        reasoner_std['query_graph']['edges']['e{}'.format(0)] = {
                                                       'predicate':'biolink:ChemicalToGeneAssociation',
                                                       'subject': 'n0',
                                                       'object': 'n1'
                                                     }

        handler = ReasonerStdHandler(source_ara='default',
                                     dict_query=reasoner_std,
                                     max_results=5)
        queries = handler.buildChpQueries()
        queries = handler.runChpQueries()
        reasoner_std_final = handler.constructDecoratedKG()

        KG = reasoner_std_final["message"]['knowledge_graph']
        res = reasoner_std_final["message"]['results']

        res_pretty = json.dumps(reasoner_std_final, indent=2)
        print(res_pretty)

    def test_sparse_drug(self):

        """ queries the drug Gemzar which has a sparsity issue. Used to throw an error.
            Should be handled.
        """

        # empty response
        reasoner_std = { "query_graph": {},
                         "knowledge_graph": {},
                         "results": []
                       }
        # empty query graph
        reasoner_std["query_graph"] = { "edges": {},
                                        "nodes": {}
                                      }
        # empty knowledge graph
        reasoner_std["knowledge_graph"] = { "edges": {},
                                            "nodes": {}
                                          }

        # add in evidence drug
        drug = ('GEMZAR', 'CHEMBL:CHEMBL888')
        reasoner_std['query_graph']['nodes']['n{}'.format('0')] = {
                                                      'category':'biolink:Drug',
                                                      'id':'{}'.format(drug[1])
        }

        # add in gene node (to be filled by contribution analysis
        reasoner_std['query_graph']['nodes']['n{}'.format('1')] = {
                                                      'category':'biolink:Gene',
                                                   }

        # link genes/drugs to disease
        reasoner_std['query_graph']['edges']['e{}'.format(0)] = {
                                                       'predicate':'biolink:chemical_to_gene_association',
                                                       'subject': 'n0',
                                                       'object': 'n1'
                                                     }

        handler = ReasonerStdHandler(source_ara='default',
                                     dict_query=reasoner_std,
                                     max_results=5)
        queries = handler.buildChpQueries()
        queries = handler.runChpQueries()
        reasoner_std_final = handler.constructDecoratedKG()

        KG = reasoner_std_final["message"]['knowledge_graph']
        res = reasoner_std_final["message"]['results']

        res_pretty = json.dumps(reasoner_std_final, indent=2)
        print(res_pretty)

if __name__ == '__main__':
    unittest.main()
