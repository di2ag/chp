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

    """ Unit tests for the wildcard_handler
        1. test_drug - queries Cyclophosphamide which has no sparsity issues
        2. test_sparse_drug - queries the drug Gemzar which has a sparsity issue.
               Used to throw an error. Should be handled.
    """

    def test_drug(self):

        """ queries Cyclophosphamide which has no sparsity issue. Should return
            top 5 genes and a probability of survival
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
        reasoner_std['query_graph']['nodes']['n0'] = {
                                                      'category':'biolink:Drug',
                                                      'id':'{}'.format(drug[1])
                                                     }

        # add in gene node (to be filled by contribution analysis
        reasoner_std['query_graph']['nodes']['n1'] = {
                                                      'category':'biolink:Gene'
                                                     }

        #add in disease node
        disease = ('Breast_Cancer', 'MONDO:0007254')
        reasoner_std['query_graph']['nodes']['n2'] = {
                                                      'category':'biolink:Disease',
                                                      'id':'{}'.format(disease[1])
                                                     }

        # add target survival node
        phenotype = ('Survival_Time', 'EFO:0000714')
        reasoner_std['query_graph']['nodes']['n3'] = {
                                                      'category': 'biolink:PhenotypicFeature',
                                                      'id': '{}'.format(phenotype[1]),
                                                     }

        # link disease to target survival node
        reasoner_std['query_graph']['edges']['e0'] = {
                                                      'predicate':'biolink:DiseaseToPhenotypicFeatureAssociation',
                                                      'subject':'n2',
                                                      'object':'n3',
                                                      'properties': {
                                                                     'qualifier': '>=',
                                                                     'days': 940
                                                                    }
                                                     }

        # link genes/drugs to disease
        reasoner_std['query_graph']['edges']['e1'] = {
                                                      'predicate':'biolink:GeneToDiseaseAssociation',
                                                      'subject': 'n1',
                                                      'object': 'n2'
                                                     }
        reasoner_std['query_graph']['edges']['e2'] = {
                                                      'predicate':'biolink:ChemicalToDiseaseOrPhenotypicFeatureAssociation',
                                                      'subject': 'n0',
                                                      'object': 'n2'
                                                     }

        handler = ReasonerStdHandler(source_ara='default',
                                     dict_query=reasoner_std,
                                     max_results=5)
        queries = handler.buildChpQueries()
        queries = handler.runChpQueries()
        reasoner_std_final = handler.constructDecoratedKG()

        KG = reasoner_std_final["message"]['knowledge_graph']

        res_pretty = json.dumps(reasoner_std_final, indent=2)
        print(res_pretty)

        # extract probability
        for _, edge in KG['edges'].items():
            if edge['predicate'] == 'biolink:DiseaseToPhenotypicFeatureAssociation':
                p_survival = edge['has_confidence_level']
                #gene_contribs = edge['Description']
                break
        print("probability of survival:",p_survival)

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
        reasoner_std['query_graph']['nodes']['n0'] = {
                                                      'category':'biolink:Drug',
                                                      'id':'{}'.format(drug[1])
                                                     }

        # add in gene node (to be filled by contribution analysis
        reasoner_std['query_graph']['nodes']['n1'] = {
                                                      'category':'biolink:Gene'
                                                     }

        #add in disease node
        disease = ('Breast_Cancer', 'MONDO:0007254')
        reasoner_std['query_graph']['nodes']['n2'] = {
                                                      'category':'biolink:Disease',
                                                      'id':'{}'.format(disease[1])
                                                     }

        # add target survival node
        phenotype = ('Survival_Time', 'EFO:0000714')
        reasoner_std['query_graph']['nodes']['n3'] = {
                                                      'category': 'biolink:PhenotypicFeature',
                                                      'id': '{}'.format(phenotype[1]),
                                                     }

        # link disease to target survival node
        reasoner_std['query_graph']['edges']['e0'] = {
                                                      'predicate':'biolink:DiseaseToPhenotypicFeatureAssociation',
                                                      'subject':'n2',
                                                      'object':'n3',
                                                      'properties': {
                                                                     'qualifier': '>=',
                                                                     'days': 940
                                                                    }
                                                     }

        # link genes/drugs to disease
        reasoner_std['query_graph']['edges']['e1'] = {
                                                      'predicate':'biolink:GeneToDiseaseAssociation',
                                                      'subject': 'n1',
                                                      'object': 'n2'
                                                     }
        reasoner_std['query_graph']['edges']['e2'] = {
                                                      'predicate':'biolink:ChemicalToDiseaseOrPhenotypicFeatureAssociation',
                                                      'subject': 'n0',
                                                      'object': 'n2'
                                                     }

        handler = ReasonerStdHandler(source_ara='default',
                                     dict_query=reasoner_std,
                                     max_results=5)
        queries = handler.buildChpQueries()
        queries = handler.runChpQueries()
        reasoner_std_final = handler.constructDecoratedKG()

        KG = reasoner_std_final["message"]['knowledge_graph']

        res_pretty = json.dumps(reasoner_std_final, indent=2)
        print(res_pretty)

        # extract probability
        for _, edge in KG['edges'].items():
            if edge['predicate'] == 'biolink:DiseaseToPhenotypicFeatureAssociation':
                p_survival = edge['has_confidence_level']
                #gene_contribs = edge['Description']
                break
        print("probability of survival:",p_survival)


if __name__ == '__main__':
    unittest.main()
