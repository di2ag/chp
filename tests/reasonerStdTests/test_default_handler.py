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

class testDefaultHandler(unittest.TestCase):

    """ Unit tests for default_handler. Uses 8 test cases:
        1. Testing negative result normalization
        2. Normal request with two genes and a drug
        3. Test with no evidence
        4. Test with no evidence, but omitting KG and Results (should be handled by handler)
        5. Test with one gene
        6. Test with one drug
        7. Test with no target (should crash)
        8. Test with no disease
        9. Test default survival
        10. check if the contribution property returns contributions
    """

    def test_negative(self):

        """ Testing negative result normalization
        """

        # empty response
        reasoner_std = { "query_graph": dict(),
                         "knowledge_graph": dict(),
                         "results": []
                       }
        # empty query graph
        reasoner_std["query_graph"] = { "edges": dict(),
                                        "nodes": dict()
                                      }
        # empty knowledge graph
        reasoner_std["knowledge_graph"] = { "edges": dict(),
                                            "nodes": dict()
                                          }
        # add in evidence gene
        gene1 = ('X', 'ENSEMBL:ENSG00000007174')
        reasoner_std['query_graph']['nodes']['n0'] = { 'category':'biolink:Gene',
                                                       'id':'{}'.format(gene1[1])
                                                     }
        # add in evidence drug
        drug = ('Y', 'CHEMBL:CHEMBL554')
        reasoner_std['query_graph']['nodes']['n1'] = { 'category':'biolink:Drug',
                                                       'id':'{}'.format(drug[1])
                                                     }
        # add in disease node
        disease = ('Breast_Cancer', 'MONDO:0007254')
        reasoner_std['query_graph']['nodes']['n2'] = { 'category':'biolink:Disease',
                                                       'id':'{}'.format(disease[1])
                                                     }
        # add target survival node
        phenotype = ('Survival_Time', 'EFO:0000714')
        reasoner_std['query_graph']['nodes']['n3'] = { 'category': 'biolink:PhenotypicFeature',
                                                       'id': '{}'.format(phenotype[1]),
                                                     }
        # link genes/drugs to disease
        reasoner_std['query_graph']['edges']['e0'] = { 'predicate':'biolink:GeneToDiseaseAssociation',
                                                       'subject': 'n0',
                                                       'object': 'n2'
                                                     }
        reasoner_std['query_graph']['edges']['e1'] = { 'predicate':'biolink:ChemicalToDiseaseOrPhenotypicFeatureAssociation',
                                                       'subject': 'n1',
                                                       'object': 'n2'
                                                     }
        # link disease to target
        reasoner_std['query_graph']['edges']['e2'] = { 'predicate':'biolink:DiseaseToPhenotypicFeatureAssociation',
                                                       'subject': 'n2',
                                                       'object': 'n3',
                                                       'properties': { 'qualifier':'>=',
                                                                       'days':970
                                                                     }
                                                     }
        handler = ReasonerStdHandler(source_ara='default',
                                     dict_query=reasoner_std)
        queries = handler.buildChpQueries()
        queries = handler.runChpQueries()
        reasoner_std_final = handler.constructDecoratedKG()
        KG = reasoner_std_final['message']['knowledge_graph']
        for edge_key in KG['edges'].keys():
            edge = KG['edges'][edge_key]
            if edge['predicate'] == 'biolink:DiseaseToPhenotypicFeatureAssociation':
                p_survival = edge['has_confidence_level']
        print("probability of survival:",p_survival)

    def test_normal_two_genes_and_drug(self):
        """ Normal request with two genes and a drug
        """

        # empty response
        reasoner_std = { "query_graph": dict(),
                         "knowledge_graph": dict(),
                         "results": list()
                       }
        # empty query graph
        reasoner_std["query_graph"] = { "edges": dict(),
                                        "nodes": dict()
                                      }
        # empty knowledge graph
        reasoner_std["knowledge_graph"] = { "edges": dict(),
                                            "nodes": dict()
                                          }

        # add in evidence gene
        gene1 = ('RAF1', 'ENSEMBL:ENSG00000132155')
        reasoner_std['query_graph']['nodes']['n0'] = { 'category':'biolink:Gene',
                                                       'id':'{}'.format(gene1[1])
                                                     }

        gene2 = ('BRCA1', 'ENSEMBL:ENSG00000012048')
        reasoner_std['query_graph']['nodes']['n1'] = { 'category':'biolink:Gene',
                                                       'id':'{}'.format(gene2[1])
                                                     }
        # add in evidence drug
        drug = ('CYCLOPHOSPHAMIDE', 'CHEMBL:CHEMBL88')
        reasoner_std['query_graph']['nodes']['n2'] = { 'category':'biolink:Drug',
                                                       'id':'{}'.format(drug[1])
                                                     }
        # add in disease node
        disease = ('Breast_Cancer', 'MONDO:0007254')
        reasoner_std['query_graph']['nodes']['n3'] = { 'category':'biolink:Disease',
                                                       'id':'{}'.format(disease[1])
                                                     }
        # add target survival node
        phenotype = ('Survival_Time', 'EFO:0000714')
        reasoner_std['query_graph']['nodes']['n4'] = { 'category': 'biolink:PhenotypicFeature',
                                                       'id': '{}'.format(phenotype[1]),
                                                     }
        # link genes/drugs to disease
        reasoner_std['query_graph']['edges']['e0'] = { 'predicate':'biolink:GeneToDiseaseAssociation',
                                                       'subject': 'n0',
                                                       'object': 'n3'
                                                     }
        reasoner_std['query_graph']['edges']['e1'] = { 'predicate':'biolink:GeneToDiseaseAssociation',
                                                       'subject': 'n1',
                                                       'object': 'n3'
                                                     }
        reasoner_std['query_graph']['edges']['e2'] = { 'predicate':'biolink:ChemicalToDiseaseOrPhenotypicFeatureAssociation',
                                                       'subject': 'n2',
                                                       'object': 'n3'
                                                     }
        # link disease to target
        reasoner_std['query_graph']['edges']['e3'] = { 'predicate':'biolink:DiseaseToPhenotypicFeatureAssociation',
                                                       'subject': 'n3',
                                                       'object': 'n4',
                                                       'properties': { 'qualifier':'>=',
                                                                       'days':970
                                                                     }
                                                     }

        handler = ReasonerStdHandler(source_ara='default',
                                     dict_query=reasoner_std)
        queries = handler.buildChpQueries()
        queries = handler.runChpQueries()
        reasoner_std_final = handler.constructDecoratedKG()
        KG = reasoner_std_final['message']['knowledge_graph']
        for edge_key in KG['edges'].keys():
            edge = KG['edges'][edge_key]
            if edge['predicate'] == 'biolink:DiseaseToPhenotypicFeatureAssociation':
                p_survival = edge['has_confidence_level']
        print("probability of survival:",p_survival)

    def test_no_evidence(self):
        """ Test with no evidence
        """
	# empty response
        reasoner_std = { "query_graph": dict(),
                         "knowledge_graph": dict(),
                         "results": list()
                       }
        # empty query graph
        reasoner_std["query_graph"] = { "edges": dict(),
                                        "nodes": dict()
                                      }
        # empty knowledge graph
        reasoner_std["knowledge_graph"] = { "edges": dict(),
                                            "nodes": dict()
                                          }

        # add in disease node
        disease = ('Breast_Cancer', 'MONDO:0007254')
        reasoner_std['query_graph']['nodes']['n0'] = { 'category':'biolink:Disease',
                                                       'id':'{}'.format(disease[1])
                                                     }
        # add target survival node
        phenotype = ('Survival_Time', 'EFO:0000714')
        reasoner_std['query_graph']['nodes']['n1'] = { 'category': 'biolink:PhenotypicFeature',
                                                       'id': '{}'.format(phenotype[1]),
                                                     }
        # link disease to target
        reasoner_std['query_graph']['edges']['e0'] = { 'predicate':'biolink:DiseaseToPhenotypicFeatureAssociation',
                                                       'subject': 'n0',
                                                       'object': 'n1',
                                                       'properties': { 'qualifier':'>=',
                                                                       'days':970
                                                                     }
                                                     }

        handler = ReasonerStdHandler(source_ara='default',
                                     dict_query=reasoner_std)
        queries = handler.buildChpQueries()
        queries = handler.runChpQueries()
        reasoner_std_final = handler.constructDecoratedKG()
        KG = reasoner_std_final['message']['knowledge_graph']
        for edge_key in KG['edges'].keys():
            edge = KG['edges'][edge_key]
            if edge['predicate'] == 'biolink:DiseaseToPhenotypicFeatureAssociation':
                p_survival = edge['has_confidence_level']
        print("probability of survival:",p_survival)

    def test_no_evidence_omitting_KG_and_results(self):
        """ Test with no evidence, but omitting KG and Results (should be handled by handler)
        """
        # empty response
        reasoner_std = { "query_graph": dict(),
                       }

        # empty query graph
        reasoner_std["query_graph"] = { "edges": dict(),
                                        "nodes": dict()
                                      }
        # add in disease node
        disease = ('Breast_Cancer', 'MONDO:0007254')
        reasoner_std['query_graph']['nodes']['n0'] = { 'category':'biolink:Disease',
                                                       'id':'{}'.format(disease[1])
                                                     }
        # add target survival node
        phenotype = ('Survival_Time', 'EFO:0000714')
        reasoner_std['query_graph']['nodes']['n1'] = { 'category': 'biolink:PhenotypicFeature',
                                                       'id': '{}'.format(phenotype[1]),
                                                     }
        # link disease to target
        reasoner_std['query_graph']['edges']['e0'] = { 'predicate':'biolink:DiseaseToPhenotypicFeatureAssociation',
                                                       'subject': 'n0',
                                                       'object': 'n1',
                                                       'properties': { 'qualifier':'>=',
                                                                       'days':970
                                                                     }
                                                     }

        handler = ReasonerStdHandler(source_ara='default',
                                     dict_query=reasoner_std)
        queries = handler.buildChpQueries()
        queries = handler.runChpQueries()
        reasoner_std_final = handler.constructDecoratedKG()
        KG = reasoner_std_final['message']['knowledge_graph']
        for edge_key in KG['edges'].keys():
            edge = KG['edges'][edge_key]
            if edge['predicate'] == 'biolink:DiseaseToPhenotypicFeatureAssociation':
                p_survival = edge['has_confidence_level']
        print("probability of survival:",p_survival)

    def test_one_gene(self):
        """ Test with one gene
        """
        # empty response
        reasoner_std = { "query_graph": dict(),
                         "knowledge_graph": dict(),
                         "results": list()
                       }
        # empty query graph
        reasoner_std["query_graph"] = { "edges": dict(),
                                        "nodes": dict()
                                      }
        # empty knowledge graph
        reasoner_std["knowledge_graph"] = { "edges": dict(),
                                            "nodes": dict()
                                          }

        # add in evidence gene
        gene1 = ('RAF1', 'ENSEMBL:ENSG00000132155')
        reasoner_std['query_graph']['nodes']['n0'] = { 'category':'biolink:Gene',
                                                       'id':'{}'.format(gene1[1])
                                                     }
        # add in disease node
        disease = ('Breast_Cancer', 'MONDO:0007254')
        reasoner_std['query_graph']['nodes']['n1'] = { 'category':'biolink:Disease',
                                                       'id':'{}'.format(disease[1])
                                                     }
        # add target survival node
        phenotype = ('Survival_Time', 'EFO:0000714')
        reasoner_std['query_graph']['nodes']['n2'] = { 'category': 'biolink:PhenotypicFeature',
                                                       'id': '{}'.format(phenotype[1]),
                                                     }
        # link genes/drugs to disease
        reasoner_std['query_graph']['edges']['e0'] = { 'predicate':'biolink:GeneToDiseaseAssociation',
                                                       'subject': 'n0',
                                                       'object': 'n1'
                                                     }
        # link disease to target
        reasoner_std['query_graph']['edges']['e1'] = { 'predicate':'biolink:DiseaseToPhenotypicFeatureAssociation',
                                                       'subject': 'n1',
                                                       'object': 'n2',
                                                       'properties': { 'qualifier':'>=',
                                                                       'days':970
                                                                     }
                                                     }

        handler = ReasonerStdHandler(source_ara='default',
                                     dict_query=reasoner_std)
        queries = handler.buildChpQueries()
        queries = handler.runChpQueries()
        reasoner_std_final = handler.constructDecoratedKG()
        KG = reasoner_std_final['message']['knowledge_graph']
        for edge_key in KG['edges'].keys():
            edge = KG['edges'][edge_key]
            if edge['predicate'] == 'biolink:DiseaseToPhenotypicFeatureAssociation':
                p_survival = edge['has_confidence_level']
        print("probability of survival:",p_survival)

    def test_one_drug(self):
        """ Test with one drug
        """
        # empty response
        reasoner_std = { "query_graph": dict(),
                         "knowledge_graph": dict(),
                         "results": list()
                       }
        # empty query graph
        reasoner_std["query_graph"] = { "edges": dict(),
                                        "nodes": dict()
                                      }
        # empty knowledge graph
        reasoner_std["knowledge_graph"] = { "edges": dict(),
                                            "nodes": dict()
                                          }

        # add in evidence drug
        drug = ('CYCLOPHOSPHAMIDE', 'CHEMBL:CHEMBL88')
        reasoner_std['query_graph']['nodes']['n0'] = { 'category':'biolink:Drug',
                                                       'id':'{}'.format(drug[1])
                                                     }
         # add in disease node
        disease = ('Breast_Cancer', 'MONDO:0007254')
        reasoner_std['query_graph']['nodes']['n1'] = { 'category':'biolink:Disease',
                                                       'id':'{}'.format(disease[1])
                                                     }
        # add target survival node
        phenotype = ('Survival_Time', 'EFO:0000714')
        reasoner_std['query_graph']['nodes']['n2'] = { 'category': 'biolink:PhenotypicFeature',
                                                       'id': '{}'.format(phenotype[1]),
                                                     }
        # link genes/drugs to disease
        reasoner_std['query_graph']['edges']['e0'] = { 'predicate':'biolink:ChemicalToDiseaseOrPhenotypicFeatureAssociation',
                                                       'subject': 'n0',
                                                       'object': 'n1'
                                                     }

        # link disease to target
        reasoner_std['query_graph']['edges']['e1'] = { 'predicate':'biolink:DiseaseToPhenotypicFeatureAssociation',
                                                       'subject': 'n1',
                                                       'object': 'n2',
                                                       'properties': { 'qualifier':'>=',
                                                                       'days':970
                                                                     }
                                                     }
        handler = ReasonerStdHandler(source_ara='default',
                                     dict_query=reasoner_std)
        queries = handler.buildChpQueries()
        queries = handler.runChpQueries()
        reasoner_std_final = handler.constructDecoratedKG()
        KG = reasoner_std_final['message']['knowledge_graph']
        for edge_key in KG['edges'].keys():
            edge = KG['edges'][edge_key]
            if edge['predicate'] == 'biolink:DiseaseToPhenotypicFeatureAssociation':
                p_survival = edge['has_confidence_level']
        print("probability of survival:",p_survival)

    def test_no_target(self):
        """ Test with no target (should crash)
        """
        with self.assertRaises(SystemExit) as se:

            # empty response
            reasoner_std = { "query_graph": dict(),
                             "knowledge_graph": dict(),
                             "results": list()
                           }
            # empty query graph
            reasoner_std["query_graph"] = { "edges": dict(),
                                            "nodes": dict()
                                          }
            # empty knowledge graph
            reasoner_std["knowledge_graph"] = { "edges": dict(),
                                                "nodes": dict()
                                              }

            handler = ReasonerStdHandler(source_ara='default',
                                         dict_query=reasoner_std)
            queries = handler.buildChpQueries()
        self.assertEqual(se.exception.code, 'Survival Node not found. Node category must be \'biolink:PhenotypicFeature\' and id must be in: EFO:0000714')

    def test_no_disease(self):
        """ Test with no disease
        """
        with self.assertRaises(SystemExit) as se:

            # empty response
            reasoner_std = { "query_graph": dict(),
                             "knowledge_graph": dict(),
                             "results": list()
                           }
            # empty query graph
            reasoner_std["query_graph"] = { "edges": dict(),
                                            "nodes": dict()
                                          }
            # empty knowledge graph
            reasoner_std["knowledge_graph"] = { "edges": dict(),
                                                "nodes": dict()
                                              }
            # add target survival node
            phenotype = ('Survival_Time', 'EFO:0000714')
            reasoner_std['query_graph']['nodes']['n2'] = { 'category': 'biolink:PhenotypicFeature',
                                                           'id': '{}'.format(phenotype[1]),
                                                         }
            handler = ReasonerStdHandler(source_ara='default',
                                         dict_query=reasoner_std)
            queries = handler.buildChpQueries()
        self.assertEqual(se.exception.code, 'Disease node not found. Node type must be \'biolink:Disease\' and curie must be in: MONDO:0007254')

    def test_default_survival(self):
        """ Test default survival
        """
        # empty response
        reasoner_std = { "query_graph": dict(),
                         "knowledge_graph": dict(),
                         "results": list()
                       }
        # empty query graph
        reasoner_std["query_graph"] = { "edges": dict(),
                                        "nodes": dict()
                                      }
        # empty knowledge graph
        reasoner_std["knowledge_graph"] = { "edges": dict(),
                                            "nodes": dict()
                                          }

        # add in evidence gene
        gene1 = ('RAF1', 'ENSEMBL:ENSG00000132155')
        reasoner_std['query_graph']['nodes']['n0'] = { 'category':'biolink:Gene',
                                                       'id':'{}'.format(gene1[1])
                                                     }
        # add in disease node
        disease = ('Breast_Cancer', 'MONDO:0007254')
        reasoner_std['query_graph']['nodes']['n1'] = { 'category':'biolink:Disease',
                                                       'id':'{}'.format(disease[1])
                                                     }
        # add target survival node
        phenotype = ('Survival_Time', 'EFO:0000714')
        reasoner_std['query_graph']['nodes']['n2'] = { 'category': 'biolink:PhenotypicFeature',
                                                       'id': '{}'.format(phenotype[1]),
                                                     }
        # link genes/drugs to disease
        reasoner_std['query_graph']['edges']['e0'] = { 'predicate':'biolink:GeneToDiseaseAssociation',
                                                       'subject': 'n0',
                                                       'object': 'n1'
                                                     }
        # link disease to target
        reasoner_std['query_graph']['edges']['e1'] = { 'predicate':'biolink:DiseaseToPhenotypicFeatureAssociation',
                                                       'subject': 'n1',
                                                       'object': 'n2',
                                                     }
        handler = ReasonerStdHandler(source_ara='default',
                                     dict_query=reasoner_std)
        queries = handler.buildChpQueries()
        queries = handler.runChpQueries()
        reasoner_std_final = handler.constructDecoratedKG()
        KG = reasoner_std_final['message']['knowledge_graph']
        for edge_key in KG['edges'].keys():
            edge = KG['edges'][edge_key]
            if edge['predicate'] == 'biolink:DiseaseToPhenotypicFeatureAssociation':
                p_survival = edge['has_confidence_level']
        print("probability of survival:",p_survival)

    def test_contribution_property(self):
        """ check if the contribution property returns contributions
        """

        # empty response
        reasoner_std = { "query_graph": dict(),
                         "knowledge_graph": dict(),
                         "results": list()
                       }
        # empty query graph
        reasoner_std["query_graph"] = { "edges": dict(),
                                        "nodes": dict()
                                      }
        # empty knowledge graph
        reasoner_std["knowledge_graph"] = { "edges": dict(),
                                            "nodes": dict()
                                          }

        # add in evidence gene
        gene1 = ('RAF1', 'ENSEMBL:ENSG00000132155')
        reasoner_std['query_graph']['nodes']['n0'] = { 'category':'biolink:Gene',
                                                       'id':'{}'.format(gene1[1])
                                                     }
        # add in disease node
        disease = ('Breast_Cancer', 'MONDO:0007254')
        reasoner_std['query_graph']['nodes']['n1'] = { 'category':'biolink:Disease',
                                                       'id':'{}'.format(disease[1])
                                                     }
        # add target survival node
        phenotype = ('Survival_Time', 'EFO:0000714')
        reasoner_std['query_graph']['nodes']['n2'] = { 'category': 'biolink:PhenotypicFeature',
                                                       'id': '{}'.format(phenotype[1]),
                                                     }
        # link genes/drugs to disease
        reasoner_std['query_graph']['edges']['e0'] = { 'predicate':'biolink:GeneToDiseaseAssociation',
                                                       'subject': 'n0',
                                                       'object': 'n1'
                                                     }
        # link disease to target
        reasoner_std['query_graph']['edges']['e1'] = { 'predicate':'biolink:DiseaseToPhenotypicFeatureAssociation',
                                                       'subject': 'n1',
                                                       'object': 'n2',
                                                       'properties': { 'qualifier':'>=',
                                                                       'days':970,
                                                                       'contributions':True
                                                                     }
                                                     }
        handler = ReasonerStdHandler(source_ara='default',
                                     dict_query=reasoner_std)
        queries = handler.buildChpQueries()
        queries = handler.runChpQueries()
        reasoner_std_final = handler.constructDecoratedKG()
        KG = reasoner_std_final['message']['knowledge_graph']
        for edge_key in KG['edges'].keys():
            edge = KG['edges'][edge_key]
            if edge['predicate'] == 'biolink:DiseaseToPhenotypicFeatureAssociation':
                p_survival = edge['has_confidence_level']
                report = edge['properties']['contributions']
        print("probability of survival:",p_survival)

if __name__ == '__main__':
    unittest.main()
