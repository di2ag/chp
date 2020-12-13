import requests
import json
import unittest

from chp.reasoner_std import ReasonerStdHandler

# Tests chp/integrator/ranking_agent handler. Uses 8 test cases:
# 1. Testing negative result normalization
# 2. Normal request with two genes and a drug
# 3. Test with no evidence
# 4. Test with no evidence, but omitting KG and Results (should be handled by handler)
# 5. Test with one gene
# 6. Test with one drug
# 7. Test with no target (should crash)
# 8. Test with no disease

class testUnsecretAgent(unittest.TestCase):

    # 1. Testing negative result normalization
    def test_negative(self):
        # empty response
        reasoner_std = { "query_graph": dict(),
                         "knowledge_graph": dict(),
                         "response": dict()
                       }
        # empty query graph
        reasoner_std["query_graph"] = { "edges": [],
                                        "nodes": []
                                      }
        # empty knowledge graph
        reasoner_std["knowledge_graph"] = { "edges": [],
                                            "nodes": []
                                          }
        # empty response graph
        reasoner_std["results"] = { "node_bindings": [],
                                    "edge_bindings": []
                                  }
        # add in evidence gene
        gene1 = ('X', 'ENSEMBL:ENSG00000007174')
        reasoner_std['query_graph']['nodes'].append({ 'id':'n{}'.format('0'),
                                                      'type':'Gene',
                                                      'curie':'{}'.format(gene1[1])
                                                   })
        # add in evidence drug
        drug = ('Y', 'CHEMBL:CHEMBL554')
        reasoner_std['query_graph']['nodes'].append({ 'id':'n{}'.format('1'),
                                                      'type':'Drug',
                                                      'curie':'{}'.format(drug[1])
                                                   })
        # add in disease node
        disease = ('Breast_Cancer', 'MONDO:0007254')
        reasoner_std['query_graph']['nodes'].append({ 'id':'n{}'.format('2'),
                                                      'type':'disease',
                                                      'curie':'{}'.format(disease[1])
                                                   })
        # add target survival node
        phenotype = ('Survival_Time', 'EFO:0000714')
        reasoner_std['query_graph']['nodes'].append({ 'id': 'n{}'.format('3'),
                                                      'type': 'PhenotypicFeature',
                                                      'curie': '{}'.format(phenotype[1]),
                                                   })
        # link genes/drugs to disease
        reasoner_std['query_graph']['edges'].append({ 'id':'e{}'.format('0'),
                                                      'type':'gene_to_disease_association',
                                                      'source_id':'n{}'.format('0'),
                                                      'target_id':'n{}'.format('2')
                                                    })
        reasoner_std['query_graph']['edges'].append({ 'id':'e{}'.format('1'),
                                                      'type':'chemical_to_disease_or_phenotypic_feature_association',
                                                      'source_id':'n{}'.format('1'),
                                                      'target_id':'n{}'.format('2')
                                                   })
        # link disease to target
        reasoner_std['query_graph']['edges'].append({ 'id':'e{}'.format('2'),
                                                      'type':'disease_to_phenotype_association',
                                                      'value':970,
                                                      'source_id':'n{}'.format('2'),
                                                      'target_id':'n{}'.format('3')
                                                   })
        handler = ReasonerStdHandler(source_ara='unsecret',
                                     dict_query=reasoner_std)
        queries = handler.buildChpQueries()
        queries = handler.runChpQueries()
        reasoner_std_final = handler.constructDecoratedKG()
        KG = reasoner_std_final['knowledge_graph']
        for edge in KG['edges']:
            if edge['type'] == 'disease_to_phenotype_association':
                p_survival = edge['has_confidence_level']
        print("probability of survival:",p_survival)

    # 2. Normal request with two genes and a drug
    def test_normal_two_genes_and_drug(self):
        # empty response
        reasoner_std = { "query_graph": dict(),
                         "knowledge_graph": dict(),
                         "response": dict()
                       }
        # empty query graph
        reasoner_std["query_graph"] = { "edges": [],
                                        "nodes": []
                                      }
        # empty knowledge graph
        reasoner_std["knowledge_graph"] = { "edges": [],
                                            "nodes": []
                                          }
        # empty response graph
        reasoner_std["results"] = { "node_bindings": [],
                                    "edge_bindings": []
                                  }
        # add in evidence gene
        gene1 = ('TP53', 'ENSEMBL:ENSG00000141510')
        reasoner_std['query_graph']['nodes'].append({ 'id':'n{}'.format('0'),
                                                      'type':'Gene',
                                                      'curie':'{}'.format(gene1[1])
                                                   })
        # add in second evidence gene
        gene2 = ('CDH1', 'ENSEMBL:ENSG00000039068')
        reasoner_std['query_graph']['nodes'].append({ 'id':'n{}'.format('1'),
                                                      'type':'Gene',
                                                      'curie':'{}'.format(gene2[1])
                                                   })
        # add in evidence drug
        drug = ('CYCLOPHOSPHAMIDE', 'CHEMBL:CHEMBL88')
        reasoner_std['query_graph']['nodes'].append({ 'id':'n{}'.format('2'),
                                                      'type':'Drug',
                                                      'curie':'{}'.format(drug[1])
                                                   })
        # add in disease node
        disease = ('Breast_Cancer', 'MONDO:0007254')
        reasoner_std['query_graph']['nodes'].append({ 'id':'n{}'.format('3'),
                                                      'type':'disease',
                                                      'curie':'{}'.format(disease[1])
                                                   })
        # add target survival node
        phenotype = ('Survival_Time', 'EFO:0000714')
        reasoner_std['query_graph']['nodes'].append({ 'id': 'n{}'.format('4'),
                                                      'type': 'PhenotypicFeature',
                                                      'curie': '{}'.format(phenotype[1]),
                                                   })
        # link genes/drugs to disease
        reasoner_std['query_graph']['edges'].append({ 'id':'e{}'.format('0'),
                                                      'type':'gene_to_disease_association',
                                                      'source_id':'n{}'.format('0'),
                                                      'target_id':'n{}'.format('3')
                                                    })
        reasoner_std['query_graph']['edges'].append({ 'id':'e{}'.format('1'),
                                                      'type':'gene_to_disease_association',
                                                      'source_id':'n{}'.format('1'),
                                                      'target_id':'n{}'.format('3')
                                                   })
        reasoner_std['query_graph']['edges'].append({ 'id':'e{}'.format('2'),
                                                      'type':'chemical_to_disease_or_phenotypic_feature_association',
                                                      'source_id':'n{}'.format('2'),
                                                      'target_id':'n{}'.format('3')
                                                   })
        # link disease to target
        reasoner_std['query_graph']['edges'].append({ 'id':'e{}'.format('3'),
                                                      'type':'disease_to_phenotype_association',
                                                      'value':1300,
                                                      'source_id':'n{}'.format('3'),
                                                      'target_id':'n{}'.format('4')
                                                   })
        handler = ReasonerStdHandler(source_ara='unsecret',
                                     dict_query=reasoner_std)
        queries = handler.buildChpQueries()
        queries = handler.runChpQueries()
        reasoner_std_final = handler.constructDecoratedKG()
        KG = reasoner_std_final['knowledge_graph']
        for edge in KG['edges']:
            if edge['type'] == 'disease_to_phenotype_association':
                p_survival = edge['has_confidence_level']
        print("probability of survival:",p_survival)

    # 3. Test with no evidence
    def test_no_evidence(self):
        with self.assertRaises(SystemExit) as se:
            # empty response
            reasoner_std = { "query_graph": dict(),
                             "knowledge_graph": dict(),
                             "response": dict()
                           }
            # empty query graph
            reasoner_std["query_graph"] = { "edges": [],
                                            "nodes": []
                                          }
            # empty knowledge graph
            reasoner_std["knowledge_graph"] = { "edges": [],
                                                "nodes": []
                                              }
            # empty response graph
            reasoner_std["results"] = { "node_bindings": [],
                                        "edge_bindings": []
                                      }
            # add in disease node
            disease = ('Breast_Cancer', 'MONDO:0007254')
            reasoner_std['query_graph']['nodes'].append({ 'id':'n{}'.format('0'),
                                                          'type':'disease',
                                                          'curie':'{}'.format(disease[1])
                                                       })
            # add target survival node
            phenotype = ('Survival_Time', 'EFO:0000714')
            reasoner_std['query_graph']['nodes'].append({ 'id': 'n{}'.format('1'),
                                                          'type': 'PhenotypicFeature',
                                                          'curie': '{}'.format(phenotype[1]),
                                                       })
            # link disease to target
            reasoner_std['query_graph']['edges'].append({ 'id':'e{}'.format('0'),
                                                          'type':'disease_to_phenotype_association',
                                                          'value':970,
                                                          'source_id':'n{}'.format('0'),
                                                          'target_id':'n{}'.format('1')
                                                       })
            handler = ReasonerStdHandler(source_ara='unsecret',
                                         dict_query=reasoner_std)
            queries = handler.buildChpQueries()
            queries = handler.runChpQueries()
            reasoner_std_final = handler.constructDecoratedKG()
            KG = reasoner_std_final['knowledge_graph']
            for edge in KG['edges']:
                if edge['type'] == 'disease_to_phenotype_association':
                    p_survival = edge['has_confidence_level']
            print("probability of survival:",p_survival)
        self.assertEqual(se.exception.code, 'Needs at least 1 gene')

    # 4. Test with no evidence, but omitting KG and Results (should be handled by handler)
    def test_no_evidence_omitting_KG_and_results(self):
        with self.assertRaises(SystemExit) as se:
            # empty response
            reasoner_std = { "query_graph": dict(),
                           }
            # empty query graph
            reasoner_std["query_graph"] = { "edges": [],
                                            "nodes": []
                                          }
            # add in disease node
            disease = ('Breast_Cancer', 'MONDO:0007254')
            reasoner_std['query_graph']['nodes'].append({ 'id':'n{}'.format('0'),
                                                          'type':'disease',
                                                          'curie':'{}'.format(disease[1])
                                                       })
            # add target survival node
            phenotype = ('Survival_Time', 'EFO:0000714')
            reasoner_std['query_graph']['nodes'].append({ 'id': 'n{}'.format('1'),
                                                          'type': 'PhenotypicFeature',
                                                          'curie': '{}'.format(phenotype[1]),
                                                       })
            # link disease to target
            reasoner_std['query_graph']['edges'].append({ 'id':'e{}'.format('0'),
                                                          'type':'disease_to_phenotype_association',
                                                          'value':970,
                                                          'source_id':'n{}'.format('0'),
                                                          'target_id':'n{}'.format('1')
                                                       })
            handler = ReasonerStdHandler(source_ara='unsecret',
                                         dict_query=reasoner_std)
            queries = handler.buildChpQueries()
            queries = handler.runChpQueries()
            reasoner_std_final = handler.constructDecoratedKG()
            KG = reasoner_std_final['knowledge_graph']
            for edge in KG['edges']:
                if edge['type'] == 'disease_to_phenotype_association':
                    p_survival = edge['has_confidence_level']
            print("probability of survival:",p_survival)
        self.assertEqual(se.exception.code, 'Needs at least 1 gene')

    # 5. Test with one gene
    def test_one_gene(self):
        # empty response
        reasoner_std = { "query_graph": dict(),
                         "knowledge_graph": dict(),
                         "response": dict()
                       }
        # empty query graph
        reasoner_std["query_graph"] = { "edges": [],
                                        "nodes": []
                                      }
        # empty knowledge graph
        reasoner_std["knowledge_graph"] = { "edges": [],
                                            "nodes": []
                                          }
        # empty response graph
        reasoner_std["results"] = { "node_bindings": [],
                                    "edge_bindings": []
                                  }
        # add in evidence gene
        gene = ('RAF1', 'ENSEMBL:ENSG00000132155')
        reasoner_std['query_graph']['nodes'].append({ 'id':'n{}'.format('0'),
                                                      'type':'Gene',
                                                      'curie':'{}'.format(gene[1])
                                                   })
        # add in disease node
        disease = ('Breast_Cancer', 'MONDO:0007254')
        reasoner_std['query_graph']['nodes'].append({ 'id':'n{}'.format('1'),
                                                      'type':'disease',
                                                      'curie':'{}'.format(disease[1])
                                                   })
        # add target survival node
        phenotype = ('Survival_Time', 'EFO:0000714')
        reasoner_std['query_graph']['nodes'].append({ 'id': 'n{}'.format('2'),
                                                      'type': 'PhenotypicFeature',
                                                      'curie': '{}'.format(phenotype[1]),
                                                   })
        # link genes
        reasoner_std['query_graph']['edges'].append({ 'id':'e{}'.format('0'),
                                                      'type':'gene_to_disease_association',
                                                      'source_id':'n{}'.format('0'),
                                                      'target_id':'n{}'.format('1')
                                                    })
        # link disease to target
        reasoner_std['query_graph']['edges'].append({ 'id':'e{}'.format('1'),
                                                      'type':'disease_to_phenotype_association',
                                                      'value':970,
                                                      'source_id':'n{}'.format('1'),
                                                      'target_id':'n{}'.format('2')
                                                   })
        handler = ReasonerStdHandler(source_ara='unsecret',
                                     dict_query=reasoner_std)
        queries = handler.buildChpQueries()
        queries = handler.runChpQueries()
        reasoner_std_final = handler.constructDecoratedKG()
        KG = reasoner_std_final['knowledge_graph']
        for edge in KG['edges']:
            if edge['type'] == 'disease_to_phenotype_association':
                p_survival = edge['has_confidence_level']
        print("probability of survival:",p_survival)

    # 6. Test with one drug
    def test_one_drug(self):
        with self.assertRaises(SystemExit) as se:
            # empty response
            reasoner_std = { "query_graph": dict(),
                             "knowledge_graph": dict(),
                             "response": dict()
                           }
            # empty query graph
            reasoner_std["query_graph"] = { "edges": [],
                                            "nodes": []
                                          }
            # empty knowledge graph
            reasoner_std["knowledge_graph"] = { "edges": [],
                                                "nodes": []
                                              }
            # empty response graph
            reasoner_std["results"] = { "node_bindings": [],
                                        "edge_bindings": []
                                      }
            # add in evidence drug
            drug = ('CYCLOPHOSPHAMIDE', 'CHEMBL:CHEMBL88')
            reasoner_std['query_graph']['nodes'].append({ 'id':'n{}'.format('0'),
                                                          'type':'Drug',
                                                          'curie':'{}'.format(drug[1])
                                                       })    
            # add in disease node
            disease = ('Breast_Cancer', 'MONDO:0007254')
            reasoner_std['query_graph']['nodes'].append({ 'id':'n{}'.format('1'),
                                                          'type':'disease',
                                                          'curie':'{}'.format(disease[1])
                                                       })
            # add target survival node
            phenotype = ('Survival_Time', 'EFO:0000714')
            reasoner_std['query_graph']['nodes'].append({ 'id': 'n{}'.format('2'),
                                                          'type': 'PhenotypicFeature',
                                                          'curie': '{}'.format(phenotype[1]),
                                                       })
            # link drug to disease
            reasoner_std['query_graph']['edges'].append({ 'id':'e{}'.format('0'),
                                                          'type':'chemical_to_disease_or_phenotypic_feature_association',
                                                          'source_id':'n{}'.format('0'),
                                                          'target_id':'n{}'.format('1')
                                                       })
            # link disease to target
            reasoner_std['query_graph']['edges'].append({ 'id':'e{}'.format('1'),
                                                          'type':'disease_to_phenotype_association',
                                                          'value':970,
                                                          'source_id':'n{}'.format('1'),
                                                          'target_id':'n{}'.format('2')
                                                       })
            handler = ReasonerStdHandler(source_ara='unsecret',
                                         dict_query=reasoner_std)
            queries = handler.buildChpQueries()
            queries = handler.runChpQueries()
            reasoner_std_final = handler.constructDecoratedKG()
            KG = reasoner_std_final['knowledge_graph']
            for edge in KG['edges']:
                if edge['type'] == 'disease_to_phenotype_association':
                    p_survival = edge['has_confidence_level']
            print("probability of survival:",p_survival)
        self.assertEqual(se.exception.code, 'Needs at least 1 gene')


    # 7. Test with no target (should crash)
    def test_no_target(self):
        with self.assertRaises(SystemExit) as se:
            # empty response
            reasoner_std = { "query_graph": dict(),
                             "knowledge_graph": dict(),
                             "response": dict()
                           }
            # empty query graph
            reasoner_std["query_graph"] = { "edges": [],
                                            "nodes": []
                                          }
            # empty knowledge graph
            reasoner_std["knowledge_graph"] = { "edges": [],
                                                "nodes": []
                                              }
            # empty response graph
            reasoner_std["results"] = { "node_bindings": [],
                                        "edge_bindings": []
                                      }
            handler = ReasonerStdHandler(source_ara='exploring',
                                         dict_query=reasoner_std)
            queries = handler.buildChpQueries()
        self.assertEqual(se.exception.code, 'Survival Node not found. Node type muse be \'PhenotypicFeature\' and curie must be in: EFO:0000714')

    # 8. Test with no disease
    def test_no_disease(self):
        with self.assertRaises(SystemExit) as se:
            # empty response
            reasoner_std = { "query_graph": dict(),
                             "knowledge_graph": dict(),
                             "response": dict()
                           }
            # empty query graph
            reasoner_std["query_graph"] = { "edges": [],
                                            "nodes": []
                                          }
            # empty knowledge graph
            reasoner_std["knowledge_graph"] = { "edges": [],
                                                "nodes": []
                                              }
            # empty response graph
            reasoner_std["results"] = { "node_bindings": [],
                                        "edge_bindings": []
                                      }
            # add target survival node
            phenotype = ('Survival_Time', 'EFO:0000714')
            reasoner_std['query_graph']['nodes'].append({ 'id': 'n{}'.format('1'),
                                                          'type': 'PhenotypicFeature',
                                                          'curie': '{}'.format(phenotype[1]),
                                                       })
            handler = ReasonerStdHandler(source_ara='unsecret',
                                         dict_query=reasoner_std)
            queries = handler.buildChpQueries()
        self.assertEqual(se.exception.code, 'Disease node not found. Node type must be \'disease\' and curie must be in: MONDO:0007254')

if __name__ == '__main__':
    unittest.main()
