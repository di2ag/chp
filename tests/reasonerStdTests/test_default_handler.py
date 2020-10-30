import requests
import json
import unittest

from chp.reasoner_std import ReasonerStdHandler

# Tests chp/integrator/exploring_agent handler. Uses 8 test cases:
# 1. Testing negative result normalization
# 2. Normal request with two genes and a drug
# 3. Test with no evidence
# 4. Test with no evidence, but omitting KG and Results (should be handled by handler)
# 5. Test with one gene
# 6. Test with one drug
# 7. Test with no target (should crash)
# 8. Test with no disease
# 9. Test default survival
# 10. check if the contribution property returns contributions

class testExploringAgent(unittest.TestCase):

    # 1. Testing negative result normalization
    def test_negative(self):

        # empty response
        reasoner_std = { "query_graph": dict(),
                         "knowledge_graph": dict(),
                         "results": dict()
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
        reasoner_std['query_graph']['nodes'].append({'type':'gene',
                                                      'id': 'n0',
                                                      'curie':'{}'.format(gene1[1])
                                                     })
        # add in evidence drug
        drug = ('Y', 'CHEMBL:CHEMBL554')
        reasoner_std['query_graph']['nodes'].append({'type':'chemical_substance',
                                                     'id': 'n1',
                                                     'curie':'{}'.format(drug[1])
                                                     })
        # add in disease node
        disease = ('Breast_Cancer', 'MONDO:0007254')
        reasoner_std['query_graph']['nodes'].append({'type':'disease',
                                                     'id': 'n2',
                                                     'curie':'{}'.format(disease[1])
                                                     })
        # add target survival node
        phenotype = ('Survival_Time', 'EFO:0000714')
        reasoner_std['query_graph']['nodes'].append({'type': 'phenotypic_feature',
                                                     'id': 'n3',
                                                     'curie': '{}'.format(phenotype[1]),
                                                     })
        # link genes/drugs to disease
        reasoner_std['query_graph']['edges'].append({ 'type':'gene_to_disease_association',
                                                       'source_id': 'n0',
                                                       'target_id': 'n2',
                                                       'id': 'e0'
                                                     })
        reasoner_std['query_graph']['edges'].append({ 'type':'chemical_to_disease_or_phenotypic_feature_association',
                                                       'source_id': 'n1',
                                                       'target_id': 'n2',
                                                       'id': 'e1'
                                                     })
        # link disease to target
        reasoner_std['query_graph']['edges'].append({ 'type':'disease_to_phenotypic_feature_association',
                                                       'source_id': 'n2',
                                                       'target_id': 'n3',
                                                       'id': 'e3',
                                                       'properties': { 'qualifier':'>=',
                                                                       'value':970
                                                                     }
                                                     })
        handler = ReasonerStdHandler(source_ara='default',
                                     dict_query=reasoner_std)
        queries = handler.buildChpQueries()
        queries = handler.runChpQueries()
        reasoner_std_final = handler.constructDecoratedKG()
        KG = reasoner_std_final['message']['knowledge_graph']
        for edge in KG['edges']:
            if edge['type'] == 'disease_to_phenotypic_feature_association':
                p_survival = edge['has_confidence_level']
        print("probability of survival:",p_survival)

    # 2. Normal request with two genes and a drug
    def test_normal_two_genes_and_drug(self):
        # empty response
        reasoner_std = { "query_graph": dict(),
                         "knowledge_graph": dict(),
                         "results": dict()
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
        gene1 = ('RAF1', 'ENSEMBL:ENSG00000132155')
        reasoner_std['query_graph']['nodes'].append({ 'type':'gene',
                                                      'id': 'n0',
                                                       'curie':'{}'.format(gene1[1])
                                                     })

        gene2 = ('BRCA1', 'ENSEMBL:ENSG00000012048')
        reasoner_std['query_graph']['nodes'].append({ 'type':'gene',
                                                      'id': 'n1',
                                                       'curie':'{}'.format(gene2[1])
                                                     })
        # add in evidence drug
        drug = ('CYCLOPHOSPHAMIDE', 'CHEMBL:CHEMBL88')
        reasoner_std['query_graph']['nodes'].append({ 'type':'chemical_substance',
                                                      'id': 'n2',
                                                       'curie':'{}'.format(drug[1])
                                                     })
        # add in disease node
        disease = ('Breast_Cancer', 'MONDO:0007254')
        reasoner_std['query_graph']['nodes'].append({ 'type':'disease',
                                                      'id': 'n3',
                                                       'curie':'{}'.format(disease[1])
                                                     })
        # add target survival node
        phenotype = ('Survival_Time', 'EFO:0000714')
        reasoner_std['query_graph']['nodes'].append({ 'type': 'phenotypic_feature',
                                                      'id': 'n4',
                                                       'curie': '{}'.format(phenotype[1]),
                                                     })
        # link genes/drugs to disease
        reasoner_std['query_graph']['edges'].append({ 'type':'gene_to_disease_association',
                                                      'id': 'e0',
                                                       'source_id': 'n0',
                                                       'target_id': 'n3'
                                                     })
        reasoner_std['query_graph']['edges'].append({ 'type':'gene_to_disease_association',
                                                      'id': 'e1',
                                                       'source_id': 'n1',
                                                       'target_id': 'n3'
                                                     })
        reasoner_std['query_graph']['edges'].append({ 'type':'chemical_to_disease_or_phenotypic_feature_association',
                                                      'id': 'e2',
                                                       'source_id': 'n2',
                                                       'target_id': 'n3'
                                                     })
        # link disease to target
        reasoner_std['query_graph']['edges'].append({ 'type':'disease_to_phenotypic_feature_association',
                                                      'id': 'e3',
                                                       'source_id': 'n3',
                                                       'target_id': 'n4',
                                                       'properties': { 'qualifier':'>=',
                                                                       'value':970
                                                                     }
                                                     })


        handler = ReasonerStdHandler(source_ara='default',
                                     dict_query=reasoner_std)
        queries = handler.buildChpQueries()
        queries = handler.runChpQueries()
        reasoner_std_final = handler.constructDecoratedKG()
        KG = reasoner_std_final['message']['knowledge_graph']
        for edge in KG['edges']:
            if edge['type'] == 'disease_to_phenotypic_feature_association':
                p_survival = edge['has_confidence_level']
        print("probability of survival:",p_survival)

    # 3. Test with no evidence
    def test_no_evidence(self):
        # empty response
        reasoner_std = { "query_graph": dict(),
                         "knowledge_graph": dict(),
                         "results": dict()
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
        reasoner_std['query_graph']['nodes'].append({ 'type':'disease',
                                                      'id': 'n0',
                                                       'curie':'{}'.format(disease[1])
                                                     })
        # add target survival node
        phenotype = ('Survival_Time', 'EFO:0000714')
        reasoner_std['query_graph']['nodes'].append({ 'type': 'phenotypic_feature',
                                                      'id': 'n1',
                                                       'curie': '{}'.format(phenotype[1]),
                                                     })
        # link disease to target
        reasoner_std['query_graph']['edges'].append({ 'type':'disease_to_phenotypic_feature_association',
                                                      'id': '00',
                                                       'source_id': 'n0',
                                                       'target_id': 'n1',
                                                       'properties': { 'qualifier':'>=',
                                                                       'value':970
                                                                     }
                                                     })

        handler = ReasonerStdHandler(source_ara='default',
                                     dict_query=reasoner_std)
        queries = handler.buildChpQueries()
        queries = handler.runChpQueries()
        reasoner_std_final = handler.constructDecoratedKG()
        KG = reasoner_std_final['message']['knowledge_graph']
        for edge in KG['edges']:
            if edge['type'] == 'disease_to_phenotypic_feature_association':
                p_survival = edge['has_confidence_level']
        print("probability of survival:",p_survival)

    # 4. Test with no evidence, but omitting KG and Results (should be handled by handler)
    def test_no_evidence_omitting_KG_and_results(self):
        # empty response
        reasoner_std = { "query_graph": dict(),
                         "knowledge_graph": dict(),
                         "results": dict()
                       }
        # empty query graph
        reasoner_std["query_graph"] = { "edges": [],
                                        "nodes": []
                                      }

        # add in disease node
        disease = ('Breast_Cancer', 'MONDO:0007254')
        reasoner_std['query_graph']['nodes'].append({ 'type':'disease',
                                                      'id': 'n0',
                                                       'curie':'{}'.format(disease[1])
                                                     })
        # add target survival node
        phenotype = ('Survival_Time', 'EFO:0000714')
        reasoner_std['query_graph']['nodes'].append({ 'type': 'phenotypic_feature',
                                                      'id': 'n1',
                                                       'curie': '{}'.format(phenotype[1]),
                                                     })
        # link disease to target
        reasoner_std['query_graph']['edges'].append({ 'type':'disease_to_phenotypic_feature_association',
                                                      'id': '00',
                                                       'source_id': 'n0',
                                                       'target_id': 'n1',
                                                       'properties': { 'qualifier':'>=',
                                                                       'value':970
                                                                     }
                                                     })

        handler = ReasonerStdHandler(source_ara='default',
                                     dict_query=reasoner_std)
        queries = handler.buildChpQueries()
        queries = handler.runChpQueries()
        reasoner_std_final = handler.constructDecoratedKG()
        KG = reasoner_std_final['message']['knowledge_graph']
        for edge in KG['edges']:
            if edge['type'] == 'disease_to_phenotypic_feature_association':
                p_survival = edge['has_confidence_level']
        print("probability of survival:",p_survival)

    # 5. Test with one gene
    def test_one_gene(self):
        # empty response
        reasoner_std = { "query_graph": dict(),
                         "knowledge_graph": dict(),
                         "results": dict()
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
        gene1 = ('RAF1', 'ENSEMBL:ENSG00000132155')
        reasoner_std['query_graph']['nodes'].append({ 'type':'gene',
                                                      'id': 'n0',
                                                       'curie':'{}'.format(gene1[1])
                                                     })
        # add in disease node
        disease = ('Breast_Cancer', 'MONDO:0007254')
        reasoner_std['query_graph']['nodes'].append({ 'type':'disease',
                                                      'id': 'n1',
                                                       'curie':'{}'.format(disease[1])
                                                     })
        # add target survival node
        phenotype = ('Survival_Time', 'EFO:0000714')
        reasoner_std['query_graph']['nodes'].append({ 'type': 'phenotypic_feature',
                                                      'id': 'n2',
                                                       'curie': '{}'.format(phenotype[1]),
                                                     })
        # link genes/drugs to disease
        reasoner_std['query_graph']['edges'].append({ 'type':'gene_to_disease_association',
                                                      'id': 'e0',
                                                       'source_id': 'n0',
                                                       'target_id': 'n1'
                                                     })
        # link disease to target
        reasoner_std['query_graph']['edges'].append({ 'type':'disease_to_phenotypic_feature_association',
                                                      'id': 'e1',
                                                       'source_id': 'n1',
                                                       'target_id': 'n2',
                                                       'properties': { 'qualifier':'>=',
                                                                       'value':970
                                                                     }
                                                     })

        handler = ReasonerStdHandler(source_ara='default',
                                     dict_query=reasoner_std)
        queries = handler.buildChpQueries()
        queries = handler.runChpQueries()
        reasoner_std_final = handler.constructDecoratedKG()
        KG = reasoner_std_final['message']['knowledge_graph']
        for edge in KG['edges']:
            if edge['type'] == 'disease_to_phenotypic_feature_association':
                p_survival = edge['has_confidence_level']
        print("probability of survival:",p_survival)

    # 6. Test with one drug
    def test_one_drug(self):
        # empty response
        reasoner_std = { "query_graph": dict(),
                         "knowledge_graph": dict(),
                         "results": dict()
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
        reasoner_std['query_graph']['nodes'].append({ 'type':'chemical_substance',
                                                      'id': 'n0',
                                                       'curie':'{}'.format(drug[1])
                                                     })
         # add in disease node
        disease = ('Breast_Cancer', 'MONDO:0007254')
        reasoner_std['query_graph']['nodes'].append({ 'type':'disease',
                                                      'id': 'n1',
                                                       'curie':'{}'.format(disease[1])
                                                     })
        # add target survival node
        phenotype = ('Survival_Time', 'EFO:0000714')
        reasoner_std['query_graph']['nodes'].append({ 'type': 'phenotypic_feature',
                                                      'id': 'n2',
                                                       'curie': '{}'.format(phenotype[1]),
                                                     })
        # link genes/drugs to disease
        reasoner_std['query_graph']['edges'].append({ 'type':'chemical_to_disease_or_phenotypic_feature_association',
                                                      'id': 'e0',
                                                       'source_id': 'n0',
                                                       'target_id': 'n1'
                                                     })

        # link disease to target
        reasoner_std['query_graph']['edges'].append({ 'type':'disease_to_phenotypic_feature_association',
                                                      'id': 'e1',
                                                       'source_id': 'n1',
                                                       'target_id': 'n2',
                                                       'properties': { 'qualifier':'>=',
                                                                       'value':970
                                                                     }
                                                     })
        handler = ReasonerStdHandler(source_ara='default',
                                     dict_query=reasoner_std)
        queries = handler.buildChpQueries()
        queries = handler.runChpQueries()
        reasoner_std_final = handler.constructDecoratedKG()
        KG = reasoner_std_final['message']['knowledge_graph']
        for edge in KG['edges']:
            if edge['type'] == 'disease_to_phenotypic_feature_association':
                p_survival = edge['has_confidence_level']
        print("probability of survival:",p_survival)

    # 7. Test with no target (should crash)
    def test_no_target(self):
        with self.assertRaises(SystemExit) as se:
            # empty response
            reasoner_std = { "query_graph": dict(),
                             "knowledge_graph": dict(),
                             "results": dict()
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

            handler = ReasonerStdHandler(source_ara='default',
                                         dict_query=reasoner_std)
            queries = handler.buildChpQueries()
        self.assertEqual(se.exception.code, 'Survival Node not found. Node type muse be \'PhenotypicFeature\' and curie must be in: EFO:0000714')

    # 8. Test with no disease
    def test_no_disease(self):
        with self.assertRaises(SystemExit) as se:
            # empty response
            reasoner_std = { "query_graph": dict(),
                             "knowledge_graph": dict(),
                             "results": dict()
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
            reasoner_std['query_graph']['nodes'].append({ 'type': 'phenotypic_feature',
                                                          'id': 'n2',
                                                           'curie': '{}'.format(phenotype[1]),
                                                         })
            handler = ReasonerStdHandler(source_ara='default',
                                         dict_query=reasoner_std)
            queries = handler.buildChpQueries()
        self.assertEqual(se.exception.code, 'Disease node not found. Node type must be \'disease\' and curie must be in: MONDO:0007254')

    # 9. Test default survival
    def test_default_survival(self):
        # empty response
        reasoner_std = { "query_graph": dict(),
                         "knowledge_graph": dict(),
                         "results": dict()
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
        gene1 = ('RAF1', 'ENSEMBL:ENSG00000132155')
        reasoner_std['query_graph']['nodes'].append({ 'type':'gene',
                                                      'id': 'n0',
                                                       'curie':'{}'.format(gene1[1])
                                                     })
        # add in disease node
        disease = ('Breast_Cancer', 'MONDO:0007254')
        reasoner_std['query_graph']['nodes'].append({ 'type':'disease',
                                                      'id': 'n1',
                                                       'curie':'{}'.format(disease[1])
                                                     })
        # add target survival node
        phenotype = ('Survival_Time', 'EFO:0000714')
        reasoner_std['query_graph']['nodes'].append({ 'type': 'phenotypic_feature',
                                                      'id': 'n2',
                                                       'curie': '{}'.format(phenotype[1]),
                                                     })
        # link genes/drugs to disease
        reasoner_std['query_graph']['edges'].append({ 'type':'gene_to_disease_association',
                                                      'id': 'e0',
                                                       'source_id': 'n0',
                                                       'target_id': 'n1'
                                                     })
        # link disease to target
        reasoner_std['query_graph']['edges'].append({ 'type':'disease_to_phenotypic_feature_association',
                                                      'id': 'e1',
                                                       'source_id': 'n1',
                                                       'target_id': 'n2',
                                                     })

        handler = ReasonerStdHandler(source_ara='default',
                                     dict_query=reasoner_std)
        queries = handler.buildChpQueries()
        queries = handler.runChpQueries()
        reasoner_std_final = handler.constructDecoratedKG()
        KG = reasoner_std_final['message']['knowledge_graph']
        for edge in KG['edges']:
            if edge['type'] == 'disease_to_phenotypic_feature_association':
                p_survival = edge['has_confidence_level']
        print("probability of survival:",p_survival)

    # 10. check if the contribution property returns contributions
    def test_contribution_property(self):
        # empty response
        reasoner_std = { "query_graph": dict(),
                         "knowledge_graph": dict(),
                         "results": dict()
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
        gene1 = ('RAF1', 'ENSEMBL:ENSG00000132155')
        reasoner_std['query_graph']['nodes'].append({ 'type':'gene',
                                                      'id': 'n0',
                                                       'curie':'{}'.format(gene1[1])
                                                     })
        # add in disease node
        disease = ('Breast_Cancer', 'MONDO:0007254')
        reasoner_std['query_graph']['nodes'].append({ 'type':'disease',
                                                      'id': 'n1',
                                                       'curie':'{}'.format(disease[1])
                                                     })
        # add target survival node
        phenotype = ('Survival_Time', 'EFO:0000714')
        reasoner_std['query_graph']['nodes'].append({ 'type': 'phenotypic_feature',
                                                      'id': 'n2',
                                                       'curie': '{}'.format(phenotype[1]),
                                                     })
        # link genes/drugs to disease
        reasoner_std['query_graph']['edges'].append({ 'type':'gene_to_disease_association',
                                                      'id': 'e0',
                                                       'source_id': 'n0',
                                                       'target_id': 'n1'
                                                     })
        # link disease to target
        reasoner_std['query_graph']['edges'].append({ 'type':'disease_to_phenotypic_feature_association',
                                                      'id': 'e1',
                                                       'source_id': 'n1',
                                                       'target_id': 'n2',
                                                       'properties': { 'qualifier':'>=',
                                                                       'value':970,
                                                                       'contributions':True
                                                                     }
                                                     })
        handler = ReasonerStdHandler(source_ara='default',
                                     dict_query=reasoner_std)
        queries = handler.buildChpQueries()
        queries = handler.runChpQueries()
        reasoner_std_final = handler.constructDecoratedKG()
        KG = reasoner_std_final['message']['knowledge_graph']
        for edge in KG['edges']:
            if edge['type'] == 'disease_to_phenotypic_feature_association':
                p_survival = edge['has_confidence_level']
                report = edge['Description']
        print("probability of survival:",p_survival)

if __name__ == '__main__':
    unittest.main()
