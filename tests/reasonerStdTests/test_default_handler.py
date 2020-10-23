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
        # empty response graph
        reasoner_std["results"] = [{ "node_bindings": dict(),
                                     "edge_bindings": dict()
                                  }]

        # add in evidence gene
        gene1 = ('X', 'ENSEMBL:ENSG00000007174')
        reasoner_std['query_graph']['nodes']['n0'] = { 'type':'gene',
                                                       'curie':'{}'.format(gene1[1])
                                                     }
        # add in evidence drug
        drug = ('Y', 'CHEMBL:CHEMBL554')
        reasoner_std['query_graph']['nodes']['n1'] = { 'type':'drug',
                                                       'curie':'{}'.format(drug[1])
                                                     }
        # add in disease node
        disease = ('Breast_Cancer', 'MONDO:0007254')
        reasoner_std['query_graph']['nodes']['n2'] = { 'type':'disease',
                                                       'curie':'{}'.format(disease[1])
                                                     }
        # add target survival node
        phenotype = ('Survival_Time', 'EFO:0000714')
        reasoner_std['query_graph']['nodes']['n3'] = { 'type': 'phenotypicfeature',
                                                       'curie': '{}'.format(phenotype[1]),
                                                     }
        # link genes/drugs to disease
        reasoner_std['query_graph']['edges']['e0'] = { 'type':'gene_to_disease_association',
                                                       'source_id': 'n0',
                                                       'target_id': 'n2'
                                                     }
        reasoner_std['query_graph']['edges']['e1'] = { 'type':'chemical_to_disease_or_phenotypic_feature_association',
                                                       'source_id': 'n1',
                                                       'target_id': 'n2'
                                                     }
        # link disease to target
        reasoner_std['query_graph']['edges']['e2'] = { 'type':'disease_to_phenotype_association',
                                                       'source_id': 'n2',
                                                       'target_id': 'n3',
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
            if edge['type'] == 'disease_to_phenotype_association':
                p_survival = edge['has_confidence_level']
        print("probability of survival:",p_survival)

    # 2. Normal request with two genes and a drug
    def test_normal_two_genes_and_drug(self):

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
        # empty response graph
        reasoner_std["results"] = [{ "node_bindings": dict(),
                                     "edge_bindings": dict()
                                  }]

        # add in evidence gene
        gene1 = ('RAF1', 'ENSEMBL:ENSG00000132155')
        reasoner_std['query_graph']['nodes']['n0'] = { 'type':'gene',
                                                       'curie':'{}'.format(gene1[1])
                                                     }

        gene2 = ('BRCA1', 'ENSEMBL:ENSG00000012048')
        reasoner_std['query_graph']['nodes']['n1'] = { 'type':'gene',
                                                       'curie':'{}'.format(gene2[1])
                                                     }
        # add in evidence drug
        drug = ('CYCLOPHOSPHAMIDE', 'CHEMBL:CHEMBL88')
        reasoner_std['query_graph']['nodes']['n2'] = { 'type':'drug',
                                                       'curie':'{}'.format(drug[1])
                                                     }
        # add in disease node
        disease = ('Breast_Cancer', 'MONDO:0007254')
        reasoner_std['query_graph']['nodes']['n3'] = { 'type':'disease',
                                                       'curie':'{}'.format(disease[1])
                                                     }
        # add target survival node
        phenotype = ('Survival_Time', 'EFO:0000714')
        reasoner_std['query_graph']['nodes']['n4'] = { 'type': 'phenotypicfeature',
                                                       'curie': '{}'.format(phenotype[1]),
                                                     }
        # link genes/drugs to disease
        reasoner_std['query_graph']['edges']['e0'] = { 'type':'gene_to_disease_association',
                                                       'source_id': 'n0',
                                                       'target_id': 'n3'
                                                     }
        reasoner_std['query_graph']['edges']['e1'] = { 'type':'gene_to_disease_association',
                                                       'source_id': 'n1',
                                                       'target_id': 'n3'
                                                     }
        reasoner_std['query_graph']['edges']['e2'] = { 'type':'chemical_to_disease_or_phenotypic_feature_association',
                                                       'source_id': 'n2',
                                                       'target_id': 'n3'
                                                     }
        # link disease to target
        reasoner_std['query_graph']['edges']['e3'] = { 'type':'disease_to_phenotype_association',
                                                       'source_id': 'n3',
                                                       'target_id': 'n4',
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
            if edge['type'] == 'disease_to_phenotype_association':
                p_survival = edge['has_confidence_level']
        print("probability of survival:",p_survival)

    # 3. Test with no evidence
    def test_no_evidence(self):
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
        # empty response graph
        reasoner_std["results"] = [{ "node_bindings": dict(),
                                     "edge_bindings": dict()
                                  }]

        # add in disease node
        disease = ('Breast_Cancer', 'MONDO:0007254')
        reasoner_std['query_graph']['nodes']['n0'] = { 'type':'disease',
                                                       'curie':'{}'.format(disease[1])
                                                     }
        # add target survival node
        phenotype = ('Survival_Time', 'EFO:0000714')
        reasoner_std['query_graph']['nodes']['n1'] = { 'type': 'phenotypicfeature',
                                                       'curie': '{}'.format(phenotype[1]),
                                                     }
        # link disease to target
        reasoner_std['query_graph']['edges']['e0'] = { 'type':'disease_to_phenotype_association',
                                                       'source_id': 'n0',
                                                       'target_id': 'n1',
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
            if edge['type'] == 'disease_to_phenotype_association':
                p_survival = edge['has_confidence_level']
        print("probability of survival:",p_survival)

    # 4. Test with no evidence, but omitting KG and Results (should be handled by handler)
    def test_no_evidence_omitting_KG_and_results(self):
        # empty response
        reasoner_std = { "query_graph": dict(),
                       }

        # empty query graph
        reasoner_std["query_graph"] = { "edges": dict(),
                                        "nodes": dict()
                                      }
        # empty knowledge graph
        reasoner_std["knowledge_graph"] = { "edges": dict(),
                                            "nodes": dict()
                                          }
        # empty response graph
        reasoner_std["results"] = [{ "node_bindings": dict(),
                                     "edge_bindings": dict()
                                  }]

        # add in disease node
        disease = ('Breast_Cancer', 'MONDO:0007254')
        reasoner_std['query_graph']['nodes']['n0'] = { 'type':'disease',
                                                       'curie':'{}'.format(disease[1])
                                                     }
        # add target survival node
        phenotype = ('Survival_Time', 'EFO:0000714')
        reasoner_std['query_graph']['nodes']['n1'] = { 'type': 'phenotypicfeature',
                                                       'curie': '{}'.format(phenotype[1]),
                                                     }
        # link disease to target
        reasoner_std['query_graph']['edges']['e0'] = { 'type':'disease_to_phenotype_association',
                                                       'source_id': 'n0',
                                                       'target_id': 'n1',
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
            if edge['type'] == 'disease_to_phenotype_association':
                p_survival = edge['has_confidence_level']
        print("probability of survival:",p_survival)

    # 5. Test with one gene
    def test_one_gene(self):
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
        # empty response graph
        reasoner_std["results"] = [{ "node_bindings": dict(),
                                     "edge_bindings": dict()
                                  }]

        # add in evidence gene
        gene1 = ('RAF1', 'ENSEMBL:ENSG00000132155')
        reasoner_std['query_graph']['nodes']['n0'] = { 'type':'gene',
                                                       'curie':'{}'.format(gene1[1])
                                                     }
        # add in disease node
        disease = ('Breast_Cancer', 'MONDO:0007254')
        reasoner_std['query_graph']['nodes']['n1'] = { 'type':'disease',
                                                       'curie':'{}'.format(disease[1])
                                                     }
        # add target survival node
        phenotype = ('Survival_Time', 'EFO:0000714')
        reasoner_std['query_graph']['nodes']['n2'] = { 'type': 'phenotypicfeature',
                                                       'curie': '{}'.format(phenotype[1]),
                                                     }
        # link genes/drugs to disease
        reasoner_std['query_graph']['edges']['e0'] = { 'type':'gene_to_disease_association',
                                                       'source_id': 'n0',
                                                       'target_id': 'n1'
                                                     }
        # link disease to target
        reasoner_std['query_graph']['edges']['e1'] = { 'type':'disease_to_phenotype_association',
                                                       'source_id': 'n1',
                                                       'target_id': 'n2',
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
            if edge['type'] == 'disease_to_phenotype_association':
                p_survival = edge['has_confidence_level']
        print("probability of survival:",p_survival)

    # 6. Test with one drug
    def test_one_drug(self):

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
        # empty response graph
        reasoner_std["results"] = [{ "node_bindings": dict(),
                                     "edge_bindings": dict()
                                  }]

        # add in evidence drug
        drug = ('CYCLOPHOSPHAMIDE', 'CHEMBL:CHEMBL88')
        reasoner_std['query_graph']['nodes']['n0'] = { 'type':'drug',
                                                       'curie':'{}'.format(drug[1])
                                                     }
         # add in disease node
        disease = ('Breast_Cancer', 'MONDO:0007254')
        reasoner_std['query_graph']['nodes']['n1'] = { 'type':'disease',
                                                       'curie':'{}'.format(disease[1])
                                                     }
        # add target survival node
        phenotype = ('Survival_Time', 'EFO:0000714')
        reasoner_std['query_graph']['nodes']['n2'] = { 'type': 'phenotypicfeature',
                                                       'curie': '{}'.format(phenotype[1]),
                                                     }
        # link genes/drugs to disease
        reasoner_std['query_graph']['edges']['e0'] = { 'type':'chemical_to_disease_or_phenotypic_feature_association',
                                                       'source_id': 'n0',
                                                       'target_id': 'n1'
                                                     }

        # link disease to target
        reasoner_std['query_graph']['edges']['e1'] = { 'type':'disease_to_phenotype_association',
                                                       'source_id': 'n1',
                                                       'target_id': 'n2',
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
            if edge['type'] == 'disease_to_phenotype_association':
                p_survival = edge['has_confidence_level']
        print("probability of survival:",p_survival)

    # 7. Test with no target (should crash)
    def test_no_target(self):
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
            # empty response graph
            reasoner_std["results"] = [{ "node_bindings": dict(),
                                         "edge_bindings": dict()
                                      }]

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
            # empty response graph
            reasoner_std["results"] = [{ "node_bindings": dict(),
                                         "edge_bindings": dict()
                                      }]
            # add target survival node
            phenotype = ('Survival_Time', 'EFO:0000714')
            reasoner_std['query_graph']['nodes']['n2'] = { 'type': 'phenotypicfeature',
                                                           'curie': '{}'.format(phenotype[1]),
                                                         }
            handler = ReasonerStdHandler(source_ara='default',
                                         dict_query=reasoner_std)
            queries = handler.buildChpQueries()
        self.assertEqual(se.exception.code, 'Disease node not found. Node type must be \'disease\' and curie must be in: MONDO:0007254')

    # 9. Test default survival
    def test_default_survival(self):
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
        # empty response graph
        reasoner_std["results"] = [{ "node_bindings": dict(),
                                     "edge_bindings": dict()
                                  }]

        # add in evidence gene
        gene1 = ('RAF1', 'ENSEMBL:ENSG00000132155')
        reasoner_std['query_graph']['nodes']['n0'] = { 'type':'gene',
                                                       'curie':'{}'.format(gene1[1])
                                                     }
        # add in disease node
        disease = ('Breast_Cancer', 'MONDO:0007254')
        reasoner_std['query_graph']['nodes']['n1'] = { 'type':'disease',
                                                       'curie':'{}'.format(disease[1])
                                                     }
        # add target survival node
        phenotype = ('Survival_Time', 'EFO:0000714')
        reasoner_std['query_graph']['nodes']['n2'] = { 'type': 'phenotypicfeature',
                                                       'curie': '{}'.format(phenotype[1]),
                                                     }
        # link genes/drugs to disease
        reasoner_std['query_graph']['edges']['e0'] = { 'type':'gene_to_disease_association',
                                                       'source_id': 'n0',
                                                       'target_id': 'n1'
                                                     }
        # link disease to target
        reasoner_std['query_graph']['edges']['e1'] = { 'type':'disease_to_phenotype_association',
                                                       'source_id': 'n1',
                                                       'target_id': 'n2',
                                                     }
        handler = ReasonerStdHandler(source_ara='default',
                                     dict_query=reasoner_std)
        queries = handler.buildChpQueries()
        queries = handler.runChpQueries()
        reasoner_std_final = handler.constructDecoratedKG()
        KG = reasoner_std_final['message']['knowledge_graph']
        for edge_key in KG['edges'].keys():
            edge = KG['edges'][edge_key]
            if edge['type'] == 'disease_to_phenotype_association':
                p_survival = edge['has_confidence_level']
        print("probability of survival:",p_survival)

    # 10. check if the contribution property returns contributions
    def test_contribution_property(self):
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
        # empty response graph
        reasoner_std["results"] = [{ "node_bindings": dict(),
                                     "edge_bindings": dict()
                                  }]

        # add in evidence gene
        gene1 = ('RAF1', 'ENSEMBL:ENSG00000132155')
        reasoner_std['query_graph']['nodes']['n0'] = { 'type':'gene',
                                                       'curie':'{}'.format(gene1[1])
                                                     }
        # add in disease node
        disease = ('Breast_Cancer', 'MONDO:0007254')
        reasoner_std['query_graph']['nodes']['n1'] = { 'type':'disease',
                                                       'curie':'{}'.format(disease[1])
                                                     }
        # add target survival node
        phenotype = ('Survival_Time', 'EFO:0000714')
        reasoner_std['query_graph']['nodes']['n2'] = { 'type': 'phenotypicfeature',
                                                       'curie': '{}'.format(phenotype[1]),
                                                     }
        # link genes/drugs to disease
        reasoner_std['query_graph']['edges']['e0'] = { 'type':'gene_to_disease_association',
                                                       'source_id': 'n0',
                                                       'target_id': 'n1'
                                                     }
        # link disease to target
        reasoner_std['query_graph']['edges']['e1'] = { 'type':'disease_to_phenotype_association',
                                                       'source_id': 'n1',
                                                       'target_id': 'n2',
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
            if edge['type'] == 'disease_to_phenotype_association':
                p_survival = edge['has_confidence_level']
                report = edge['Description']
        print("probability of survival:",p_survival)

if __name__ == '__main__':
    unittest.main()
