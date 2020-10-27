import requests
import json
import unittest
import logging

from chp.reasoner_std import ReasonerStdHandler

# Tests chp/integrator/exploring_agent handler. Uses 8 test cases:
# 1. Testing negative result normalization
# 2. Normal request with two genes and a drug
# 3. Test with no evidence
# 4. Test with no evidence, but omitting KG and Results (should be handled by handler)
# 5. Test with one gene
# 6. Test with one drug

API_ADDRESS = 'http://127.0.0.1:8000'

logger = logging.getLogger(__name__)
logger.warning('Make sure your django server is running at: {}.'.format(API_ADDRESS))

def submitQuery(query):
    payload = {'query': query}

    r = requests.post('{}/submitQuery/'.format(API_ADDRESS), json=payload)
    chp_res = json.loads(r.content)
    return chp_res

class testExploringAgent(unittest.TestCase):

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

        # set ara source
        reasoner_std['reasoner_id'] = 'exploring'

        # submit query
        chp_res = submitQuery(reasoner_std)

        QG = chp_res['query_graph']
        KG = chp_res['knowledge_graph']
        res = chp_res['results']

        # extract probability
        KG_result_edge = res['edge_bindings'][0]['kg_id']
        for edge in KG['edges']:
            if edge['id'] == KG_result_edge:
                p_survival = edge['has_confidence_level']

        # probability of surival given QG specification
        print("Probability of survival > 1000 days is:", p_survival)

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
        gene1 = ('RAF1', 'ENSEMBL:ENSG00000132155')
        reasoner_std['query_graph']['nodes'].append({ 'id':'n{}'.format('0'),
                                                      'type':'Gene',
                                                      'curie':'{}'.format(gene1[1])
                                                   })
        # add in second evidence gene
        gene2 = ('BRCA1', 'ENSEMBL:ENSG00000012048')
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
                                                      'value':970,
                                                      'source_id':'n{}'.format('3'),
                                                      'target_id':'n{}'.format('4')
                                                   })
        # set ara source
        reasoner_std['reasoner_id'] = 'exploring'

        # submit query
        chp_res = submitQuery(reasoner_std)

        QG = chp_res['query_graph']
        KG = chp_res['knowledge_graph']
        res = chp_res['results']

        for edge in KG['edges']:
            if edge['type'] == 'disease_to_phenotype_association':
                p_survival = edge['has_confidence_level']
        print("probability of survival:",p_survival)

    # 3. Test with no evidence
    def test_no_evidence(self):
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
        # set ara source
        reasoner_std['reasoner_id'] = 'exploring'

        # submit query
        chp_res = submitQuery(reasoner_std)

        QG = chp_res['query_graph']
        KG = chp_res['knowledge_graph']
        res = chp_res['results']

        for edge in KG['edges']:
            if edge['type'] == 'disease_to_phenotype_association':
                p_survival = edge['has_confidence_level']
        print("probability of survival:",p_survival)

    # 4. Test with no evidence, but omitting KG and Results (should be handled by handler)
    def test_no_evidence_omitting_KG_and_results(self):
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
        # set ara source
        reasoner_std['reasoner_id'] = 'exploring'

        # submit query
        chp_res = submitQuery(reasoner_std)

        QG = chp_res['query_graph']
        KG = chp_res['knowledge_graph']
        res = chp_res['results']

        for edge in KG['edges']:
            if edge['type'] == 'disease_to_phenotype_association':
                p_survival = edge['has_confidence_level']
        print("probability of survival:",p_survival)

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
        # set ara source
        reasoner_std['reasoner_id'] = 'exploring'

        # submit query
        chp_res = submitQuery(reasoner_std)

        QG = chp_res['query_graph']
        KG = chp_res['knowledge_graph']
        res = chp_res['results']

        for edge in KG['edges']:
            if edge['type'] == 'disease_to_phenotype_association':
                p_survival = edge['has_confidence_level']
        print("probability of survival:",p_survival)

    # 6. Test with one drug
    def test_one_drug(self):
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
        # set ara source
        reasoner_std['reasoner_id'] = 'exploring'

        # submit query
        chp_res = submitQuery(reasoner_std)

        QG = chp_res['query_graph']
        KG = chp_res['knowledge_graph']
        res = chp_res['results']

        for edge in KG['edges']:
            if edge['type'] == 'disease_to_phenotype_association':
                p_survival = edge['has_confidence_level']
        print("probability of survival:",p_survival)

class testRankingAgent(unittest.TestCase):

    # 1. Normal gene only request
    def test_gene_request(self):

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
        # add in evidence genes
        gene = ('RAF1', 'ENSEMBL:ENSG00000132155')
        reasoner_std['query_graph']['nodes'].append({ 'id':'n{}'.format('0'),
                                                      'type':'Gene',
                                                      'curie':'{}'.format(gene[1])
                                                   })
        # add target survival node
        reasoner_std['query_graph']['nodes'].append({ 'id': 'n{}'.format('1'),
                                                      'type': 'PhenotypicFeature',
                                                      'curie': 'UBERON:0000071',
                                                   })
        # link evidence to target survival node
        reasoner_std['query_graph']['edges'].append({ 'id':'e{}'.format('0'),
                                                      'type':'causes',
                                                      'onset_qualifier':1000,
                                                      'source_id':'n{}'.format('0'),
                                                      'target_id':'n{}'.format('1')
                                                   })
        # set ara source
        reasoner_std['reasoner_id'] = 'ranking'

        # submit query
        chp_res = submitQuery(reasoner_std)

        QG = chp_res['query_graph']
        KG = chp_res['knowledge_graph']
        res = chp_res['results']

        # extract probability
        KG_result_edge = res['edge_bindings'][0]['kg_id']
        for edge in KG['edges']:
            if edge['id'] == KG_result_edge:
                p_survival = edge['has_confidence_level']
        print("probability of survival:",p_survival)

    # Normal drug only request
    def test_drug_request(self):
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
        # add target survival node
        reasoner_std['query_graph']['nodes'].append({ 'id': 'n{}'.format('1'),
                                                      'type': 'PhenotypicFeature',
                                                      'curie': 'UBERON:0000071',
                                                   })
        # link evidence to target survival node
        reasoner_std['query_graph']['edges'].append({ 'id':'e{}'.format('0'),
                                                      'type':'causes',
                                                      'onset_qualifier':1000,
                                                      'source_id':'n{}'.format('0'),
                                                      'target_id':'n{}'.format('1')
                                                   })
        # set ara source
        reasoner_std['reasoner_id'] = 'ranking'

        # submit query
        chp_res = submitQuery(reasoner_std)

        QG = chp_res['query_graph']
        KG = chp_res['knowledge_graph']
        res = chp_res['results']

        # extract probability
        KG_result_edge = res['edge_bindings'][0]['kg_id']
        for edge in KG['edges']:
            if edge['id'] == KG_result_edge:
                p_survival = edge['has_confidence_level']
        print("probability of survival:",p_survival)


    # Normal gene only request without survival time specification
    def test_use_default_survival_request(self):

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
        # add in evidence genes
        gene = ('RAF1', 'ENSEMBL:ENSG00000132155')
        reasoner_std['query_graph']['nodes'].append({ 'id':'n{}'.format('0'),
                                                      'type':'Gene',
                                                      'curie':'{}'.format(gene[1])
                                                   })
        # add target survival node
        reasoner_std['query_graph']['nodes'].append({ 'id': 'n{}'.format('1'),
                                                      'type': 'PhenotypicFeature',
                                                      'curie': 'UBERON:0000071',
                                                   })
        # link evidence to target survival node
        reasoner_std['query_graph']['edges'].append({ 'id':'e{}'.format('0'),
                                                      'type':'causes',
                                                      'source_id':'n{}'.format('0'),
                                                      'target_id':'n{}'.format('1')
                                                   })
        # set ara source
        reasoner_std['reasoner_id'] = 'ranking'

        # submit query
        chp_res = submitQuery(reasoner_std)

        QG = chp_res['query_graph']
        KG = chp_res['knowledge_graph']
        res = chp_res['results']

        # extract probability
        KG_result_edge = res['edge_bindings'][0]['kg_id']
        for edge in KG['edges']:
            if edge['id'] == KG_result_edge:
                p_survival = edge['has_confidence_level']
        print("probability of survival:",p_survival)

if __name__ == '__main__':
    unittest.main()
