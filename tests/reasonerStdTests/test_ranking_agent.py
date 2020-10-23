import requests
import json
import unittest

from chp.reasoner_std import ReasonerStdHandler

# Tests chp/integrator/ranking_agent handler. Uses 5 test cases:
# 1. Normal gene only request
# 2. Normal drug only request
# 3. Tests with too much evidence (should return raised error)
# 4. Test with no evidence (should return raised error)
# 5. Normal gene only request without survival time specification - This
#    was at the request of Patrick Wang of Ranking Agent. Default survival 
#    time (if no survival is specified) is 970.

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
        handler = ReasonerStdHandler(source_ara='ranking',
                                     dict_query=reasoner_std)
        queries = handler.buildChpQueries()
        queries = handler.runChpQueries()
        reasoner_std_final = handler.constructDecoratedKG()

        KG = reasoner_std_final['knowledge_graph']
        res = reasoner_std_final['results']

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
        handler = ReasonerStdHandler(source_ara='ranking',
                                     dict_query=reasoner_std)
        queries = handler.buildChpQueries()
        queries = handler.runChpQueries()
        reasoner_std_final = handler.constructDecoratedKG()

        KG = reasoner_std_final['knowledge_graph']
        res = reasoner_std_final['results']

        # extract probability
        KG_result_edge = res['edge_bindings'][0]['kg_id']
        for edge in KG['edges']:
            if edge['id'] == KG_result_edge:
                p_survival = edge['has_confidence_level']
        print("probability of survival:",p_survival)

    # test with too much evidence
    def test_too_much_evidence(self):
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
            # add in evidence gene
            gene1 = ('RAF1', 'ENSEMBL:ENSG00000132155')
            reasoner_std['query_graph']['nodes'].append({ 'id':'n{}'.format('0'),
                                                          'type':'Gene',
                                                          'curie':'{}'.format(gene1[1])
                                                       })
            # add in second evidence gene
            gene2 = ('TTN', 'ENSEMBL:ENSG00000155657')
            reasoner_std['query_graph']['nodes'].append({ 'id':'n{}'.format('1'),
                                                          'type':'Gene',
                                                          'curie':'{}'.format(gene2[1])
                                                       })

            # add target survival node
            reasoner_std['query_graph']['nodes'].append({ 'id': 'n{}'.format('2'),
                                                          'type': 'PhenotypicFeature',
                                                          'curie': 'UBERON:0000071',
                                                       })
            # link evidence to target survival node
            reasoner_std['query_graph']['edges'].append({ 'id':'e{}'.format('0'),
                                                          'type':'causes',
                                                          'onset_qualifier':1000,
                                                          'source_id':'n{}'.format('0'),
                                                          'target_id':'n{}'.format('2')
                                                       })
            # link evidence to target survival node
            reasoner_std['query_graph']['edges'].append({ 'id':'e{}'.format('1'),
                                                          'type':'causes',
                                                          'onset_qualifier':1000,
                                                          'source_id':'n{}'.format('1'),
                                                          'target_id':'n{}'.format('2')
                                                       })
            handler = ReasonerStdHandler(source_ara='ranking',
                                         dict_query=reasoner_std)
            queries = handler.buildChpQueries()
        self.assertEqual(se.exception.code, 'More than 1 piece of evidence')

    # test with no evidence
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

            # add target survival node
            reasoner_std['query_graph']['nodes'].append({ 'id': 'n{}'.format('2'),
                                                          'type': 'PhenotypicFeature',
                                                          'curie': 'UBERON:0000071',
                                                       })
            handler = ReasonerStdHandler(source_ara='ranking',
                                         dict_query=reasoner_std)
            queries = handler.buildChpQueries()

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
        handler = ReasonerStdHandler(source_ara='ranking',
                                     dict_query=reasoner_std)
        queries = handler.buildChpQueries()
        queries = handler.runChpQueries()
        reasoner_std_final = handler.constructDecoratedKG()

        KG = reasoner_std_final['knowledge_graph']
        res = reasoner_std_final['results']

        # extract probability
        KG_result_edge = res['edge_bindings'][0]['kg_id']
        for edge in KG['edges']:
            if edge['id'] == KG_result_edge:
                p_survival = edge['has_confidence_level']
        print("probability of survival:",p_survival)

if __name__ == '__main__':
    unittest.main()

