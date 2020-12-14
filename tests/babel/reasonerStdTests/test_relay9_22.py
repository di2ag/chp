import requests
import json
import unittest

from chp.reasoner_std import ReasonerStdHandler

# Tests chp/integrator/relay9-22
# 1.

class testRelay9_22(unittest.TestCase):

    # 1.
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

        # add in evidence drug
        drug = ('CYCLOPHOSPHAMIDE', 'CHEMBL:CHEMBL88')
        reasoner_std['query_graph']['nodes'].append({ 'id':'n{}'.format('0'),
                                                      'type':'Drug',
                                                      'curie':'{}'.format(drug[1])
                                                   })

        # add target survival node
        phenotype = ('Survival_Time', 'EFO:0000714')
        reasoner_std['query_graph']['nodes'].append({ 'id': 'n{}'.format('1'),
                                                      'type': 'PhenotypicFeature',
                                                      'curie': '{}'.format(phenotype[1]),
                                                   })

        # link evidence to target survival node
        reasoner_std['query_graph']['edges'].append({ 'id':'e{}'.format('0'),
                                                      'type':'chemical_to_disease_or_phenotypic_feature_association',
                                                      'value':1000,
                                                      'source_id':'n{}'.format('0'),
                                                      'target_id':'n{}'.format('1')
                                                   })


        handler = ReasonerStdHandler(source_ara='relay9_22',
                                     dict_query=reasoner_std)
        queries = handler.buildChpQueries()
        queries = handler.runChpQueries()
        reasoner_std_final = handler.constructDecoratedKG()

        KG = reasoner_std_final['knowledge_graph']
        res = reasoner_std_final['results']

        # extract probability
        KG_result_edge = res['edge_bindings'][0]['kg_id']
        for edge in KG['edges']:
            if edge['type'] == 'chemical_to_disease_or_phenotypic_feature_association':
                p_survival = edge['has_confidence_level']
                gene_contribs = edge['Description']
        print("probability of survival:",p_survival)
        print(gene_contribs)

if __name__ == '__main__':
    unittest.main()
