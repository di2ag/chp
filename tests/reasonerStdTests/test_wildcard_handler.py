import requests
import json
import unittest

from chp.reasoner_std import ReasonerStdHandler

# Tests chp/integrator/wildcard_handler.py
# 1.

class testWildCardHandler(unittest.TestCase):

    # 1.
    def test_negative(self):
        # empty response
        reasoner_std = { "query_graph": {},
                         "knowledge_graph": {},
                         "response": {}
                       }
        # empty query graph
        reasoner_std["query_graph"] = { "edges": {},
                                        "nodes": {}
                                      }
        # empty knowledge graph
        reasoner_std["knowledge_graph"] = { "edges": {},
                                            "nodes": {}
                                          }
        # empty response graph
        reasoner_std["results"] = [{ "node_bindings": {},
                                    "edge_bindings": {}
                                  }]

        # add in evidence drug
        drug = ('CYCLOPHOSPHAMIDE', 'CHEMBL:CHEMBL88')
        reasoner_std['query_graph']['nodes']['n{}'.format('0')] = {
                                                      'type':'chemical_substance',
                                                      'curie':'{}'.format(drug[1])
        }

        # add in gene node (to be filled by contribution analysis
        reasoner_std['query_graph']['nodes']['n{}'.format('1')] = {
                                                      'type':'gene',
                                                   }

        #add in disease node
        disease = ('Breast_Cancer', 'MONDO:0007254')
        reasoner_std['query_graph']['nodes']['n{}'.format('2')] = {
                                                      'type':'disease',
                                                      'curie':'{}'.format(disease[1])
                                                   }

        # add target survival node
        phenotype = ('Survival_Time', 'EFO:0000714')
        reasoner_std['query_graph']['nodes']['n{}'.format('3')] = {
                                                      'type': 'phenotypic_feature',
                                                      'curie': '{}'.format(phenotype[1]),
                                                   }

        # link disease to target survival node
        reasoner_std['query_graph']['edges']['e{}'.format('0')] = {
                                                      'type':'disease_to_phenotypic_feature_association',
                                                      'source_id':'n{}'.format('2'),
                                                      'target_id':'n{}'.format('3'),
                                                     'properties': {
                                                         'qualifier': '>=',
                                                         'value': 940
                                                     }
                                                   }

        # link genes/drugs to disease
        reasoner_std['query_graph']['edges']['e{}'.format(1)] = {
                                                       'type':'gene_to_disease_association',
                                                       'source_id': 'n1',
                                                       'target_id': 'n2'
                                                     }
        reasoner_std['query_graph']['edges']['e{}'.format(2)] = {
                                                       'type':'chemical_to_disease_or_phenotypic_feature_association',
                                                       'source_id': 'n0',
                                                       'target_id': 'n2'
                                                     }


        json_formatted_str = json.dumps(reasoner_std, indent=2)
        print(json_formatted_str)

        handler = ReasonerStdHandler(source_ara='default',
                                     dict_query=reasoner_std,
                                     max_results=30)
        queries = handler.buildChpQueries()
        queries = handler.runChpQueries()
        reasoner_std_final = handler.constructDecoratedKG()

        KG = reasoner_std_final['knowledge_graph']
        res = reasoner_std_final['results']
        print(json.dumps(KG, indent=2))
        print(json.dumps(res, indent=2))

        # extract probability
        for _, edge in KG['edges'].items():
            if edge['type'] == 'disease_to_phenotypic_feature_association':
                p_survival = edge['has_confidence_level']
                #gene_contribs = edge['Description']
                break
        print("probability of survival:",p_survival)

        for qg_id, edge_bind in res[1]['edge_bindings'].items():
            kg_id = edge_bind["kg_id"]
            kg_edge = KG["edges"][kg_id]
            if kg_edge["type"] == 'gene_to_disease_association':
                weight = kg_edge["weight"]
                kg_gene_id = kg_edge["source_id"]
                kg_node = KG["nodes"][kg_gene_id]
                gene_info = (kg_node["name"], kg_gene_id, weight)
                break
        print(gene_info)

if __name__ == '__main__':
    unittest.main()
