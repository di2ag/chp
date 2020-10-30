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
                                                      'type':'drug',
                                                      'curie':'{}'.format(drug[1])
                                                   })

        # add in gene node (to be filled by contribution analysis
        reasoner_std['query_graph']['nodes'].append({ 'id':'n{}'.format('1'),
                                                      'type':'gene',
                                                   })

        #add in disease node
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

        # link disease to target survival node
        reasoner_std['query_graph']['edges'].append({ 'id':'e{}'.format('0'),
                                                      'type':'disease_to_phenotypic_association',
                                                      'value':1000,
                                                      'source_id':'n{}'.format('2'),
                                                      'target_id':'n{}'.format('3')
                                                   })

        # link genes/drugs to disease
        reasoner_std['query_graph']['edges'].append({ 'id': 'e{}'.format(1),
                                                       'type':'gene_to_disease_association',
                                                       'source_id': 'n1',
                                                       'target_id': 'n2'
                                                     })
        reasoner_std['query_graph']['edges'].append({  'id': 'e{}'.format(2),
                                                       'type':'chemical_to_disease_or_phenotypic_feature_association',
                                                       'source_id': 'n0',
                                                       'target_id': 'n2'
                                                     })


        json_formatted_str = json.dumps(reasoner_std, indent=2)
        print(json_formatted_str)

        handler = ReasonerStdHandler(source_ara='default',
                                     dict_query=reasoner_std)
        queries = handler.buildChpQueries()
        queries = handler.runChpQueries()
        reasoner_std_final = handler.constructDecoratedKG()

        KG = reasoner_std_final['knowledge_graph']
        res = reasoner_std_final['results']
        print(json.dumps(KG, indent=2))
        print(json.dumps(res, indent=2))

        # extract probability
        for edge in KG['edges']:
            if edge['type'] == 'disease_to_phenotypic_association':
                p_survival = edge['has_confidence_level']
                #gene_contribs = edge['Description']
                break
        print("probability of survival:",p_survival)
        #print(gene_contribs)

        #print(json.dumps(KG['edges'], indent=2))
        # extract a relative contribution
        # create maps
        qg_edge_map = {edge["id"]: edge for edge in reasoner_std_final['query_graph']['edges']}
        qg_node_map = {node["id"]: node for node in reasoner_std_final['query_graph']['nodes']}
        kg_edge_map = {edge["id"]: edge for edge in reasoner_std_final['knowledge_graph']['edges']}
        kg_node_map = {node["id"]: node for node in reasoner_std_final['knowledge_graph']['nodes']}

        #print(res["edge_bindings"])
        #print(len(res["edge_bindings"]))

        for edge_bind in res['edge_bindings']:
            qg_id = edge_bind["qg_id"]
            kg_id = edge_bind["kg_id"]
            kg_edge = kg_edge_map[kg_id]
            if kg_edge["type"] == 'gene_to_disease_association':
                weight = kg_edge["weight"]
                kg_gene_id = kg_edge["source_id"]
                kg_node = kg_node_map[kg_gene_id]
                gene_info = (kg_node["name"], kg_node["curie"], weight)
                break
        print(gene_info)

if __name__ == '__main__':
    unittest.main()
