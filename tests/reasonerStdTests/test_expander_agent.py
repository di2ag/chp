import requests
import json
import unittest

from chp.reasoner_std import ReasonerStdHandler

# Tests chp/integrator/expander_agent.py
# 1.

class testExpanderAgent(unittest.TestCase):

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
                                                      'source_id':'n{}'.format('0'),
                                                      'target_id':'n{}'.format('1')
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

        handler = ReasonerStdHandler(source_ara='expander',
                                     dict_query=reasoner_std)
        queries = handler.buildChpQueries()
        queries = handler.runChpQueries()
        reasoner_std_final = handler.constructDecoratedKG()

        KG = reasoner_std_final['knowledge_graph']
        res = reasoner_std_final['results']

        # extract probability
        for edge in KG['edges']:
            if edge['type'] == 'disease_to_phenotypic_association':
                p_survival = edge['has_confidence_level']
                #gene_contribs = edge['Description']
                break
        print("probability of survival:",p_survival)
        #print(gene_contribs)

        print(json.dumps(KG['edges'], indent=2))
        # extract a relative contribution
        for node_bind in res['node_bindings'][0]:
            if node_bind['qg_id'] == 'n1':
                print('Got here')
                for node in KG['nodes']:
                    if node['id'] == node_bind['kg_id']:
                        top_gene = [node['name'], node['curie']]
                        break
                for edge in KG['edges']:
                    if edge['source_id'] == node_bind['kg_id']:
                        top_gene.append(edge['weight'])
        print(top_gene)

if __name__ == '__main__':
    unittest.main()
