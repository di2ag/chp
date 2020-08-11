import itertools
import tqdm
import numpy as np
import pickle
import os
import sys
from biothings_explorer.hint import Hint
from biothings_explorer.user_query_dispatcher import FindConnection
from biothings_explorer.export.reasoner import ReasonerConverter

from chp.query import Query
from chp.reasoner_std import ReasonerStdHandler

#import itertools
#import tqdm
#import numpy as np
#import pickle
#import os
#import sys
#from biothings_explorer.hint import Hint
#from biothings_explorer.user_query_dispatcher import FindConnection
#from biothings_explorer.export.reasoner import ReasonerConverter

#from chp.query import Query
#from chp.reasoner_std import ReasonerStdHandler

#query_path = '/home/ghyde/bkb-pathway-provider/tests/reasonerStdTests/sample_query1.pk'

import requests
import json

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
                                              'name':'{}'.format(gene[0]),
                                              'curie':'{}'.format(gene[1])
                                           })

# add in evidence drug
#drug = ('CYCLOPHOSPHAMIDE', 'CHEMBL:CHEMBL88')
#reasoner_std['query_graph']['nodes'].append({ 'id':'n{}'.format('0'),
#                                              'type':'Drug',
#                                              'name':'{}'.format(drug[0]),
#                                              'curie':'{}'.format(drug[1])
#                                           })

# add target survival node
reasoner_std['query_graph']['nodes'].append({ 'id': 'n{}'.format('1'),
                                              'type': 'Death',
                                              'curie': 'UBERON:0000071',
                                              'operator': '>=',
                                              'value': '1000'
                                           })

# link evidence to target survival node
reasoner_std['query_graph']['edges'].append({ 'id':'e{}'.format('0'),
                                              'type':'causes',
                                              'curie':['RO:0002410'],
                                              'source_id':'n{}'.format('0'),
                                              'target_id':'n{}'.format('1')
                                           })

handler = ReasonerStdHandler(source_ara='ranking',
                             dict_query=reasoner_std)

queries = handler.buildChpQueries()
queries = handler.runChpQueries()
reasoner_std_final = handler.constructDecoratedKG()

QG = reasoner_std_final['query_graph']
KG = reasoner_std_final['knowledge_graph']
res = reasoner_std_final['results']

# extract probability
KG_result_node = res['node_bindings'][0]['kg_id']
for node in KG['nodes']:
    if node['id'] == KG_result_node:
        p_survival = node['has_confidence_level']

# probability of surival given QG specification
print("Probability of survival > 1000 days is:", p_survival)



#reasoner_std['reasoner_id'] = 'ranking'
#payload = {'query': reasoner_std}
#r = requests.post('http://chp.thayer.dartmouth.edu/submitQuery/', json=payload)
#print(r.text)

