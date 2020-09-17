import itertools
import tqdm
import numpy as np
import pickle
import os
#from biothings_explorer.hint import Hint
#from biothings_explorer.user_query_dispatcher import FindConnection

from chp.core.query import Query
from chp.core.reasoner_std import ReasonerStdHandler

query_path = '/home/ghyde/bkb-pathway-provider/tests/reasonerStdTests/sample_query.pk'

if os.path.exists(query_path):
    with open(query_path, 'rb') as f_:
        reasoner_std = pickle.load(f_)
else:
    ht = Hint()

    breast_cancer = ht.query('breast cancer')['Disease'][0]
    tamoxifin = ht.query('tamoxifen')['ChemicalSubstance'][0]

    fc = FindConnection(input_obj=breast_cancer, output_obj=tamoxifin, intermediate_nodes='Gene')
    fc.connect(verbose=True)

    reasoner_std = fc.to_reasoner_std()

    #-- Add target node and edge to Query Graph
    reasoner_std['query_graph']['nodes'].append({'id': 't1',
                                                 'type': 'PhenotypicFeature',
                                                 'curie': 'CHPDART:SURVIVAL',
                                                 'operator': '>=',
                                                 'value': '365'})
    reasoner_std['probability_targets'] = [('Survival', '>=', '365')]
    for i, node in enumerate(reasoner_std['query_graph']['nodes']):
        if node['id'][0] != 't':
            reasoner_std['query_graph']['edges'].append({'id': 'et{}'.format(i),
                                                         'source_id': node['id'],
                                                         'target_id': 't1'})

    with open('sample_query.pk', 'wb') as f_:
        pickle.dump(reasoner_std, f_)

handler = ReasonerStdHandler(source_ara='explorer',
                             dict_query=reasoner_std)

queries = handler.buildChpQueries()
queries = handler.runChpQueries()
reasoner_std_final = handler.constructDecoratedKG()

'''
#-- Get all options to fill query graph edges
options = dict()
for edge in reasoner_std['query_graph']['edges']:
    options[edge['id']] = []
for result in reasoner_std['results']['edge_bindings']:
    for qg_id, kg_ids in result.items():
        options[qg_id].append(kg_ids[0])

#-- Hash edges by ID
edges = {edge['id']: edge for edge in reasoner_std['knowledge_graph']['edges']}
nodes = {node['id']: node for node in reasoner_std['knowledge_graph']['nodes']}

#-- Get answers defined as a set of edges such that all nodes are consistent, i.e. complete subgraph that answers query.
answers = []
total = np.prod([len(option) for option in options.values()])
for answer in tqdm.tqdm(iterable=itertools.product(*list(options.values())), total=total):
    for i, kg_id in enumerate(answer):
        source_id = edges[kg_id]['source_id']
        target_id = edges[kg_id]['target_id']
        if i == 0:
            prev_target_id = target_id
            continue
        elif prev_target_id != source_id:
            break
        elif i == len(answer) - 1:
            answers.append(answer)

#-- Find target nodes
target_edges = []
for edge_id in options:
    if edge_id[0] == 't':
        target_edges.append(edge_id)

#-- Build out queries for each answer.
#-- Hash edges by ID
qnodes = {node['id']: node for node in reasoner_std['query_graph']['nodes']}
qedges = {edge['id']: edge for edge in reasoner_std['query_graph']['edges']}
queries = []
for answer in answers:
    #-- Get target
    targets = []
    meta_targets = []
    target_source_id = None
    for edge in target_edges:
        target_id = qedges[edge]['target_id']
        if target_source_id is None:
            target_source_id = qedges[edge]['source_id']
        elif target_source_id != qedges[edge]['source_id']:
            raise ValueError('Source of Probability target edge must be the same for each targets.')
        target_node = qnodes[target_id]
        ont, feature = target_node['curie'].split(':')
        if target_node['type'] == 'PhenotypicFeature':
            if ont == 'CHPDART':
                if feature == 'SURVIVAL':
                    meta_targets.append(('Survival_Time',
                                         target_node['operator'],
                                         target_node['value']))
                else:
                    raise NotImplementedError
            else:
                raise NotImplementedError
        else:
            raise NotImplementedError
    #-- Collect Evidence starting from target edges
    answer_edges_by_target_id = {edges[edge]['target_id']: edges[edge] for edge in answer}
    target_node = nodes[target_source_id]
    evidence = {}
    meta_evidence = []
    disease = None
    while True:
        #-- Add evidence
        if target_node['type'] == 'ChemicalSubstance':
            meta_evidence.append(('Patient_Drug(s)', '==', target_node['id']))
        elif target_node['type'] == 'Gene':
            evidence['_mut_{}'.format(target_node['name'])] = 'True'
        elif target_node['type'] == 'Disease':
            disease = target_node['name']
            if disease != 'BREAST CANCER':
                raise NotImplementedError
        else:
            raise NotImplementedError
        try:
            target_id = answer_edges_by_target_id[target_node['id']]['source_id']
        except KeyError:
            #-- If the current target node does not have an incoming edge
            break
        #-- Next target node is
        target_node = nodes[target_id]
    queries.append(Query(evidence=evidence,
                         meta_evidence=meta_evidence,
                         targets=targets,
                         meta_targets=meta_targets))


#TODO:
Now that we can read in reasoner_std query and extract answer subgraphs, let's parse query graph to determine
dataset to use (e.g. breast cancer), target probability to calculate (e.g. survival), and associated evidence.
'''
