import os
import sys
import pickle
from operator import ge, le, eq
import itertools
import tqdm
import copy
from concurrent.futures import ProcessPoolExecutor, wait
import time

from pybkb.core.cpp_base.reasoning import revision as cpp_revision
from pybkb.core.cpp_base.reasoning import updating as cpp_updating
from pybkb.core.python_base.reasoning import updating as py_updating
from pybkb.core.common.bayesianKnowledgeBase import BKB_component, BKB_I_node, BKB_S_node

class Reasoner:
    def __init__(self, fused_bkb, patients, cpp_reasoning=False):
        self.fused_bkb = fused_bkb
        self.patients = patients
        self.cpp_reasoning = cpp_reasoning

        #-- Preprocess src hash values
        src_component_indices = fused_bkb.getSrcComponents()
        src_hashs = dict()
        for comp_idx in src_component_indices:
            for state_idx in range(fused_bkb.getNumberComponentINodes(comp_idx)):
                state_name = fused_bkb.getComponentINodeName(comp_idx, state_idx)
                #-- Collect src numbers and string name or hash
                src_num = int(state_name.split('_')[0][1:-1])
                try:
                    src_hashs[src_num] = int(''.join(state_name.split('_')[1:]))
                except ValueError:
                    continue
                    #src_hashs[src_num] = state_name
        #self.src_hashs = set(src_hashs)
        self.src_hashs = src_hashs

    #-- Should be source hash value followed by a dictionary of all available meta data. The file is assumed to be a pickle.
    def set_src_metadata(self, metadata_file):
        with open(metadata_file, 'rb') as m_:
            self.metadata = pickle.load(m_)
        #print('Available Metadata:')
        self.metadata_labels = list()
        self.metadata_ranges = dict()
        for hash_key in self.metadata:
            self.metadata_labels.extend(list(self.metadata[hash_key].keys()))
            for label in self.metadata[hash_key]:
                try:
                    val = float(self.metadata[hash_key][label])
                    if label in self.metadata_ranges:
                        min_, max_ = self.metadata_ranges[label]
                        if val < min_:
                            min_ = val
                        if val > max_:
                            max_ = val
                        self.metadata_ranges[label] = (min_, max_)
                    else:
                        self.metadata_ranges[label] = (val, val)
                except:
                    if label in self.metadata_ranges:
                        self.metadata_ranges[label].add(self.metadata[hash_key][label])
                    else:
                        self.metadata_ranges[label] = set([self.metadata[hash_key][label]])
        self.metadata_labels = list(set(self.metadata_labels))
        #for meta in set(self.metadata_labels):
        #    print('\t{}'.format(meta))

    def solve_query(self, query):
        print('Reasoning...')
        start_time = time.time()
        if self.cpp_reasoning:
            if query.type == 'revision':
                res = cpp_revision(query.bkb,
                               query.evidence,
                               marginal_evidence=query.marginal_evidence,
                               targets=query.targets,
                               file_prefix=query.name)
            elif query.type == 'updating':
                res = cpp_updating(query.bkb,
                               query.evidence,
                               marginal_evidence=query.marginal_evidence,
                               targets=query.targets,
                               file_prefix=query.name)
            else:
                raise ValueError('Unreconginzed reasoning type: {}.'.format(query.type))
        else:
            if query.type == 'revision':
                raise NotImplementedError('Python Revision is not currently implemented.')
            elif query.type == 'updating':
                res = py_updating(query.bkb,
                               query.evidence,
                               query.targets)
            else:
                raise ValueError('Unreconginzed reasoning type: {}.'.format(query.type))
        compute_time = time.time() - start_time
        query.result = res
        query.compute_time = compute_time
        return query

    def process_metaVariables(self, bkb, meta_variables):
        #-- Get all src components
        #src_components = bkb.getSrcComponents()

        #-- Collect source hashs that match metadata as well as population stats
        transformed_meta = dict()
        pop_stats = dict()
        for i, meta in enumerate(meta_variables):
            bkb, transformed_meta_, matched_srcs = _addDemographicOption(meta, bkb, self.src_hashs, self.metadata, option_dependencies=meta_variables[:i])
            transformed_meta.update(transformed_meta_)

        #-- Process Sources
        #-- If first piece of evidence connect to all patients.
        comp_idx = bkb.getComponentIndex('{} {} {}'.format(meta[0], meta[1], meta[2]))
        inode_true_idx = bkb.getComponentINodeIndex(comp_idx, 'True')
        inode_false_idx = bkb.getComponentINodeIndex(comp_idx, 'False')
        bkb = _addSrcConnections(comp_idx, inode_true_idx, inode_false_idx, bkb, matched_srcs, self.src_hashs)
        return transformed_meta, bkb

    def analyze_query(self, query):
        #-- Copy BKB
        bkb = copy.deepcopy(self.fused_bkb)

        meta_variables = []
        if query.meta_evidence is not None:
            meta_variables += query.meta_evidence

        if query.meta_targets is not None:
            meta_variables += query.meta_targets
            transformed_meta_targets = ['{} {} {}'.format(target[0], target[1], target[2]) for target in query.meta_targets]

        transformed_meta, bkb = self.process_metaVariables(bkb, meta_variables)
        #-- Collect Evidence
        transformed_meta_evidence = dict()
        for ev in query.meta_evidence:
            meta_name = '{} {} {}'.format(ev[0], ev[1], ev[2])
            transformed_meta_evidence[meta_name] = transformed_meta[meta_name]

        #-- Conduct Reasoning
        query.evidence.update(transformed_meta_evidence)
        query.targets.extend(transformed_meta_targets)
        query.bkb = bkb

        bkb.save('bkb-test.bkb')

        return self.solve_query(query)

def _process_operator(op):
    if op == '>=':
        return ge
    elif op == '<=':
        return le
    elif op == '==':
        return eq
    else:
        raise ValueError('Unknown Operator')

def _processDependencyHead(head, bkb):
    head_ev, state = head
    prop_head, op_str_head, val_head = head_ev
    op_head = _process_operator(op_str_head)

    comp_head_idx = bkb.getComponentIndex('{} {} {}'.format(prop_head, op_str_head, val_head))
    if state:
        i_node_head_idx = bkb.getComponentINodeIndex(comp_head_idx, 'True')
    else:
        i_node_head_idx = bkb.getComponentINodeIndex(comp_head_idx, 'False')

    return comp_head_idx, i_node_head_idx

def _processDependecyTail(tail, bkb):
    processed_tail = list()
    for tail_ in tail:
        tail_ev, state = tail_
        prop_tail, op_str_tail, val_tail = tail_ev
        op_tail = _process_operator(op_str_tail)

        comp_tail_idx = bkb.getComponentIndex('{} {} {}'.format(prop_tail, op_str_tail, val_tail))
        if state:
            i_node_tail_idx = bkb.getComponentINodeIndex(comp_tail_idx, 'True')
        else:
            i_node_tail_idx = bkb.getComponentINodeIndex(comp_tail_idx, 'False')
        processed_tail.append((comp_tail_idx, i_node_tail_idx))
    return processed_tail

def _processOptionDependency(option, option_dependencies, bkb, src_population, src_population_data):
    #-- Get consistent option combinations
    tail_product = [combo for combo in itertools.product(option_dependencies, [True, False])]
    tail_combos = list()
    for combo in itertools.combinations(tail_product, r=len(option_dependencies)):
        combo = list(combo)
        prop_set = set()
        for ev_state in combo:
            ev_, state = ev_state
            prop_, op_, val_ = ev_
            prop_set.add(prop_)
        if len(prop_set) == len(combo):
            tail_combos.append(combo)
    combos = [combo for combo in itertools.product([(option, True), (option, False)], tail_combos)]

    #-- Calculate joint probabilities from data
    counts = list()
    for combo in combos:
        head, tail = combo
        #-- Put head and tail in one list
        combo = [head] + tail
        count = 0
        for entity_name, src_name in src_population.items():
            truth = list()
            for ev_state in combo:
                ev_, state = ev_state
                prop_, op_str_, val_ = ev_
                op_ = _process_operator(op_str_)
                res = op_(src_population_data[src_name][prop_], val_)
                truth.append(res == state)
            if all(truth):
                count += 1
        counts.append(count)
    probs = [float(count / len(src_population)) for count in counts]

    #-- Setup each S-node
    for j, combo in enumerate(combos):
        head, tail = combo
        #-- Process head
        comp_head_idx, i_node_head_idx = _processDependencyHead(head, bkb)
        #-- Process Tail
        processed_tail = _processDependecyTail(tail, bkb)
        #-- Add Snode
        if probs[j] > 0:
            bkb.addSNode(BKB_S_node(init_component_index=comp_head_idx, init_state_index=i_node_head_idx, init_probability=probs[j], init_tail=processed_tail))

    return bkb

def _addDemographicOption(option, bkb, src_population, src_population_data, option_dependencies=list()):
    prop, op_str, val = option
    op = _process_operator(op_str)
    matched_srcs = set()
    pop_count_true = 0
    for entity_name, src_name in src_population.items():
        if op(src_population_data[src_name][prop], val):
            matched_srcs.add(entity_name)
            pop_count_true += 1
    prob = float(pop_count_true / len(src_population))

    #comp = BKB_component('{} {} {}'.format(prop, op_str, val))
    #inode_true = BKB_I_node(init_name='True', init_component=comp)
    #inode_false = BKB_I_node(init_name='False', init_component=comp)
    
    comp_idx = bkb.addComponent('{} {} {}'.format(prop, op_str, val))
    inode_true_idx = bkb.addComponentState(comp_idx, 'True')
    inode_false_idx = bkb.addComponentState(comp_idx, 'False')

    #-- Create option dictionary
    options_dict = {bkb.getComponentName(comp_idx): bkb.getComponentINodeName(comp_idx, inode_true_idx)}

    #-- If no chain rule dependencies, just add prior s-nodes
    if len(option_dependencies) == 0:
        snode_1 = BKB_S_node(init_component_index=comp_idx, init_state_index=inode_true_idx, init_probability=prob)
        snode_2 = BKB_S_node(init_component_index=comp_idx, init_state_index=inode_false_idx, init_probability=1-prob)
        bkb.addSNode(snode_1)
        bkb.addSNode(snode_2)

    #-- Process Dependencies
    else:
        bkb = _processOptionDependency(option, option_dependencies, bkb, src_population, src_population_data)

    return bkb, options_dict, matched_srcs

def _constructSNodesByHead(bkb):
    S_nodes_by_head = dict()
    for snode in bkb.getAllSNodes():
        (cidx, iidx) = snode.getHead()
        try:
            S_nodes_by_head[cidx]
        except KeyError:
            S_nodes_by_head[cidx] = dict()
        try:
            S_nodes_by_head[cidx][iidx].append(snode)
        except KeyError:
            S_nodes_by_head[cidx][iidx] = list()
            S_nodes_by_head[cidx][iidx].append(snode)
    return S_nodes_by_head

def _addSrcConnections(comp_idx, inode_true_idx, inode_false_idx, bkb, matched_srcs, src_population):
    #-- Attach to source nodes
    src_component_indices = bkb.getSrcComponents()
    S_nodes_by_head = _constructSNodesByHead(bkb)

    for entity_name, src_name in src_population.items():
        for src_comp_idx in src_component_indices:
            try:
                src_state_idx = bkb.findINode(src_comp_idx, str(src_name), contains=True)

                #-- Find Prior snode, capture probability and remove.
                #cidx = bkb.getComponentIndex(src_comp)
                #iidx = src_comp.getStateIndex(src_state)
                s_node = bkb.S_nodes_by_head[src_comp_idx][src_state_idx][0]
            except AttributeError:
                continue
            prob = s_node.probability
            bkb.removeSNode(s_node)

            if entity_name in matched_srcs:
                bkb.addSNode(BKB_S_node(init_component_index=src_comp_idx,
                                              init_state_index=src_state_idx,
                                              init_probability=prob,
                                              init_tail = [(comp_idx, inode_true_idx)]))
            else:
                bkb.addSNode(BKB_S_node(init_component_index=src_comp_idx,
                                              init_state_index=src_state_idx,
                                              init_probability=prob,
                                              init_tail = [(comp_idx, inode_false_idx)]))
    return bkb
