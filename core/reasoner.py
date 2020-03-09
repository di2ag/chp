import os
import sys
import math
import pickle
from operator import ge, le, eq
import itertools
import tqdm
import copy
from concurrent.futures import ProcessPoolExecutor, wait
import time

#sys.path.append('/home/cyakaboski/src/python/projects/bkb-pathway-provider/core')
sys.path.append('/home/ghyde/bkb-pathway-provider/core')
from query import Query

from pybkb.core.cpp_base.reasoning import revision as cpp_revision
from pybkb.core.cpp_base.reasoning import updating as cpp_updating
from pybkb.core.python_base.reasoning import updating as py_updating
from pybkb.core.python_base.fusion_collapse import collapse_sources
from pybkb.core.common.bayesianKnowledgeBase import BKB_S_node

class Reasoner:
    def __init__(self, fused_bkb, patient_data=None, cpp_reasoning=False):
        self.fused_bkb = fused_bkb
        self.metadata = patient_data
        self.cpp_reasoning = cpp_reasoning
        self.metadata = patient_data
        if patient_data is not None:
            try:
                self.setup()
            except:
                raise ValueError('Likely patient data is not in the right format.')

    def setup(self):
        #-- Collapse sources in fused bkb
        self.collapsed_bkb = collapse_sources(self.fused_bkb)

        #-- Preprocess src hash values
        src_component_indices = self.fused_bkb.getSrcComponents()
        src_hashs = dict()
        for comp_idx in src_component_indices:
            for state_idx in range(self.fused_bkb.getNumberComponentINodes(comp_idx)):
                state_name = self.fused_bkb.getComponentINodeName(comp_idx, state_idx)
                #-- Collect src numbers and string name or hash
                src_num = int(state_name.split('_')[0][1:-1])
                try:
                    src_hashs[src_num] = int(''.join(state_name.split('_')[1:]))
                except ValueError:
                    continue
                    #src_hashs[src_num] = state_name
        #self.src_hashs = set(src_hashs)
        self.src_hashs = src_hashs

        #-- Preprocess Patient Variant data into dictionary
        for src_hash, data_dict in self.metadata.items():
            if 'Patient_Gene_Variants' in data_dict:
                #-- Transform into directionary
                gene_variant_dict = dict()
                for gene_variant_str in data_dict['Patient_Gene_Variants']:
                    gene_variant_split = gene_variant_str.split('-')
                    gene = gene_variant_split[0]
                    variant = gene_variant_split[1]
                    gene_variant_dict[gene] = variant
                data_dict['Patient_Gene_Variant'] = gene_variant_dict

        #-- Setup Metadata labels and ranges
        self.metadata_labels = list()
        self.metadata_ranges = dict()
        for hash_key in self.metadata:
            self.metadata_labels.extend(list(self.metadata[hash_key].keys()))
            for label in self.metadata[hash_key]:
                val = self.metadata[hash_key][label]
                if type(val) == str or type(val) == int or type(val) == float:
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
                    except ValueError:
                        if label in self.metadata_ranges:
                            self.metadata_ranges[label].add(val)
                        else:
                            self.metadata_ranges[label] = set([val])
                elif type(val) == list or type(val) == tuple:
                    if label in self.metadata_ranges:
                        self.metadata_ranges[label].update(set(val))
                    else:
                        self.metadata_ranges[label] = set(val)
                elif type(val) == dict:
                    #-- collapse dict into a set
                    if label in self.metadata_ranges:
                        self.metadata_ranges[label].update(set(self.metadata[hash_key][label]))
                    else:
                        self.metadata_ranges[label] = set(list(self.metadata[hash_key][label].items()))
                else:
                    raise TypeError('Unrecongized data type: {},\n{}'.format(type(val), val))
        self.metadata_labels = list(set(self.metadata_labels))

    #-- Should be source hash value followed by a dictionary of all available meta data. The file is assumed to be a pickle.
    def set_src_metadata(self, metadata_file):
        with open(metadata_file, 'rb') as m_:
            self.metadata = pickle.load(m_)
        #-- Run setup
        self.setup()

    def solve_query(self, query, target_strategy='chain', interpolation='independence'):
        #-- I don't think we need to make a copy of the query.
        #query = copy.deepcopy(query)
        #print('Reasoning...')
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

        #-- If target strategy is topological multiply all the target probabilities together.
        if target_strategy == 'topological':
            updates = query.result.process_updates()
            #-- Convert to actual meta target probability using independence assumption.
            meta_target_updates = dict()
            for comp_name, state_dict in updates.items():
                #-- There should always be two states, i.e. True and False. If theres only one just omit it.
                if len(state_dict) < 2:
                    continue
                meta_target_name = comp_name.split('_')[0]
                if meta_target_name not in meta_target_updates:
                    meta_target_updates[meta_target_name] = state_dict
                else:
                    for state_name, log_prob in state_dict.items():
                        meta_target_updates[meta_target_name][state_name] += log_prob
            #-- Convert to probabilities and Normalize
            try:
                for meta_target, state_dict in meta_target_updates.items():
                    total_prob = 0
                    for target_state, log_prob in state_dict.items():
                        total_prob += math.exp(log_prob)
                    new_state_dict = dict()
                    for target_state, log_prob in state_dict.items():
                        new_state_dict[target_state] = float(math.exp(log_prob) / total_prob)
                    meta_target_updates[meta_target] = new_state_dict
            except OverflowError:
                print('Warning: Encountered an overflow error, so leaving as log probabilities')
            query.result.meta_target_updates = meta_target_updates
        return query

    def process_metaVariables(self, bkb, meta_evidence, meta_targets, target_strategy='chain'):
        if target_strategy == 'chain':
            meta_variables = []
            if meta_evidence is not None:
                meta_variables += meta_evidence
            if meta_targets is not None:
                meta_variables += meta_targets
        elif target_strategy == 'topological':
            if len(meta_targets) > 1:
                raise NotImplementedError('Only one target can be specified and must be a demographic target.')
            meta_variables = meta_evidence
        else:
            raise ValueError('Target Strategy must be chain or topological.')

        #-- Collect source hashs that match metadata as well as population stats
        transformed_meta = dict()
        pop_stats = dict()
        if meta_variables is not None:
            for i, meta in enumerate(meta_variables):
                bkb, transformed_meta_, matched_srcs = _addDemographicOption(meta, bkb, self.src_hashs, self.metadata, option_dependencies=meta_variables[:i])
                transformed_meta.update(transformed_meta_)

            #-- Process Sources
            #-- If first piece of evidence connect to all patients.
            comp_idx = bkb.getComponentIndex('{} {} {}'.format(meta[0], meta[1], meta[2]))
            inode_true_idx = bkb.getComponentINodeIndex(comp_idx, 'True')
            inode_false_idx = bkb.getComponentINodeIndex(comp_idx, 'False')
            bkb = _addSrcConnections(comp_idx, inode_true_idx, inode_false_idx, bkb, matched_srcs, self.src_hashs)

        #-- If topological target stradegy connect target to each bottom I node.
        if target_strategy == 'topological':
            #-- Construct link and stats from bottom I-nodes to individual meta targets.
            bkb, transformed_meta_ = _addTargetToLastTopologVariables(meta_targets[0], bkb, self.src_hashs, self.metadata)
            transformed_meta.update(transformed_meta_)
        return transformed_meta, bkb

    def solve_query_independence(self, query, genetic_evidence, target_strategy):
        bkb = query.bkb
        queries = list()
        non_gene_evidence = copy.deepcopy(query.evidence)
        for gene_key in genetic_evidence:
            del non_gene_evidence[gene_key]
        for genetic_comp, genetic_state in genetic_evidence.items():
            #-- Create individual gene evidence
            independ_evidence = copy.deepcopy(non_gene_evidence)
            independ_evidence[genetic_comp] = genetic_state
            q = Query(evidence=independ_evidence,
                      targets=query.targets,
                      type='updating')
            q.bkb = bkb
            queries.append(q)

        finished_queries = list()
        for q_ in queries:
            finished_queries.append(self.solve_query(q_, target_strategy=target_strategy))

        #-- Run all independent queries in parrellel
        #with ProcessPoolExecutor(max_workers=15) as executor:
        #    #finished_queries = list()
        #    #for q_ in executor.map(self.solve_query, [(q__, target_strategy) for q__ in queries]):
        #    #    finished_queries.append(q_)
        #    futures = list()
        #    for q_ in queries:
        #        futures.append(executor.submit(self.solve_query, q_, target_strategy))
        #    finished_queries = [q_.result() for q_ in futures]

        #-- Collect queries and calculate independent probs
        res = dict()
        for q_ in finished_queries:
            updates = q_.result.process_updates()
            for comp_name, state_dict in updates.items():
                state_keys = state_dict.keys()
                state_probs = state_dict.values()
                # skip -1 where there is no inference on a particular state_key
                if min(state_probs) <= 0:
                    continue
                if len(state_dict) < 2:
                    continue
                if comp_name in res:
                    for state_name, prob in state_dict.items():
                        #-- If probability is greater than one than its a log prob
                        if prob > 1:
                            res[comp_name][state_name] += prob
                        else:
                            res[comp_name][state_name] *= prob
                else:
                    res[comp_name] = state_dict
        query.independ_queries = queries
        query.independ_result = res
        return query

    '''
    Target strategy can be chain, where targets are located in the chain rule (legacy) or
    topological where target is attached to all nodes at the bottom of the topological sort. This
    methodolgy can only except on target then.
    '''
    def analyze_query(self, query, save_dir=None, preprocessed_bkb=None, target_strategy='chain', interpolation='standard'):
        #-- Duplicate Gene Evidence
        genetic_evidence = copy.deepcopy(query.evidence)
        if preprocessed_bkb is not None:
            #bkb = copy.deepcopy(preprocessed_bkb)

            if query.meta_targets is not None:
                transformed_meta_targets = ['{} {} {}'.format(target[0], target[1], target[2]) for target in query.meta_targets]

            if query.meta_evidence is not None:
                transformed_meta_evidence = dict()
                for ev in query.meta_evidence:
                    meta_name = '{} {} {}'.format(ev[0], ev[1], ev[2])
                    transformed_meta_evidence[meta_name] = 'True'
                    query.evidence.update(transformed_meta_evidence)

            query.targets.extend(transformed_meta_targets)
            query.bkb = preprocessed_bkb

            if save_dir is not None:
                query.save(save_dir)

            if interpolation == 'independence':
                return self.solve_query_independence(query, genetic_evidence, target_strategy)
            return self.solve_query(query)

        #-- Else use the collapsed bkb that was processed during setup.
        bkb = copy.deepcopy(self.collapsed_bkb)

        if query.meta_targets is not None:
            transformed_meta_targets = ['{} {} {}'.format(target[0], target[1], target[2]) for target in query.meta_targets]

        transformed_meta, bkb = self.process_metaVariables(bkb, query.meta_evidence, query.meta_targets, target_strategy)

        #-- Collect Evidence
        transformed_meta_evidence = dict()
        if query.meta_evidence is not None:
            for ev in query.meta_evidence:
                meta_name = '{} {} {}'.format(ev[0], ev[1], ev[2])
                transformed_meta_evidence[meta_name] = transformed_meta[meta_name]

        #-- Collect Targets differently if using topological target strategy
        if target_strategy == 'topological':
            transformed_meta_targets = dict()
            for target in query.meta_targets:
                meta_name = '{} {} {}'.format(target[0], target[1], target[2])
                for name in transformed_meta:
                    if meta_name in name:
                        transformed_meta_targets[name] = transformed_meta[name]

        query.evidence.update(transformed_meta_evidence)
        query.targets.extend(transformed_meta_targets)

        query.bkb = bkb

        if save_dir is not None:
            query.save(save_dir)
        if interpolation == 'independence':
            return self.solve_query_independence(query, genetic_evidence, target_strategy)
        return self.solve_query(query, target_strategy=target_strategy)

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
                #-- Check if the src_population item is a list of items (useful for drugs):
                if type(src_population_data[src_name][prop_]) == tuple:
                    if val_ in src_population_data[src_name][prop_]:
                        res = True
                    else:
                        res = False
                else:
                    res = op_(src_population_data[src_name][prop_], val_)
                truth.append(res == state)
            if all(truth):
                count += 1
        counts.append(count)

    probs = [float(count) / len(src_population) for count in counts]

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

def _collectBkbBottomINodes(bkb):
    #-- Collect all I-nodes that are not present in any snode tails. Probably can optimize this. Needs a cycle check.
    #print(bkb)
    heads = set()
    tails = set()
    for snode in bkb.getAllSNodes():
        heads.add(snode.getHead())
        for tail_idx in range(snode.getNumberTail()):
            tails.add(snode.getTail(tail_idx))
    #-- This picks out any inode that appears in an snode head put not in any snode tails, i.e. the bottom node.
    bottom_inodes = heads - tails
    return bottom_inodes

def _addTargetToLastTopologVariables(target, bkb, src_population, src_population_data):
    bottom_inodes = _collectBkbBottomINodes(bkb)

    prop, op_str, val = target
    op = _process_operator(op_str)

    transformed_meta = dict()
    for comp_idx, state_idx in bottom_inodes:
        comp_name = bkb.getComponentName(comp_idx)
        if 'mut-var_' == comp_name[:8]:
            #-- If this is a variant component
            gene = '_'.join(comp_name.split('_')[1:])
            variant = bkb.getComponentINodeName(comp_idx, state_idx)
            gene_variant = '{}-{}'.format(gene, variant)
            joint_count = 0
            prior_count = 0
            neg_joint_count = 0
            for src_hash, data_dict in src_population_data.items():
                try:
                    if gene_variant in data_dict['Patient_Gene_Variants']:
                        prior_count += 1
                        if op(src_population_data[src_hash][prop], val):
                            joint_count += 1
                        else:
                            neg_joint_count += 1
                except:
                    continue

            prob_true_cond = joint_count / prior_count
            prob_false_cond = neg_joint_count / prior_count

            #-- Add target Inode and attach s-node
            target_comp_idx = bkb.addComponent('{} {} {}_mut-var_{}_{}'.format(prop, op_str, val,
                                                                                gene, variant))
            if prob_true_cond > 0:
                target_true_idx = bkb.addComponentState(target_comp_idx, 'True')
                bkb.addSNode(BKB_S_node(init_component_index=target_comp_idx,
                                        init_state_index=target_true_idx,
                                        init_probability=prob_true_cond,
                                        init_tail=[(comp_idx, state_idx)]))
            if prob_false_cond > 0:
                target_false_idx = bkb.addComponentState(target_comp_idx, 'False')
                bkb.addSNode(BKB_S_node(init_component_index=target_comp_idx,
                                        init_state_index=target_false_idx,
                                        init_probability=prob_false_cond,
                                        init_tail=[(comp_idx, state_idx)]))
        elif 'mut_' == comp_name[:4]:
            #-- If this is a mutation component
            gene = '_'.join(comp_name.split('_')[1:])
            joint_count = 0
            prior_count = 0
            neg_joint_count = 0
            for src_hash, data_dict in src_population_data.items():
                try:
                    if gene in data_dict['Patient_Genes']:
                        prior_count += 1
                        if op(src_population_data[src_hash][prop], val):
                            joint_count += 1
                        else:
                            neg_joint_count += 1
                except:
                    continue
            prob_true_cond = joint_count / prior_count
            prob_false_cond = neg_joint_count / prior_count

            #-- Add target Inode and attach s-node
            target_comp_idx = bkb.addComponent('{} {} {}_mut_{}'.format(prop, op_str, val, gene))
            if prob_true_cond > 0:
                target_true_idx = bkb.addComponentState(target_comp_idx, 'True')
                bkb.addSNode(BKB_S_node(init_component_index=target_comp_idx,
                                        init_state_index=target_true_idx,
                                        init_probability=prob_true_cond,
                                        init_tail=[(comp_idx, state_idx)]))
            if prob_false_cond > 0:
                target_false_idx = bkb.addComponentState(target_comp_idx, 'False')
                bkb.addSNode(BKB_S_node(init_component_index=target_comp_idx,
                                        init_state_index=target_false_idx,
                                        init_probability=prob_false_cond,
                                        init_tail=[(comp_idx, state_idx)]))
        transformed_meta[bkb.getComponentName(target_comp_idx)] = 'True'
    return bkb, transformed_meta

def _addDemographicOption(option, bkb, src_population, src_population_data, option_dependencies=list()):
    prop, op_str, val = option
    op = _process_operator(op_str)
    matched_srcs = set()
    pop_count_true = 0
    for entity_name, src_name in src_population.items():
        #print(src_population_data[src_name][prop])
        if type(src_population_data[src_name][prop]) == tuple:
            if val in src_population_data[src_name][prop]: #op(set(val).issubset(set(src_population_data[src_name][prop]))):
                matched_srcs.add(entity_name)
                pop_count_true += 1
                res = True
            else:
                res = False
        else:
            if op(src_population_data[src_name][prop], val):
                matched_srcs.add(entity_name)
                pop_count_true += 1
    prob = float(pop_count_true / len(src_population))
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

def _linkSource(src_comp, src_state, non_src_comp, non_src_state, other_non_src_tails, prob,
                comp_idx, inode_true_idx, inode_false_idx, bkb, matched_srcs, src_population):
    #-- Process the source name collection
    src_state_name = bkb.getComponentINodeName(src_comp, src_state)
    src_split = src_state_name.split('_')
    src_nums_str = src_split[0][1:-1]
    src_name_str = ''.join(src_split[1:])
    src_nums = [int(num) for num in src_nums_str.split(',')]
    src_names = src_name_str.split(',')

    #-- Delete the source node and snodes
    #-- Actually leaving states in for speed.
    #bkb.removeComponentState(src_comp, src_state)

    true_srcs = list()
    false_srcs = list()
    for src_num, src_name in zip(src_nums, src_names):
        if src_num in matched_srcs:
            true_srcs.append((src_num, src_name))
        else:
            false_srcs.append((src_num, src_name))

    true_prob = len(true_srcs) / (len(true_srcs) + len(false_srcs))

    #-- Make new true and false collections
    true_state_name = '[{}]_{}'.format(','.join([str(src[0]) for src in true_srcs]), ','.join([src[1] for src in true_srcs]))
    false_state_name = '[{}]_{}'.format(','.join([str(src[0]) for src in false_srcs]), ','.join([src[1] for src in false_srcs]))

    #-- Now add the I nodes
    if true_prob > 0:
        src_true_idx = bkb.addComponentState(src_comp, true_state_name)
        #-- Now connect with S nodes
        bkb.addSNode(BKB_S_node(init_component_index=non_src_comp,
                                init_state_index=non_src_state,
                                init_probability=1,
                                init_tail=[(src_comp, src_true_idx)] + other_non_src_tails))
        bkb.addSNode(BKB_S_node(init_component_index=src_comp,
                                init_state_index=src_true_idx,
                                init_probability=prob*true_prob,
                                init_tail=[(comp_idx, inode_true_idx)]))
    if (1 - true_prob) > 0:
        src_false_idx = bkb.addComponentState(src_comp, false_state_name)
        #-- Now connect with S nodes
        bkb.addSNode(BKB_S_node(init_component_index=non_src_comp,
                                init_state_index=non_src_state,
                                init_probability=1,
                                init_tail=[(src_comp, src_false_idx)] + other_non_src_tails))
        bkb.addSNode(BKB_S_node(init_component_index=src_comp,
                                init_state_index=src_false_idx,
                                init_probability=prob*(1 - true_prob),
                                init_tail=[(comp_idx, inode_false_idx)]))

    return bkb

def _addSrcConnections(comp_idx, inode_true_idx, inode_false_idx, bkb, matched_srcs, src_population):
    #-- Attach to source nodes
    src_comp_indices = bkb.getSrcComponents()
    non_src_comp_indices = set(bkb.getAllComponentIndices()) - set(src_comp_indices)
    S_nodes_by_head = bkb.constructSNodesByHead()

    for non_src_comp in tqdm.tqdm(non_src_comp_indices, desc='Linking Sources', leave=False):
        for non_src_state in bkb.getAllComponentINodeIndices(non_src_comp):
            snodes = S_nodes_by_head[non_src_comp][non_src_state]
            for snode in snodes:
                src_tails = list()
                non_src_tails = list()
                for tail_idx in range(snode.getNumberTail()):
                    tail_comp, tail_state = snode.getTail(tail_idx)
                    #-- If this is a source component
                    if tail_comp in src_comp_indices:
                        src_tails.append((tail_comp, tail_state))
                    else:
                        non_src_tails.append((tail_comp, tail_state))
                for src_tail_comp, src_tail_state in src_tails:
                    prior_src_snode = S_nodes_by_head[src_tail_comp][src_tail_state][0]
                    prob = prior_src_snode.probability
                    bkb = _linkSource(src_tail_comp, src_tail_state, non_src_comp, non_src_state, non_src_tails, prob,
                                      comp_idx, inode_true_idx, inode_false_idx, bkb, matched_srcs, src_population)
                    #-- delete the prior snode
                    try:
                        bkb.removeSNode(prior_src_snode)
                    except:
                        pass
                if len(src_tails) > 0:
                    #-- Delete other S node
                    try:
                        bkb.removeSNode(snode)
                    except:
                        pass
    return bkb
