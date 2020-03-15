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
import zlib

sys.path.append('/home/cyakaboski/src/python/projects/bkb-pathway-provider/core')
#sys.path.append('/home/ghyde/bkb-pathway-provider/core')
from query import Query

from pybkb.core.cpp_base.reasoning import revision as cpp_revision
from pybkb.core.cpp_base.reasoning import updating as cpp_updating
from pybkb.core.python_base.reasoning import updating as py_updating
from pybkb.core.python_base.reasoning import checkMutex
from pybkb.core.python_base.fusion import fuse
from pybkb.core.python_base.fusion_collapse import collapse_sources
from pybkb.core.common.bayesianKnowledgeBase import BKB_S_node
from pybkb.core.common.bayesianKnowledgeBase import bayesianKnowledgeBase as BKB

CACHED_BKB_DIR = '/home/public/data/ncats/cachedCollapsedBkb'
ILLEGAL_SOURCE_STRINGS = ['PatientX', 'noGeneEvidence', 'geneEvidence']

class Reasoner:
    def __init__(self, fused_bkb=None, collapsed_bkb=None,  patient_data=None, cpp_reasoning=False, gene_var_direct=None, max_new_ev=None):
        self.fused_bkb = fused_bkb
        self.metadata = patient_data
        self.cpp_reasoning = cpp_reasoning
        self.metadata = patient_data
        self.collapsed_bkb = collapsed_bkb
        if patient_data is not None:
            try:
                self.setup()
            except:
                raise ValueError('Likely patient data is not in the right format.')
        self.gene_var_direct = gene_var_direct
        self.max_new_ev = max_new_ev

    def getCollapsedBKB(self, fused_bkb):
        #-- First check to see if this bkb has already been collapsed and saved.
        collapsed_bkb_hash_name = zlib.adler32(fused_bkb.to_str().encode('utf-8'))
        collapsed_bkb_path = os.path.join(CACHED_BKB_DIR, '{}.bkb'.format(collapsed_bkb_hash_name))
        if os.path.exists(collapsed_bkb_path):
            print('Loaded collapsed BKB from memory.')
            collapsed_bkb = BKB()
            collapsed_bkb.load(collapsed_bkb_path)
        else:
            collapsed_bkb = collapse_sources(fused_bkb)
            #-- Save collapsed bkb using fused_bkb hash as file name.
            collapsed_bkb.save(collapsed_bkb_path)
        return collapsed_bkb

    def setup(self):
        #-- Collapse sources in fused bkb
        if self.collapsed_bkb is None:
            if self.fused_bkb is None:
                raise ValueError('You must pass either a fused or already collapsed bkb.')
            self.collapsed_bkb = self.getCollapsedBKB(self.fused_bkb)
        #self.collapsed_bkb.save('collapsed.bkb')

        #-- Preprocess src hash values
        src_component_indices = self.fused_bkb.getSrcComponents()
        src_hashs = dict()
        src_hashs_inverse = dict()
        for comp_idx in src_component_indices:
            for state_idx in range(self.fused_bkb.getNumberComponentINodes(comp_idx)):
                state_name = self.fused_bkb.getComponentINodeName(comp_idx, state_idx)
                #-- Collect src numbers and string name or hash
                src_num = int(state_name.split('_')[0][1:-1])
                try:
                    src_hash = int(''.join(state_name.split('_')[1:]))
                    src_hashs[src_num] = src_hash
                    src_hashs_inverse[src_hash] = src_num
                except ValueError:
                    continue
                    #src_hashs[src_num] = state_name
        #self.src_hashs = set(src_hashs)
        self.src_hashs = src_hashs
        self.src_hashs_inverse = src_hashs_inverse

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
        #print(checkMutex(query.bkb))
        #query.bkb.makeGraph()
        #print(query.evidence)
        #input()
        #print(query.targets)
        #input()
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
            #print(updates)
            #input()
            #-- Convert to actual meta target probability using independence assumption.
            meta_target_updates = dict()
            for comp_name, state_dict in updates.items():
                #-- There should always be two states, i.e. True and False. If theres only one just omit it.
                if min(state_dict.values()) < 0:
                    continue
                if len(state_dict) < 2:
                    continue
                meta_target_name = comp_name.split('_')[0]
                if meta_target_name not in meta_target_updates:
                    meta_target_updates[meta_target_name] = state_dict
                elif min(state_dict.values()) > 0:
                    for state_name, prob in state_dict.items():
                        if meta_target_updates[meta_target_name][state_name] < 0:
                            meta_target_updates[meta_target_name][state_name] = prob
                        else:
                            meta_target_updates[meta_target_name][state_name] *= prob
                #-- Normalize
                total_prob = 0
                for state_name, unnorm_prob in meta_target_updates[meta_target_name].items():
                    total_prob += unnorm_prob
                new_state_dict = dict()
                for state_name, unnorm_prob in meta_target_updates[meta_target_name].items():
                    new_state_dict[state_name] = unnorm_prob / total_prob
                meta_target_updates[meta_target_name] = new_state_dict
            if len(meta_target_updates) > 0:
                query.result.meta_target_updates = meta_target_updates
        return query

    def process_metaVariables(self, bkb, meta_evidence, meta_targets, target_strategy='chain', num_gene_evidence=0):
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
        elif target_strategy == 'explicit':
            meta_variables = []
            if meta_evidence is not None:
                meta_variables += meta_evidence
            if meta_targets is not None:
                meta_variables += meta_targets

            #-- First make a chained demographic bkb
            demo_chain_bkb = BKB()
            for i, meta in enumerate(meta_variables):
                demo_chained_bkb, _, _ = _addDemographicOption(meta, demo_chain_bkb, self.src_hashs, self.metadata, option_dependencies=meta_variables[:i], include_src_tags=True)
            #print('Fuck YOu')
            #demo_chained_bkb.makeGraph()
            #print(checkMutex(demo_chained_bkb))
            demo_dummy_bkb = _makeDemoDummyBKB(meta_variables)
            #demo_dummy_bkb.makeGraph()
            #-- Fuse demographically chained bkb with dummy BKB
            demo_fused_bkb = fuse([demo_chained_bkb, demo_dummy_bkb],
                                  [1, 1],
                                  ['noGeneEvidence', 'geneEvidence'])
            demo_fused_bkb.save('demo_fused_t1.bkb')
            #-- Clean demo fusion graph, i.e. delete weird inodes
            S_nodes_by_head = demo_fused_bkb.constructSNodesByHead()
            src_comp_indices = demo_fused_bkb.getSrcComponents()
            for src_comp in src_comp_indices:
                for src_state in demo_fused_bkb.getAllComponentINodeIndices(src_comp):
                    if len(S_nodes_by_head[src_comp][src_state]) > 0:
                        #print('Here')
                        for snode in S_nodes_by_head[src_comp][src_state][1:]:
                            demo_fused_bkb.removeSNode(snode)
            #print(checkMutex(demo_fused_bkb))
            #demo_fused_bkb.makeGraph()
            #input()

            #-- Construct geneEvidence helper
            helper_evidence = _getExplicitHelperEvidence(demo_fused_bkb, num_gene_evidence)
            '''
            if num_gene_evidence == 0:
                src_comp_indices = demo_fused_bkb.getSrcComponents()
                for src_comp in src_comp_indices:
                        src_state = demo_fused_bkb.findINode(src_comp, 'noGeneEvidence', contains=True)
                        if src_state != -1:
                            helper_evidence[demo_fused_bkb.getComponentName(src_comp)] = demo_fused_bkb.getComponentINodeName(src_comp, src_state)
            else:
                src_comp_indices = demo_fused_bkb.getSrcComponents()
                for src_comp in src_comp_indices:
                        src_state = demo_fused_bkb.findINode(src_comp, 'geneEvidence', contains=True)
                        if src_state != -1:
                            helper_evidence[demo_fused_bkb.getComponentName(src_comp)] = demo_fused_bkb.getComponentINodeName(src_comp, src_state)
            print(helper_evidence)
            '''
            #-- Merge demographic fused bkb into gene bkb
            for comp_name, state_name in demo_fused_bkb.getINodeNames():
                comp_idx = bkb.addComponent(comp_name)
                bkb.addComponentState(comp_idx, state_name)
            for snode in demo_fused_bkb.getAllSNodes():
                head_comp, head_state = snode.getHead()
                head_comp_name = demo_fused_bkb.getComponentName(head_comp)
                head_state_name = demo_fused_bkb.getComponentINodeName(head_comp, head_state)
                head_comp_idx = bkb.getComponentIndex(head_comp_name)
                head_state_idx = bkb.getComponentINodeIndex(head_comp_idx, head_state_name)
                tail = list()
                for tail_idx in range(snode.getNumberTail()):
                    tail_comp, tail_state = snode.getTail(tail_idx)
                    tail_comp_name = demo_fused_bkb.getComponentName(tail_comp)
                    tail_state_name = demo_fused_bkb.getComponentINodeName(tail_comp, tail_state)
                    tail_comp_idx = bkb.getComponentIndex(tail_comp_name)
                    tail_state_idx = bkb.getComponentINodeIndex(tail_comp_idx, tail_state_name)
                    tail.append((tail_comp_idx, tail_state_idx))
                prob = snode.probability
                bkb.addSNode(BKB_S_node(init_component_index=head_comp_idx,
                                        init_state_index=head_state_idx,
                                        init_probability=prob,
                                        init_tail=tail))
            #bkb.makeGraph()



            #-- Process demo options explicity
            bkb = _addDemographicOptionsExplicitly(meta_variables, bkb, self.src_hashs, self.src_hashs_inverse, self.metadata)
            #-- Construct transformed meta variables
            transformed_meta = dict()
            for meta in meta_variables:
                transformed_meta['{} {} {}'.format(meta[0], meta[1], meta[2])] = 'True'

            #-- Update with gene helper
            transformed_meta.update(helper_evidence)
            print(transformed_meta)
            return transformed_meta, bkb
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
            bkb = _addSrcConnections(comp_idx, inode_true_idx, inode_false_idx, bkb, matched_srcs, self.src_hashs, self.src_hashs_inverse)

        #-- If topological target stradegy connect target to each bottom I node.
        if target_strategy == 'topological':
            #-- Construct link and stats from bottom I-nodes to individual meta targets.
            bkb, transformed_meta_ = _addTargetToLastTopologVariables(meta_targets[0], bkb, self.src_hashs, self.metadata)
            transformed_meta.update(transformed_meta_)
        return transformed_meta, bkb

    def solve_query_independence(self, query, genetic_evidence, target_strategy, parallel=False):
        if len(genetic_evidence) == 0:
            return self.solve_query(query, target_strategy)
        bkb = query.bkb
        non_gene_evidence = copy.deepcopy(query.evidence)
        for gene_key in genetic_evidence:
            del non_gene_evidence[gene_key]

        def make_independent_query(base_query, bkb, non_gene_evidence, genetic_comp, genetic_state, parallel):
            independ_evidence = copy.deepcopy(non_gene_evidence)
            independ_evidence[genetic_comp] = genetic_state
            #print(independ_evidence)
            #input('Stopped')
            q = Query(evidence=independ_evidence,
                      targets=copy.deepcopy(base_query.targets),
                      type='updating')

            #if not parallel:
            #    q.bkb = copy.deepcopy(bkb)
            #else:
            q.bkb = bkb
            q.patient_data = self.metadata
            return q

        queries = list()
        for genetic_comp, genetic_state in tqdm.tqdm(genetic_evidence.items(), desc='Setting up independent runs', leave=False):
            queries.append(make_independent_query(query, bkb, non_gene_evidence, genetic_comp, genetic_state, parallel))

        if not parallel:
            finished_queries = list()
            for q_ in tqdm.tqdm(queries, desc='Solving Independent Queries', leave=False):
                finished_queries.append(self.solve_query(q_, target_strategy=target_strategy))

        else:
            start_time = time.time()
            #-- Run all independent queries in parrellel
            with ProcessPoolExecutor(max_workers=30) as executor:
                #finished_queries = list()
                #for q_ in executor.map(self.solve_query, [(q__, target_strategy) for q__ in queries]):
                #    finished_queries.append(q_)
                futures = list()
                for q_ in queries:
                    futures.append(executor.submit(self.solve_query, q_, copy.deepcopy(target_strategy)))
                finished_queries = [q_.result() for q_ in futures]
            print('Multiprocessing time: {}'.format(time.time() - start_time))

        #-- Collect queries and calculate independent probs
        res = dict()
        for q_ in finished_queries:
            updates = q_.result.process_updates()
            #print(updates)
            #input('Stopped')
            for comp_name, state_dict in updates.items():
                state_keys = state_dict.keys()
                state_probs = state_dict.values()
                # We were skipping -1 where there is no inference on a particular state_key but now lets just make it really small
                if min(state_probs) <= 0:
                    state_probs_ = list()
                    for key, prob_ in zip(state_keys, state_probs):
                        if prob_ < 0:
                            state_probs_.append(1e-10)
                            state_dict[key] = 1e-10
                        else:
                            state_probs_.append(prob_)
                            state_dict[key] = prob_
                    state_probs = state_probs_
                if len(state_dict) < 2:
                    continue
                if comp_name in res:
                    sumStateProbs = sum(state_probs)
                    for state_name, prob in state_dict.items():
                        ##-- If probability is greater than one than its a log prob
                        #if prob > 1:
                        #    res[comp_name][state_name] += prob
                        #else:
                        res[comp_name][state_name] *= (prob/sumStateProbs)
                else:
                    res[comp_name] = state_dict
                resProbs = res[comp_name].values()
                sumResProbs = sum(resProbs)
                for state_key in res[comp_name].keys():
                    res[comp_name][state_key] /= sumResProbs
        query.independ_queries = queries
        query.independ_result = res
        #for q in query.independ_queries:
        #    q.getReport()
        #    input('Stopped')
        return query

    '''
    Target strategy can be chain, where targets are located in the chain rule (legacy) or
    topological where target is attached to all nodes at the bottom of the topological sort. This
    methodolgy can only except on target then.
    '''
    def analyze_query(self, query, save_dir=None, preprocessed_bkb=None, target_strategy='chain', interpolation='independence', check_mutex=False):
        #-- Set up query parameters
        query.patient_data = self.metadata
        query.target_strategy = target_strategy
        query.interpolation = interpolation

        if self.gene_var_direct is not None and self.max_new_ev is not None:
                    query.gene_var_direct = self.gene_var_direct
                    query.max_new_ev = self.max_new_ev
                    needed = query.checkAndAdjustEvidence()

        #-- Check if there is any genetic evidence
        num_gene_evidence = len(query.evidence)

        #-- Make a bkb query hash.
        query_bkb_hash = zlib.adler32(''.join([self.collapsed_bkb.to_str(),
                                               str(query.meta_evidence),
                                               str(query.meta_targets),
                                               target_strategy,
                                               interpolation]).encode('utf-8'))
        query_bkb_path = os.path.join(CACHED_BKB_DIR, '{}.bkb'.format(query_bkb_hash))

        #-- Duplicate Gene Evidence
        genetic_evidence = copy.deepcopy(query.evidence)
        #-- See if we can find a preprocessed bkb.
        if preprocessed_bkb is not None:
            query.bkb = preprocessed_bkb
        elif os.path.exists(query_bkb_path):
            print('Loaded query BKB from memory.')
            query.bkb = BKB()
            query.bkb.load(query_bkb_path)

        #-- If we have a preprocessed passed or from memory BKB
        if query.bkb is not None:
            if query.meta_targets is not None:
                transformed_meta_targets = ['{} {} {}'.format(target[0], target[1], target[2]) for target in query.meta_targets]

            if query.meta_evidence is not None:
                transformed_meta_evidence = dict()
                for ev in query.meta_evidence:
                    meta_name = '{} {} {}'.format(ev[0], ev[1], ev[2])
                    transformed_meta_evidence[meta_name] = 'True'
                    query.evidence.update(transformed_meta_evidence)

            #--Collect all helper evidence
            if target_strategy == 'explicit':
                helper_evidence = _getExplicitHelperEvidence(query.bkb, num_gene_evidence)
                query.evidence.update(helper_evidence)

            query.targets.extend(transformed_meta_targets)

            if check_mutex:
                print('No Mutex Issues:', checkMutex(query.bkb))
            if save_dir is not None:
                query.save(save_dir)

            if interpolation == 'independence':
                return self.solve_query_independence(query, genetic_evidence, target_strategy)
            return self.solve_query(query)

        #-- Else use the collapsed bkb that was processed during setup.
        bkb = copy.deepcopy(self.collapsed_bkb)

        if query.meta_targets is not None:
            transformed_meta_targets = ['{} {} {}'.format(target[0], target[1], target[2]) for target in query.meta_targets]

        transformed_meta, bkb = self.process_metaVariables(bkb, query.meta_evidence, query.meta_targets, target_strategy, num_gene_evidence=num_gene_evidence)

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
        #--Collect all helper evidence
        if target_strategy == 'explicit':
            helper_evidence = dict()
            for comp, state in transformed_meta.items():
                if comp not in transformed_meta_evidence and comp not in transformed_meta_targets:
                    helper_evidence[comp] = state
            query.evidence.update(helper_evidence)

        query.evidence.update(transformed_meta_evidence)
        query.targets.extend(transformed_meta_targets)

        query.bkb = bkb
        if check_mutex:
            print('No Mutex Issues:', checkMutex(bkb))
        bkb.save(query_bkb_path)
        #bkb.save('collapsed_and_link.bkb')
        #input('Saved')
        if save_dir is not None:
            query.save(save_dir)
        if interpolation == 'independence':
            return self.solve_query_independence(query, genetic_evidence, target_strategy)
        return self.solve_query(query, target_strategy=target_strategy)

def _getExplicitHelperEvidence(bkb, num_gene_evidence):
    #-- Construct geneEvidence helper
    helper_evidence = dict()
    if num_gene_evidence == 0:
        src_comp_indices = bkb.getSrcComponents()
        for src_comp in src_comp_indices:
                src_state = bkb.findINode(src_comp, 'noGeneEvidence', contains=True)
                if src_state != -1:
                    helper_evidence[bkb.getComponentName(src_comp)] = bkb.getComponentINodeName(src_comp, src_state)
    else:
        src_comp_indices = bkb.getSrcComponents()
        for src_comp in src_comp_indices:
                src_state = bkb.findINode(src_comp, 'geneEvidence', contains=True)
                if src_state != -1:
                    helper_evidence[bkb.getComponentName(src_comp)] = bkb.getComponentINodeName(src_comp, src_state)
    return helper_evidence


def _makeDemoDummyBKB(meta_variables):
    bkb = BKB()
    for meta in meta_variables:
        comp_idx = bkb.addComponent('{} {} {}'.format(meta[0], meta[1], meta[2]))
        stateTrue_idx = bkb.addComponentState(comp_idx, 'True')
        stateFalse_idx = bkb.addComponentState(comp_idx, 'False')
        #-- Add snodes
        bkb.addSNode(BKB_S_node(init_component_index=comp_idx,
                                init_state_index=stateTrue_idx, 
                                init_probability=1))
        bkb.addSNode(BKB_S_node(init_component_index=comp_idx,
                                init_state_index=stateFalse_idx, 
                                init_probability=1))
    return bkb

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

def _getDemographicCombos(meta_variables):
    #-- Get demographic options
    demo_product = [combo for combo in itertools.product(meta_variables, [True, False])]
    #-- Get only legal demographic combos, i.e. we don't want a combo like [(Age >= 50 == True), (Age >= 50 = False), (Surival >= 365 = True)]...
    demo_combos_list = list()
    for combo in itertools.combinations(demo_product, r=len(meta_variables)):
        combo = list(combo)
        prop_set = set()
        for ev_state in combo:
            ev_, state = ev_state
            prop_, op_, val_ = ev_
            prop_set.add(prop_)
        if len(prop_set) == len(combo):
            demo_combos_list.append(tuple(combo))

    #-- Instiante a demographic combination dict to capture which sources match each demo combo.
    demo_combos = {combo: list() for combo in demo_combos_list}
    return demo_combos

def _linkSrcToDemographicCombinations(demo_combos, bkb, non_src_head, src, non_src_tail, prior_prob, src_hashs, src_hashs_inverse, src_population_data, processed_meta_variable_priors):
    #-- Copy democombos dictionary
    demo_combos = copy.deepcopy(demo_combos)
    
    src_comp, src_state = src
    #-- Process the source name collection
    src_state_name = bkb.getComponentINodeName(src_comp, src_state)
    src_split = src_state_name.split('_')
    src_nums_str = src_split[0][1:-1]
    src_name_str = ''.join(src_split[1:])
    src_names_str = src_name_str.split(',')
    src_names = [int(src_hash) for src_hash in src_names_str]
    src_nums = [int(num) for num in src_nums_str.split(',')]
    #-- If this source component has format [num]_hash1,hash2, ... then it used fusion hack.
    if len(src_nums) != len(src_names):
        #-- Go through and process src nums by looking at hashes.
        src_nums = [src_hashs_inverse[int(src_hash)] for src_hash in src_names_str]

    #-- Collect all sources that match the respective meta variable combinations
    for src_name in src_names:
        for combo in demo_combos:
            truth = list()
            for meta_var in combo:
                var_, state = meta_var
                prop_, op_str_, val_ = var_
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
                demo_combos[combo].append(src_name)
    #print(demo_combos)
    #input()
    #-- Build all the joint snodes
    for combo, srcs in demo_combos.items():
        if len(srcs) == 0:
            continue
        #-- Remake the source collection component
        src_nums_ = [src_hashs_inverse[src_hash] for src_hash in srcs]
        src_state_name = '[{}]_{}'.format(','.join([str(src_num) for src_num in src_nums]),
                                          ','.join([str(src_hash) for src_hash in srcs]))
        new_src_state = bkb.addComponentState(src_comp, src_state_name)
        new_state_prob = float(len(srcs) / len(src_names))
        #-- Add new source collection prior
        bkb.addSNode(BKB_S_node(init_component_index=src_comp,
                                init_state_index=new_src_state,
                                init_probability=new_state_prob*prior_prob))
        #-- Get the demographic combination tail
        combo_tail = list()
        for meta_var in combo:
            var_, state = meta_var
            var_comp = bkb.addComponent('{} {} {}'.format(var_[0], var_[1], var_[2]))
            if state:
                var_state = bkb.addComponentState(var_comp, 'True')
            else:
                var_state = bkb.addComponentState(var_comp, 'False')
            #-- Link src to each demographic option
            #bkb.addSNode(BKB_S_node(init_component_index=var_comp,
            #                        init_state_index=var_state,
            #                        init_probability=1,
            #                        init_tail=[src]))
            #processed_meta_variable_priors.add((var_comp, var_state))
            combo_tail.append((var_comp, var_state))
        #-- Make Joint S-node
        complete_tail = combo_tail + [(src_comp, new_src_state)] + non_src_tail
        non_src_head_comp, non_src_head_state = non_src_head
        bkb.addSNode(BKB_S_node(init_component_index=non_src_head_comp,
                                init_state_index=non_src_head_state,
                                init_probability=1,
                                init_tail=complete_tail))
    return bkb, processed_meta_variable_priors

def _addDemographicOptionsExplicitly(meta_variables, bkb, src_hashs, src_hashs_inverse, src_population_data):
    src_comp_indices = bkb.getSrcComponents()
    non_src_comp_indices = set(bkb.getAllComponentIndices()) - set(src_comp_indices)
    S_nodes_by_head = bkb.constructSNodesByHead()

    illegal_source_strings = ILLEGAL_SOURCE_STRINGS + [meta[0] for meta in meta_variables]

    processed_meta_variable_priors = set()
    #-- Construct demographic combo dictionary
    demo_combos = _getDemographicCombos(meta_variables)

    for non_src_comp in tqdm.tqdm(non_src_comp_indices, desc='Linking Sources', leave=False):
        for non_src_state in bkb.getAllComponentINodeIndices(non_src_comp):
            #print(bkb.getComponentName(non_src_comp), bkb.getComponentINodeName(*(non_src_comp, non_src_state)))
            try:
                snodes = S_nodes_by_head[non_src_comp][non_src_state]
            except KeyError:
                print('Warning: {} = {} has no incoming S nodes.'.format(bkb.getComponentName(non_src_comp), bkb.getComponentINodeName(non_src_comp, non_src_state)))
                continue
            for snode in snodes:
                non_src_tail = list()
                for tail_idx in range(snode.getNumberTail()):
                    tail_comp, tail_state = snode.getTail(tail_idx)
                    if tail_comp in src_comp_indices:
                        src_ = (tail_comp, tail_state)
                    else:
                        non_src_tail.append((tail_comp, tail_state))
                #-- Check for illegal source before linking. 
                found_illegal = False
                for illegal_src_val in illegal_source_strings:
                    if illegal_src_val in bkb.getComponentINodeName(*src_):
                        found_illegal = True
                        break
                if found_illegal:
                    continue
                prior_snode = S_nodes_by_head[src_[0]][src_[1]][0]
                prior_prob = prior_snode.probability
                bkb, processed_meta_variable_priors = _linkSrcToDemographicCombinations(demo_combos,
                                                                                        bkb,
                                                                                        (non_src_comp,
                                                                                         non_src_state),
                                                                                        src_,
                                                                                        non_src_tail,
                                                                                        prior_prob,
                                                                                        src_hashs,
                                                                                        src_hashs_inverse,
                                                                                        src_population_data,
                                                                                        processed_meta_variable_priors)
                #-- delete the prior snode
                try:
                    bkb.removeSNode(prior_snode)
                except:
                    pass
                if not found_illegal:
                    #-- Delete other S node
                    try:
                        bkb.removeSNode(snode)
                    except:
                        pass
    return bkb

def _processOptionDependency(option, option_dependencies, bkb, matched_srcs, src_population, src_population_data, include_src_tags):
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
        count_joint = 0
        count_prior = 0
        for entity_name, src_name in src_population.items():
            truth_joint = list()
            truth_prior = list()
            for k, ev_state in enumerate(combo):
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
                truth_joint.append(res == state)
                if k > 0:
                    truth_prior.append(res == state)
            if all(truth_joint):
                count_joint += 1
            if all(truth_prior):
                count_prior += 1
        counts.append((count_joint, count_prior))
    probs_joint = [float(count[0]) / len(src_population) for count in counts]
    probs_prior = [float(count[1]) / len(src_population) for count in counts]
    probs_cond = []
    for idx, prob_join in enumerate(probs_joint):
        if probs_prior[idx] == 0:
            probs_cond.append(0)
        else:
            probs_cond.append( prob_join / probs_prior[idx] )
    #probs_cond = [probs_joint[i] / probs_prior[i] for i in range(len(counts))]

    processed_src_tags = set()
    #-- Setup each S-node
    for j, combo in enumerate(combos):
        head, tail = combo
        #-- Process head
        comp_head_idx, i_node_head_idx = _processDependencyHead(head, bkb)
        #-- Process Tail
        processed_tail = _processDependecyTail(tail, bkb)
        if include_src_tags:
            bkb, src_tag, processed_src_tags = _processSrcTags(bkb, comp_head_idx, i_node_head_idx, matched_srcs, src_population, processed_src_tags)
            src_tag = [src_tag]
        else:
            src_tag = list()
        #-- Add Snode
        if probs_cond[j] > 0:
            bkb.addSNode(BKB_S_node(init_component_index=comp_head_idx, init_state_index=i_node_head_idx, init_probability=probs_cond[j], init_tail=processed_tail + src_tag))
    return bkb

def _processSrcTags(bkb, comp_head, state_head, matched_srcs, src_population, processed_src_tags):
    #print(matched_srcs)
    #print(src_population)
    #input()
    #print('Processing: {} {}'.format(comp_head, state_head))
    #input()
    if bkb.getComponentINodeName(comp_head, state_head) == 'True':
        src_name = '[{}]_{}'.format(','.join([str(src_num) for src_num in matched_srcs]),
                                    ','.join([str(src_population[src_num]) for src_num in matched_srcs]))
    else:
        src_name = '[{}]_{}'.format(','.join([str(src_num) for src_num in set(src_population.keys()) - matched_srcs]),
                                     ','.join([str(src_population[src_num]) for src_num in set(src_population.keys()) - matched_srcs]))

    src_tag_comp = bkb.addComponent('Collection[{} = {}]'.format(bkb.getComponentName(comp_head), bkb.getComponentINodeName(comp_head, state_head)))
    src_tag_state = bkb.addComponentState(src_tag_comp, src_name)
    #-- Add tag priors
    if (src_tag_comp, src_tag_state) not in processed_src_tags:
        bkb.addSNode(BKB_S_node(init_component_index=src_tag_comp,
                                init_state_index=src_tag_state,
                                init_probability=1))
        processed_src_tags.add((src_tag_comp, src_tag_state))
    return bkb, (src_tag_comp, src_tag_state), processed_src_tags
'''
    if other_tail is None:
        snode_1 = BKB_S_node(init_component_index=comp_idx, init_state_index=inode_true_idx, init_probability=prob, init_tail=[(src_tag_true_comp, src_tag_true_state)])
        snode_2 = BKB_S_node(init_component_index=comp_idx, init_state_index=inode_false_idx, init_probability=1-prob, init_tail=[(src_tag_false_comp, src_tag_false_state)])
    else:
        snode_1 = BKB_S_node(init_component_index=comp_idx, init_state_index=inode_true_idx, init_probability=prob, init_tail=[(src_tag_true_comp, src_tag_true_state)] + other_tail)
        snode_2 = BKB_S_node(init_component_index=comp_idx, init_state_index=inode_false_idx, init_probability=1-prob, init_tail=[(src_tag_false_comp, src_tag_false_state)] + other_tail)

        bkb.addSNode(snode_1)
        bkb.addSNode(snode_2)
    return bkb
'''
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
    for comp_idx, state_idx in tqdm.tqdm(bottom_inodes, desc='Implementing Topological Strategy', leave=False):
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

def _addDemographicOption(option, bkb, src_population, src_population_data, option_dependencies=list(), include_src_tags=False):
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
    #bkb.makeGraph()

    #-- Create option dictionary
    options_dict = {bkb.getComponentName(comp_idx): bkb.getComponentINodeName(comp_idx, inode_true_idx)}

    #-- If no chain rule dependencies, just add prior s-nodes
    if len(option_dependencies) == 0:
        if include_src_tags:
            bkb, src_tag_true, _ = _processSrcTags(bkb, comp_idx, inode_true_idx, matched_srcs, src_population, set())
            bkb, src_tag_false, _ = _processSrcTags(bkb, comp_idx, inode_false_idx, matched_srcs, src_population, set())
            src_tag_true = [src_tag_true]
            src_tag_false = [src_tag_false]
        else:
            src_tag_true = list()
            src_tag_false = list()
        snode_1 = BKB_S_node(init_component_index=comp_idx, init_state_index=inode_true_idx, init_probability=prob, init_tail=src_tag_true)
        snode_2 = BKB_S_node(init_component_index=comp_idx, init_state_index=inode_false_idx, init_probability=1-prob, init_tail=src_tag_false)
        bkb.addSNode(snode_1)
        bkb.addSNode(snode_2)

    #-- Process Dependencies
    else:
        bkb = _processOptionDependency(option, option_dependencies, bkb, matched_srcs, src_population, src_population_data, include_src_tags)

    return bkb, options_dict, matched_srcs

def _linkSource(src_comp, src_state, non_src_comp, non_src_state, other_non_src_tails, prob,
                comp_idx, inode_true_idx, inode_false_idx, bkb, matched_srcs, src_hashs, src_hashs_inverse, processed_srcs):
    #-- Process the source name collection
    src_state_name = bkb.getComponentINodeName(src_comp, src_state)
    src_split = src_state_name.split('_')
    src_nums_str = src_split[0][1:-1]
    src_name_str = ''.join(src_split[1:])
    src_names = src_name_str.split(',')
    src_nums = [int(num) for num in src_nums_str.split(',')]
    #-- If this source component has format [num]_hash1,hash2, ... then it used fusion hack.
    if len(src_nums) != len(src_names):
        #-- Go through and process src nums by looking at hashes.
        src_nums = [src_hashs_inverse[int(src_hash)] for src_hash in src_names]

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
                                init_probability=prob,
                                init_tail=[(src_comp, src_true_idx)] + other_non_src_tails))
        if (src_comp, src_true_idx) not in processed_srcs:
            bkb.addSNode(BKB_S_node(init_component_index=src_comp,
                                    init_state_index=src_true_idx,
                                    init_probability=prob*true_prob,
                                    init_tail=[(comp_idx, inode_true_idx)]))
        processed_srcs.add((src_comp, src_true_idx))
    if (1 - true_prob) > 0:
        src_false_idx = bkb.addComponentState(src_comp, false_state_name)
        #-- Now connect with S nodes
        bkb.addSNode(BKB_S_node(init_component_index=non_src_comp,
                                init_state_index=non_src_state,
                                init_probability=prob,
                                init_tail=[(src_comp, src_false_idx)] + other_non_src_tails))
        if (src_comp, src_false_idx) not in processed_srcs:
            bkb.addSNode(BKB_S_node(init_component_index=src_comp,
                                    init_state_index=src_false_idx,
                                    init_probability=prob*(1 - true_prob),
                                    init_tail=[(comp_idx, inode_false_idx)]))
        processed_srcs.add((src_comp, src_false_idx))
    return bkb, processed_srcs

def _addSrcConnections(comp_idx, inode_true_idx, inode_false_idx, bkb, matched_srcs, src_hashes, src_hashs_inverse):
    #-- Attach to source nodes
    src_comp_indices = bkb.getSrcComponents()
    non_src_comp_indices = set(bkb.getAllComponentIndices()) - set(src_comp_indices)
    S_nodes_by_head = bkb.constructSNodesByHead()

    processed_srcs = set()

    for non_src_comp in tqdm.tqdm(non_src_comp_indices, desc='Linking Sources', leave=False):
        for non_src_state in bkb.getAllComponentINodeIndices(non_src_comp):
            try:
                snodes = S_nodes_by_head[non_src_comp][non_src_state]
            except KeyError:
                print('Warning: {} = {} has no incoming S nodes.'.format(bkb.getComponentName(non_src_comp), bkb.getComponentINodeName(non_src_comp, non_src_state)))
                continue
            for snode in snodes:
                prob = snode.probability
                src_tails = list()
                non_src_tails = list()
                for tail_idx in range(snode.getNumberTail()):
                    tail_comp, tail_state = snode.getTail(tail_idx)
                    #-- If this is a source component
                    if tail_comp in src_comp_indices:
                        src_tails.append((tail_comp, tail_state))
                    else:
                        non_src_tails.append((tail_comp, tail_state))
                had_illegal = False
                for src_tail_comp, src_tail_state in src_tails:
                    #-- If source state contains illegal values such as PatientX, etc, don't link.
                    src_name = bkb.getComponentINodeName(src_tail_comp, src_tail_state)
                    found_illegal = False
                    for illegal_src_val in ILLEGAL_SOURCE_STRINGS:
                        if illegal_src_val in src_name:
                            found_illegal = True
                            had_illegal = True
                            break
                    if found_illegal:
                        continue
                    prior_src_snode = S_nodes_by_head[src_tail_comp][src_tail_state][0]
                    #prob = prior_src_snode.probability
                    bkb, processed_srcs = _linkSource(src_tail_comp, src_tail_state, non_src_comp, non_src_state, non_src_tails, prob,
                                      comp_idx, inode_true_idx, inode_false_idx, bkb, matched_srcs, src_hashes, src_hashs_inverse, processed_srcs)
                    #-- delete the prior snode
                    try:
                        bkb.removeSNode(prior_src_snode)
                    except:
                        pass
                if len(src_tails) > 0 and not had_illegal:
                    #-- Delete other S node
                    try:
                        bkb.removeSNode(snode)
                    except:
                        pass
    return bkb
