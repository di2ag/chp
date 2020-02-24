import os
import sys
import pickle
import operator
import tqdm
import copy
from concurrent.futures import ProcessPoolExecutor, wait
import time

from pybkb.core.cpp_base.reasoning import revision, updating
from pybkb.core.common.bayesianKnowledgeBase import BKB_component, BKB_I_node, BKB_S_node

class Reasoner:
    def __init__(self, fused_bkb, patients):
        self.fused_bkb = fused_bkb
        self.patients = patients

        #-- Preprocess src hash values
        src_components = fused_bkb.getSrcComponents()
        src_hashs = dict()
        for component in src_components:
            for idx in range(component.getNumberStates()):
                state = component.getState(idx)
                #-- Collect src numbers and string name or hash
                src_num = int(state.name.split('_')[0][1:-1])
                try:
                    src_hashs[src_num] = int(''.join(state.name.split('_')[1:]))
                except ValueError:
                    src_hashs[src_num] = state.name
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
        if query.type == 'revision':
            res = revision(query.bkb,
                           query.evidence,
                           marginal_evidence=query.marginal_evidence,
                           targets=query.targets,
                           file_prefix=query.name)
        elif query.type == 'updating':
            res = updating(query.bkb,
                           query.evidence,
                           marginal_evidence=query.marginal_evidence,
                           targets=query.targets,
                           file_prefix=query.name)
        else:
            raise ValueError('Unreconginzed reasoning type: {}.'.format(query.type))
        compute_time = time.time() - start_time
        query.result = res
        query.compute_time = compute_time
        return query

    def process_metaVariables(self, bkb, meta_variables):
        #-- Get all src components
        src_components = bkb.getSrcComponents()

        #-- Collect source hashs that match metadata as well as population stats
        transformed_meta = dict()
        pop_stats = dict()
        for meta in meta_variables:
            matched_srcs = list()
            prop, cond, val = meta
            if cond == '==':
                op = operator.eq
            elif cond == '>=':
                op = operator.ge
                val = float(val)
            elif cond == '<=':
                op = operator.le
                val = float(val)
            else:
                raise ValueError('Unknown Condition')
            count_true = 0
            count_false = 0
            for src_num, src_hash in self.src_hashs.items():
                try:
                    if op(type(val)(self.metadata[src_hash][prop]), val):
                        matched_srcs.append(src_num)
                        count_true += 1
                    else:
                        count_false += 1
                except ValueError:
                    continue
                except KeyError:
                    continue
            total = count_true + count_false
            pop_stats[meta] = (count_true/total, count_false/total)

            #-- Build Augmented BKB with attach metadata using frequency priors
            #-- Construct new I-node for each meta evidence.
            meta_component = BKB_component('{} {} {}'.format(prop, cond, val))
            inode_true = BKB_I_node(init_name='True', init_component=meta_component)
            inode_false = BKB_I_node(init_name='False', init_component=meta_component)
            #meta_component.addINode(inode_true)
            #meta_component.addINode(inode_false)

            #-- Add component to BKB
            bkb.addComponent(meta_component)
            bkb.addComponentState(meta_component, inode_true)
            bkb.addComponentState(meta_component, inode_false)

            #print(meta_component.states)
            #print(meta_component.getStateIndex(inode_true))

            #-- Create prior S-nodes
            bkb.addSNode(BKB_S_node(init_component=meta_component, init_state=inode_true, init_probability=pop_stats[meta][0]))
            bkb.addSNode(BKB_S_node(init_component=meta_component, init_state=inode_false, init_probability=pop_stats[meta][0]))

            #-- Link source nodes
            for component in tqdm.tqdm(src_components, desc='Linking {} Source Nodes'.format(prop)):
                for idx in range(component.getNumberStates()):
                    state = component.getState(idx)
                    cidx = bkb.getComponentIndex(component)
                    iidx = component.getStateIndex(state)
                    #-- There should only be a single s-node as these are source reliabilities.
                    src_s_node = bkb.S_nodes_by_head[cidx][iidx][0]
                    src_reliab = src_s_node.probability
                    bkb.removeSNode(src_s_node)

                    src_num = int(state.name.split('_')[0][1:-1])
                    if src_num in matched_srcs:
                        bkb.addSNode(BKB_S_node(init_component=component, init_state=state, init_probability=src_reliab, init_tail=[(meta_component, inode_true)]))
                        bkb.addSNode(BKB_S_node(init_component=component, init_state=state, init_probability=0, init_tail=[(meta_component, inode_false)]))
                    else:
                        bkb.addSNode(BKB_S_node(init_component=component, init_state=state, init_probability=0, init_tail=[(meta_component, inode_true)]))
                        bkb.addSNode(BKB_S_node(init_component=component, init_state=state, init_probability=src_reliab, init_tail=[(meta_component, inode_false)]))

            #-- Transform meta evidence to standard BKB evidence using new I_node.
            transformed_meta[meta_component.name] = inode_true.name

        return transformed_meta, bkb

    def analyze_query(self, query):
        #-- Copy BKB
        bkb = copy.deepcopy(self.fused_bkb)

        if query.meta_evidence is not None:
            print('Processing Demographic Evidence...')
            transformed_meta_evidence, bkb = self.process_metaVariables(bkb, query.meta_evidence)

        if query.meta_targets is not None:
            print('Processing Demographic Targets...')
            transformed_meta_targets, bkb = self.process_metaVariables(bkb, query.meta_targets)

        #-- Conduct Reasoning
        query.evidence.update(transformed_meta_evidence)
        for meta_name in transformed_meta_targets:
            query.targets.append(meta_name)
        query.bkb = bkb
        return self.solve_query(query)
