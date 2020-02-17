import os
import sys
import pickle
import operator
import tqdm
import copy
from concurrent.futures import ProcessPoolExecutor, wait

from pybkb.core.cpp_base.reasoning import revision, updating

class Reasoner:
    def __init__(self, fused_bkb, patients):
        self.fused_bkb = fused_bkb
        self.patients = patients
        #-- Preprocess patient metadata
        src_components = fused_bkb.getSrcComponents()
        src_hashs = []
        for component in src_components:
            src_hashs.extend([int(state.name.split('_')[-1]) for state in component.states])
        self.src_hashs = set(src_hashs)

    #-- Should be source hash value followed by a dictionary of all available meta data. The file is assumed to be a pickle.
    def set_src_metadata(self, metadata_file):
        with open(metadata_file, 'rb') as m_:
            self.metadata = pickle.load(m_)
        print('Available Metadata:')
        self.metadata_labels = list()
        for hash_key in self.metadata:
            self.metadata_labels.extend(list(self.metadata[hash_key].keys()))
        for meta in set(self.metadata_labels):
            print('\t{}'.format(meta))


    def solve_query(self, query):
        if query.type == 'revision':
            res = revision(self.fused_bkb,
                           query.evidence,
                           marginal_evidence=query.marginal_evidence,
                           targets=query.targets,
                           file_prefix=query.name)
        elif query.type == 'updating':
            res = updating(self.fused_bkb,
                           query.evidence,
                           marginal_evidence=query.marginal_evidence,
                           targets=query.targets,
                           file_prefix=query.name)
        else:
            raise ValueError('Unreconginzed reasoning type: {}.'.format(query.type))

        return res

    def analyze_query(self, query):
        if query.meta_evidence is not None:
            #-- Collect source hashs that match metadata
            matched_srcs = list()
            for meta in query.meta_evidence:
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
                for hash_key in self.src_hashs:
                    try:
                        if op(type(val)(self.metadata[hash_key][prop]), val):
                            matched_srcs.append(hash_key)
                    except ValueError:
                        continue
            src_evidence = list()
            src_components = self.fused_bkb.getSrcComponents()
            print(matched_srcs)
            for src_hash in tqdm.tqdm(set(matched_srcs)):
                evid_ = dict()
                for component in src_components:
                    i_node = self.fused_bkb.findINode(component, str(src_hash), contains=True)
                    if i_node != -1:
                        evid_[component.name] = i_node.name
                if len(evid_.keys()) > 0:
                    src_evidence.append(evid_)

        #-- Conduct Reasoning
        prob = 0
        futures = []
        pool = ProcessPoolExecutor(30)
        print('Starting Pool Workers.')
        for i, src_evid in enumerate(src_evidence):
            new_query = copy.copy(query)
            new_query.evidence.update(src_evid)
            new_query.name = 'q' + str(i)
            futures.append(pool.submit(self.solve_query, new_query))

        results = list()
        done = wait(futures)[0]
        for fut in done:
            res = fut.result()
            prob += res.result[0]['Prob']
            results.append(res)

        return results
