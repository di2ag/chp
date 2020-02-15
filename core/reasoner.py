import os
import sys

from pybkb.core.cpp_base.reasoning import revision, updating

class Reasoner:
    def __init__(self, fused_bkb):
        self.fused_bkb = fused_bkb


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
        result = self.solve_query(query)
        if result is not None:
            result.summary()
        return result
