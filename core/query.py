import os
import sys


class Query:
    def __init__(self, evidence, targets=None, marginal_evidence=None, type='revision', name='query0', meta_evidence=None):
        self.evidence = evidence
        self.targets = targets
        self.marginal_evidence = marginal_evidence
        self.type = type
        self.name = name
        self.meta_evidence = meta_evidence

    def build_bkb_query(self):
        raise NotImplementedError

