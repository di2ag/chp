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
        self.result = None
        self.bkb = None
        self.compute_time = -1

    def getReport(self):
        string = '---- Query Details -----\n'
        string += 'Demographic Evidence:\n'
        for evid in self.meta_evidence:
            string += '\t{} {} {}\n'.format(evid[0], evid[1], evid[2])
        string += 'Evidence:\n'
        for rvName, stateName in self.evidence.items():
            string += '\t{} = {}\n'.format(rvName, stateName)
        string += 'Targets:\n'
        for target in self.targets:
            string += '\t{}\n'.format(target)
        print(string)
        if self.result is not None:
            self.result.summary()
            print('Computed in {} sec.'.format(self.compute_time))
        else:
            print('No results found.')
