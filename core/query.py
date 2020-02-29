import os
import sys
import pickle


class Query:
    def __init__(self, evidence=dict(), targets=list(), marginal_evidence=None, type='revision', name='query0', meta_evidence=None, meta_targets=None):
        self.evidence = evidence
        self.targets = targets
        self.marginal_evidence = marginal_evidence
        self.type = type
        self.name = name
        self.meta_evidence = meta_evidence
        self.meta_targets = meta_targets
        self.result = None
        self.bkb = None
        self.compute_time = -1

    def save(self, directory):
        #-- Save a pickle file
        with open(os.path.join(directory, '{}.pk'.format(self.name)), 'wb') as pickle_file:
            pickle.dump(file=pickle_file, obj=self)

        #-- Save out each piece in seperate files
        with open(os.path.join(directory, '{}.evid'.format(self.name)), 'w') as f_:
            for comp_name, state_name in self.evidence.items():
                f_.write('{},{}\n'.format(comp_name, state_name))

        with open(os.path.join(directory, '{}.targ'.format(self.name)), 'w') as f_:
            for comp_name in self.targets:
                f_.write('{}\n'.format(comp_name))

        #-- Save out bkb
        self.bkb.save('{}.bkb'.format(self.name))

    def getReport(self):
        string = '---- Query Details -----\n'
        string += 'Demographic Evidence:\n'
        if self.meta_evidence is not None:
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
