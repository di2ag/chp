import os
import sys
import pickle
import json

class Query:
    def __init__(self, evidence=dict(), targets=list(), marginal_evidence=None, type='updating', name='query0', meta_evidence=None, meta_targets=None):
        self.evidence = evidence
        self.targets = targets
        self.marginal_evidence = marginal_evidence
        self.type = type
        self.name = name
        self.meta_evidence = meta_evidence
        self.meta_targets = meta_targets
        self.result = None
        self.bkb = None
        self.independ_queries = None
        self.independ_result = None
        self.compute_time = -1

    def save(self, directory, only_json=False):
        if not only_json:
            #-- Save a pickle file
            pickle_filename = os.path.join(directory, '{}.pk'.format(self.name))
            with open(pickle_filename, 'wb') as pickle_file:
                pickle.dump(file=pickle_file, obj=self)

            #-- Save out each piece in seperate files
            with open(os.path.join(directory, '{}.evid'.format(self.name)), 'w') as f_:
                for comp_name, state_name in self.evidence.items():
                    f_.write('{},{}\n'.format(comp_name, state_name))

            with open(os.path.join(directory, '{}.targ'.format(self.name)), 'w') as f_:
                for comp_name in self.targets:
                    f_.write('{}\n'.format(comp_name))

            #-- Save out bkb
            if self.bkb is not None:
                self.bkb.save('{}.bkb'.format(self.name))

        #-- Save out JSON query info
        json_dict = {'evidence': self.evidence,
                     'targets': self.targets,
                     'type': self.type,
                     'meta_evidence': self.meta_evidence,
                     'meta_targets': self.meta_targets}

        #-- Potentially save out JSON Results
        if self.result is not None:
            inode_contrib = self.result.process_inode_contributions()
            result_dict = {'result': {'Updates': self.result.process_updates(),
                                      'Contributions': {' '.join(target): df.to_dict()
                                                        for target, df in self.result.contribs_to_dataframes(inode_contrib).items()},
                                      'Explanations': self.getExplanations()}}
            json_dict.update(result_dict)
        json_file = os.path.join(directory, '{}.json'.format(self.name))
        with open(json_file, 'w') as f_:
            json.dump(json_dict, f_)

        if only_json:
            return json_file
        else:
            return pickle_filename, json_file

    def read(self, query_file, file_format='pickle'):
        if file_format == 'pickle':
            with open(query_file, 'rb') as qf_:
                return pickle.load(qf_)
        elif file_format == 'json':
            with open(query_file, 'r') as qf_:
                query_dict = json.load(qf_)
            return Query(**query_dict)
        else:
            raise ValueError('Unrecognized file format: {}'.format(file_format))

    def getExplanations(self):
        explain_dict = dict()
        if self.independ_result is not None:
            explain_dict['Assumptions'] = 'Query assumes independence between genetic evidence.'
        else:
            explain_dict['Assumptions'] = 'Query does not assume independence between genetic evidence.'
        inode_dict = self.result.process_inode_contributions()
        explain_dict['Sensitivity'] = list()
        for target, contrib_dict in inode_dict.items():
            target_str = ' '.join(target)
            most_sig_inodes = list()
            max_contrib = -1
            for inode, contrib in contrib_dict.items():
                inode_str = ' '.join(inode)
                if contrib > max_contrib:
                    most_sig_inodes = [inode_str]
                    max_contrib = contrib
                elif contrib == max_contrib:
                    most_sig_inodes.append(inode_str)
                else:
                    continue
            contrib_explain = 'The most sensitivity variables for {} are {}'.format(target_str,
                                                                                    ', '.join(most_sig_inodes))
            explain_dict['Sensitivity'].append(contrib_explain)
        explain_dict['MostSignificantPatients'] = ['Unknown']
        return explain_dict

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
        elif self.independ_result is not None:
            print('---- Results Using Independence Assumption -----')
            for update, state_dict in self.independ_result.items():
                print('\t{}'.format(update))
                for state, prob in state_dict.items():
                    print('\t\t{} = {}'.format(state, prob))
        else:
            print('No results found.')
