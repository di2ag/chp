import os
import sys
import pickle
import json

class Query:
    def __init__(self, evidence=dict(),
                 targets=list(),
                 marginal_evidence=None,
                 type='updating',
                 name='query0',
                 meta_evidence=None,
                 meta_targets=None):
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
        self.patient_data = None
        self.interpolation = None
        self.target_strategy = None

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
                                      'Explanations': self.jsonExplanations()}}
        elif self.independ_result is not None:
            result_dict = {'result': {'Updates': self.independ_result,
                                      'Contributions': None,
                                      'Explanations': self.jsonExplanations()}}

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
            explain_dict['Sensitivity'] = ['No sensitivity information for results with independence assumption yet, try correlation.']
            mostSigPatsIndepend = dict()
            for q in self.independ_queries:
                mostSigPatsIndepend.update(q.getSourcePatientAnalysis())
            explain_dict['Most Significant Patients'] = mostSigPatsIndepend
            explain_dict['Interpolation Strategy'] = 'Used gene and variant independence as interpolation strategy.'
            return explain_dict
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
                if 'Source' in inode_str or target_str in inode_str or 'Collection' in inode_str:
                    continue
                if contrib > max_contrib:
                    most_sig_inodes = [inode_str]
                    max_contrib = contrib
                elif contrib == max_contrib:
                    most_sig_inodes.append(inode_str)
                else:
                    continue
            contrib_explain = 'The most sensitive variables for {} are {}'.format(target_str,
                                                                                    ', '.join(most_sig_inodes))
            explain_dict['Sensitivity'].append(contrib_explain)
        explain_dict['Most Significant Patients'] = self.getSourcePatientAnalysis()
        if self.interpolation is None:
            explain_dict['Interpolation Strategy'] = 'No interpolation stradegy was used.'
        elif self.interpolation == 'independence':
            explain_dict['Interpolation Strategy'] = 'Independent interpolation strategy was such that all genetic pieces of evidence \
                    and the product of their probabilities were used in determing the posterior distribtion on the targets.'
        elif self.interpolation == 'correlation':
            explain_dict['Interpolation Strategy'] = dict()
            explain_dict['Interpolation Strategy']['Description'] = 'Interpolated using mutation correlations.'
            explain_dict['Interpolation Strategy']['Correlation Contribution Chain'] = self.getPatientXExplanations()
        else:
            raise ValueError('Interpolation stradegy: {} was not recognized.'.format(self.interpolation))
        return explain_dict

    def chainSearch(self, head, chain, target, snode_contribs):
        for tail_list, contrib in snode_contribs[target][head].items():
            for tail in tail_list:
                tail_comp, tail_state = tail
                if '_mut_' in tail_comp and 'Source' not in tail_comp:
                    chain.append((tail_comp, contrib))
                    self.chainSearch(tail, chain, target, snode_contribs)
        return chain

    def getPatientXExplanations(self):
        snode_contribs = self.result.process_contributions()
        patientX_chain = dict()
        for target, head_tail_contrib in snode_contribs.items():
            patientX_chain[target] = dict()
            for evid_comp, evid_state in self.evidence.items():
                chain = self.chainSearch((evid_comp, evid_state), list(), target, snode_contribs)
                if len(chain) > 0:
                    patientX_chain[target][(evid_comp, evid_state)] = chain
        return patientX_chain

    def getSourcePatientAnalysis(self):
        inode_dict = self.result.process_inode_contributions()
        data = dict()
        for target, contrib_dict in inode_dict.items():
            src_hashs = set()
            for inode, contrib in contrib_dict.items():
                comp_name, state_name = inode
                if 'Source' in comp_name:
                    src_split1 = state_name.split('_')[-1]
                    src_split2 = src_split1.split(',')
                    for src_hash_str in src_split2:
                        try:
                            src_hashs.add(int(src_hash_str))
                        except ValueError:
                            continue
            if len(src_hashs) == 0:
                data[target] = None
            else:
                data[target] = {self.patient_data[src_hash]['Patient_ID']: self.patient_data[src_hash] for src_hash in src_hashs}
        return data

    def jsonExplanations(self):
        explain = self.getExplanations()
        jsonSigPatients = dict()
        for target, pat_data_dict in explain['Most Significant Patients'].items():
            if pat_data_dict is not None:
                target_str = '{} = {}'.format(target[0], target[1])
                jsonSigPatients[target_str] = dict()
                for patient_idx, data_dict in pat_data_dict.items():
                    jsonSigPatients[target_str][patient_idx] = dict()
                    for info_name, data in data_dict.items():
                        if type(data) == tuple:
                            data = list(data)
                        jsonSigPatients[target_str][patient_idx][info_name] = data
        explain['Most Significant Patients'] = jsonSigPatients
        return explain

    def printExplanations(self):
        explain = self.getExplanations()
        string = 'Assumptions: {}\n'.format(explain['Assumptions'])
        string += 'Sensitivity:\n'
        for sense in explain['Sensitivity']:
            string += '\t{}.\n'.format(sense)
        string += 'Most Significant Patient Information:\n'
        for target, pat_data_dict in explain['Most Significant Patients'].items():
            if pat_data_dict is not None:
                string += '{}\n'.format(target)
                for patient_id, data_dict in pat_data_dict.items():
                    string += '\t{}:\n'.format(patient_id)
                    for info_name, data in data_dict.items():
                        if type(data) == list or type(data) == tuple:
                            data_str = ','.join([str(val) for val in data])
                            if len(data_str) > 100:
                                data_str = data_str[:100] + '...'
                                string += '\t\t{} = {}\n'.format(info_name, data_str)
                        elif type(data) == dict:
                            continue
                        else:
                            string += '\t\t{} = {}\n'.format(info_name, data)
        string += 'Interpolation Strategy Details:\n'
        if type(explain['Interpolation Strategy']) == dict:
            for label, info in explain['Interpolation Strategy'].items():
                if type(info) == dict:
                    string += '{}:\n'.format(label)
                    for target, evid_chain in info.items():
                        string += '\tTarget: {} = {}\n'.format(target[0], target[1])
                        for evid, chain in evid_chain.items():
                            string += '\t\tEvidence: {} = {}\n'.format(evid[0], evid[1])
                            string += '\t\tChain: {}\n'.format(chain)
                else:
                    string += '{}: {}\n'.format(label, info)
        else:
            string += str(explain['Interpolation Strategy']) + '\n'
        return string

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
            print('------ Explanations ------')
            print(self.printExplanations())
        elif self.independ_result is not None:
            print('---- Results Using Independence Assumption -----')
            print('Total Result:')
            for update, state_dict in self.independ_result.items():
                print('\t{}'.format(update))
                for state, prob in state_dict.items():
                    print('\t\t{} = {}'.format(state, prob))
            print('------ Explanations ------')
            print(self.printExplanations())
            print('--------Individual Query Reports -------')
            for q in self.independ_queries:
                q.getReport()
        else:
            print('No results found.')
