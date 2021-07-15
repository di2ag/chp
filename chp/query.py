'''
Source code developed by DI2AG.
Thayer School of Engineering at Dartmouth College
Authors:    Dr. Eugene Santos, Jr
            Mr. Chase Yakaboski,
            Mr. Gregory Hyde,
            Dr. Keum Joo Kim
'''


import os
import sys
import pickle
import json
import csv
import io
import contextlib
import pandas as pd
import random

class Query:
    def __init__(self,
                 evidence=None,
                 targets=None,
                 marginal_evidence=None,
                 reasoning_type='updating',
                 name='query0',
                 dynamic_evidence=None,
                 dynamic_targets=None,
                 meta_evidence=None,
                 meta_targets=None,
                 ):
        """ BKB query class for CHP.

        Note:
            All evidence and target random variables and states including dynamic and meta must be exactly
            the same as the associated I-node names in the underlying BKB.

        kwargs:
            :param evidence: A dictionary of BKB evidence which has Random Variable
            names as keys and the associated Random Variable State names as the keys value.
            :type evidence: dict, optional
            :param targets: A list of Random Variable names that are targets of the BKB reasoning task.
            :type targets: list, defaults to None
            :param marginal_evidence: Not currently supported.
            :type marginal_evidence: None
            :param reasoning_type: The type of BKB reasoning to perform: updating or revision.
            :type reasoning_type: str, defaults to 'updating'.
            :param name: Name of the query.
            :type name: str, defaults to 'query0'.
            :param dynamic_evidence: A dictionary of dynamic BKB evidence which has Random Variable
            names as keys and a dictionary containing the keys: 'op' with values that can be ['>=', '==', '<=']
            and 'value' which is the associated value to apply to the 'op'.
            :type dynamic_evidence: dict, defaults to None
            :param dynamic_targets: A dictionary of dynamic BKB targets which has Random Variable
            names as keys and a dictionary containing the keys: 'op' with values that can be ['>=', '==', '<=']
            and 'value' which is the associated value to apply to the 'op'.
            :type dynamic_targets: dict, defaults to None
            :param meta_evidence: A dictionary of BKB evidence which has Random Variable names 
            appended to an underscore as keys and the associated Random Variable state name as the key value.
            :type meta_evidence: dict, defaults to None
            :param meta_targets: A list of Random Variable names append to an underscore that represent the meta targets
            of the BKB reasoning task.
            :type meta_targets: list, defaults to None
        """
        # Initialize evidence and targets
        self.evidence = evidence
        self.targets = targets
        self.marginal_evidence = marginal_evidence
        self.dynamic_evidence = dynamic_evidence
        self.dynamic_targets = dynamic_targets
        self.meta_evidence = meta_evidence
        self.meta_targets = meta_targets
        if evidence is None:
            self.evidence = {}
        if targets is None:
            self.targets = []
        if dynamic_evidence is None:
            self.dynamic_evidence = {}
        if dynamic_targets is None:
            self.dynamic_targets = {}
        if meta_evidence is None:
            self.meta_evidence = {}
        if meta_targets is None:
            self.meta_targets = []

        # Initiallize query parameters
        self.reasoning_type = reasoning_type
        self.name = name

        # Internal attributes that are set after query is run.
        self.result = None
        self.bkb = None
        self.compute_time = -1
        self.from_joint_reasoner = False

        #TODO: Likely deprecate below parameters.
        #self.independ_queries = None
        #self.independ_result = None
        #self.patient_data = None
        #self.interpolation = None
        #self.target_strategy = None
        #self.gene_var_direct = None
        #self.max_new_ev = None

    def make_bogus_updates(self):
        """ Makes random update probabilities based on dynamic_targets. Used for
        debugging only.

        Note:
            Does not actually perform any BKB updating.

        :return: Returns random probabilities for dynamic_targets.
        :rtype: dict
        """
        bogus_updates = {}
        for target in self.dynamic_targets:
            for state in ['True', 'False']:
                comp_name = ' '.join(target)
                prob = random.random()
                if comp_name not in bogus_updates:
                    bogus_updates[comp_name] = {state: prob}
                else:
                    bogus_updates[comp_name][state] = 1 - prob
        return bogus_updates

    def compose_evidence(self):
        """ Composes all the normal, dynamic, and meta evidence together into one dictionary.
        
        :return: A dictionary of the form of {[RandomVariableName]:  [RandomVariableState], ...}.
        :rtype: dict
        """
        evidence = self.evidence
        # Compose dynamic evidence
        for rv, evidence_dict in self.dynamic_evidence.items():
            op = evidence_dict["op"]
            value = evidence_dict["value"]
            evidence[rv] = '{} {}'.format(op, value)
        # Compose meta evidence
        for rv, state in self.meta_evidence.items():
            # Add underscore at beginning of rv name.
            evidence['_' + rv] = state
        return evidence

    def compose_targets(self):
        """ Composes all the normal, dynamic, and meta targets together into one list.
        
        :return: A list of composed targets from all target types.
        :rtype: dict
        """
        targets = self.targets
        # Compose dynamic targets
        for rv in self.dynamic_targets:
            targets.append(rv)
        # Compose meta targets
        for rv in self.meta_targets:
            # Add underscore at beginning of rv name.
            targets.append('_' + rv)
        return targets

    def add_evidence(self, random_variable, state):
        """ Add a single piece of evidence to the query.

        :param random_variable: A Random Variable name found in the BKB understudy.
        :type random_variable: str
        :param state: A Random Variable State name found in the BKB understudy.
        :type state: str
        """
        self.evidence[random_variable] = state

    def add_dynamic_evidence(self, random_variable, op, value):
        """ Add a single piece of dynamic evidence to the query.

        :param random_variable: A Random Variable name found in the BKB understudy.
        :type random_variable: str
        :param op: One of: '>=, '==', '<='.
        :type state: str
        :param value: The value to use with the operator.
        :type value: str
        """
        self.dynamic_evidence[random_variable] = {
                "op": op,
                "value": value,
                }

    def add_meta_evidence(self, random_variable, state):
        """ Add a single piece of meta evidence to the query.

        Note:
            Meta evidence is usually found in the BKB after some form of interpolation is done.

        :param random_variable: A Random Variable name found in the BKB understudy.
        :type random_variable: str
        :param state: A Random Variable State name found in the BKB understudy.
        :type state: str
        """
        self.meta_evidence[random_variable] = state
    
    def add_target(self, random_variable):
        """ Add a single target to the query.

        :param random_variable: A Random Variable name found in the BKB understudy.
        :type random_variable: str
        """
        self.targets.append(random_variable)

    def add_dynamic_target(self, random_variable, op, value):
        """ Add a single piece of dynamic targets to the query.

        :param random_variable: A Random Variable name found in the BKB understudy.
        :type random_variable: str
        :param op: One of: '>=, '==', '<='.
        :type state: str
        :param value: The value to use with the operator.
        :type value: str
        """
        self.dynamic_targets[random_variable] = {
                "op": op,
                "value": value,
                }

    def add_meta_target(self, random_variable):
        """ Add a single piece of meta targets to the query.

        Note:
            Meta targets is usually found in the BKB after some form of interpolation is done.

        :param random_variable: A Random Variable name found in the BKB understudy.
        :type random_variable: str
        :param state: A Random Variable State name found in the BKB understudy.
        :type state: str
        """
        self.meta_targets.append(random_variable)

""" Code to depreciate in next version.

    def save(self, directory, only_json=False):
        if not only_json:
            #-- Save a pickle file
            pickle_filename = os.path.join(directory, '{}.pk'.format(self.name))
            with open(pickle_filename, 'wb') as pickle_file:
                pickle.dump(file=pickle_file, obj=self)

            #-- Save out each piece in seperate files
            if self.evidence is not None:
                with open(os.path.join(directory, '{}.evid'.format(self.name)), 'w') as f_:
                    for comp_name, state_name in self.evidence.items():
                        f_.write('{},{}\n'.format(comp_name, state_name))

            if self.targets is not None:
                with open(os.path.join(directory, '{}.targ'.format(self.name)), 'w') as f_:
                    for comp_name in self.targets:
                        f_.write('{}\n'.format(comp_name))

            #-- Save out bkb
            if self.bkb is not None:
                self.bkb.save('{}.bkb'.format(self.name), pickle=True)

        #-- Save out JSON query info
        json_dict = {'evidence': self.evidence,
                     'targets': self.targets,
                     'type': self.type,
                     'dynamic_evidence': self.dynamic_evidence,
                     'dynamic_targets': self.dynamic_targets}

        #-- Potentially save out JSON Results
        if self.result is not None:
            inode_contrib = self.result.process_inode_contributions(include_srcs=False)
            result_dict = {'result': {'Updates': self.result.process_updates(),
                                      'Contributions': {' '.join(target): df.to_dict()
                                                        for target, df in self.result.contribs_to_dataframes(inode_contrib).items()},
                                      'Explanations': self.jsonExplanations(),
                                      'Report': self.getReportString()}}
            json_dict.update(result_dict)
        elif self.independ_result is not None:
            result_dict = {'result': {'Updates': self.independ_result,
                                      'Contributions': {' '.join(target): df.to_dict()
                                                        for target, df in self.getIndependentContributions().items()},
                                      'Explanations': self.jsonExplanations(),
                                      'Report': self.getReportString()}}

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

    def getIndependentContributions(self):
        independent_contribs = dict()
        contribs = list()
        for q_ in self.independ_queries:
            contribs.append(q_.result.process_inode_contributions(include_srcs=False))
        target_intersect = set.intersection(*[set(contrib.keys()) for contrib in contribs])
        for target in target_intersect:
            target_dict = dict()
            for contrib in contribs:
                for inode, cont in contrib[target].items():
                    if inode in independent_contribs:
                        target_dict[inode] += cont
                        target_dict[inode] /= 2
                    else:
                        target_dict[inode] = cont
            independent_contribs[target] = target_dict
        dfs = dict()
        for target, contrib_dict in independent_contribs.items():
            data_dict = {'I-node': list(), 'Contribution': list()}
            for inode, contrib in contrib_dict.items():
                data_dict['I-node'].append('{} = {}'.format(*inode))
                data_dict['Contribution'].append(contrib)
            df = pd.DataFrame(data=data_dict)
            df.sort_values(by=['Contribution'], inplace=True, ascending=False)
            df.set_index(['I-node'], inplace=True)
            dfs[target] = df
        return dfs

    def getExplanations(self, contributions_include_srcs=True, contributions_top_n_inodes=None, contributions_ignore_prefixes=None):
        explain_dict = dict()
        if self.independ_result is not None:
            explain_dict['Assumptions'] = 'Query assumes independence between genetic evidence.'
            explain_dict['Contributions Analysis'] = ['No sensitivity information is a available in the contributions field.']
            mostSigPatsIndepend = dict()
            query_patients = dict()
            for q in self.independ_queries:
                pat_dict = q.getSourcePatientAnalysis()
                for strat, target_data in pat_dict.items():
                    for target, data in target_data.items():
                        if strat not in query_patients:
                            query_patients[strat] = dict()
                        if data is not None:
                            if target in query_patients:
                                query_patients[strat][target].append(set([pat_id for pat_id in data]))
                            else:
                                query_patients[strat][target] = [set([pat_id for pat_id in data])]
                mostSigPatsIndepend.update(q.getSourcePatientAnalysis())
            pat_intersect = dict()
            for strat, target_pats in query_patients.items():
                if strat not in pat_intersect:
                    pat_intersect[strat] = dict()
                for target, pats in target_pats.items():
                    pat_intersect[strat][target] = set.intersection(*pats)
            for strat, target_data in mostSigPatsIndepend.items():
                for target, data in target_data.items():
                    if data is not None and strat == 'Most Significant Patients':
                        unsig_pats = set(data.keys()) - pat_intersect[strat][target]
                        for pat_id in unsig_pats: del mostSigPatsIndepend[target][pat_id]
            explain_dict['Patient Analysis'] = mostSigPatsIndepend
            explain_dict['Interpolation Strategy'] = 'Used gene and variant independence as interpolation strategy.'
            return explain_dict
        else:
            explain_dict['Assumptions'] = 'Query does not assume independence between genetic evidence.'
        inode_dict = self.result.process_inode_contributions(include_srcs=contributions_include_srcs,
                                                             top_n_inodes=contributions_top_n_inodes,
                                                             ignore_prefixes=contributions_ignore_prefixes,
                                                             remove_tuples=True)
        explain_dict['Contributions Analysis'] = inode_dict
        '''
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
            explain_dict['Contributions Analysis'].append(contrib_explain)
        '''
        explain_dict['Patient Analysis'] = self.getSourcePatientAnalysis()
        if self.interpolation is None or self.interpolation == 'standard':
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
        data = {'All Involved Patients': dict(), 'Most Significant Patients': dict()}
        for target, contrib_dict in inode_dict.items():
            target_str = '{} = {}'.format(target[0], target[1])
            src_hashs = list()
            for inode, contrib in contrib_dict.items():
                comp_name, state_name = inode
                if 'Source' in comp_name:
                    src_split1 = state_name.split('_')[-1]
                    src_split2 = src_split1.split(',')
                    sources = set()
                    for src_hash_str in src_split2:
                        try:
                            sources.add(int(src_hash_str))
                        except ValueError:
                            continue
                    if len(sources) > 0:
                        src_hashs.append(sources)
            #-- All involved Patients
            src_hashs_union = set.union(*src_hashs)
            #-- Most Significant Patients
            src_hashs_intersection = set.intersection(*src_hashs)
            if len(src_hashs_union) == 0:
                data['All Involved Patients'][target_str] = None
            else:
                data['All Involved Patients'][target_str] = {self.patient_data[src_hash]['patient_id']: self.patient_data[src_hash]
                                                      for src_hash in src_hashs_union}
            if len(src_hashs_intersection) == 0:
                data['Most Significant Patients'][target_str] = None
            else:
                data['Most Significant Patients'][target_str] = {self.patient_data[src_hash]['patient_id']: self.patient_data[src_hash]
                                                      for src_hash in src_hashs_intersection}
        return data

    def jsonExplanations(self, contributions_include_srcs=True, contributions_top_n_inodes=None, contributions_ignore_prefixes=None):
        explain = self.getExplanations(contributions_include_srcs=contributions_include_srcs,
                                       contributions_top_n_inodes=contributions_top_n_inodes,
                                       contributions_ignore_prefixes=contributions_ignore_prefixes)
        jsonSigPatients = dict()
        for analysis_type, infer_pat_data_dict in explain['Patient Analysis'].items():
            if infer_pat_data_dict is not None:
                jsonSigPatients[analysis_type] = dict()
                for target_str, pat_data_dict in infer_pat_data_dict.items():
                    jsonSigPatients[analysis_type][target_str] = dict()
                    if pat_data_dict is not None:
                        for patient_idx, data_dict in pat_data_dict.items():
                            jsonSigPatients[analysis_type][target_str][patient_idx] = dict()
                            if not data_dict is None:
                                for info_name, data in data_dict.items():
                                    if type(data) == tuple:
                                        data = list(data)
                                    jsonSigPatients[analysis_type][target_str][patient_idx][info_name] = data
        explain['Patient Analysis'] = jsonSigPatients
        return explain

    def printExplanations(self):
        explain = self.getExplanations()
        string = 'Assumptions: {}\n'.format(explain['Assumptions'])
        string += 'Sensitivity:\n'
        for sense in explain['Sensitivity']:
            string += '\t{}.\n'.format(sense)
        string += 'Patient Analysis:\n'
        for strat, target_patient_data in explain['Patient Analysis'].items():
            string += '{}:\n'.format(strat)
            for target, pat_data_dict in target_patient_data.items():
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

    def checkAndAdjustEvidence(self):
        #print(self.evidence)
        allVar = True
        #list of var types we want
        vars = list()
        for key in self.evidence.keys():
            if key[0:4] != 'var_':
                allVar = False
            vars.append(key[4:])
        #print(vars)
        if allVar and len(vars) != 0:
            geneVarFreq = list()
            with open(self.gene_var_direct, 'r') as csv_file:
                reader = csv.reader(csv_file)
                for row in reader:
                    geneVarFreq.append(row[0])
            count = 0
            newEvidenceDict = dict()
            for geneVar in geneVarFreq:
                varType = geneVar.split('-')[1]+'='
                #print(varType)
                if varType in vars:
                    newEvidenceDict['mut-var_'+geneVar+'='] = 'True'
                    count += 1
                if count == self.max_new_ev:
                    break
            self.evidence = newEvidenceDict
            return True
        else:
            return False
        #print(self.evidence)

    def getReportString(self):
        #-- Redirect sys.stdout to a string memory buffer.
        f = io.StringIO()
        with contextlib.redirect_stdout(f):
            self.getReport()
        return f.getvalue()

    def getReport(self):
        string = '---- Query Details -----\n'
        string += 'Demographic Evidence:\n'
        if self.dynamic_evidence is not None:
            for evid in self.dynamic_evidence:
                string += '\t{} {} {}\n'.format(evid[0], evid[1], evid[2])
        string += 'Evidence:\n'
        if self.evidence is not None:
            for rvName, stateName in self.evidence.items():
                string += '\t{} = {}\n'.format(rvName, stateName)
        string += 'Targets:\n'
        if self.targets is not None:
            for target in self.targets:
                string += '\t{}\n'.format(target)
        print(string)
        if self.result is not None:
            self.result.summary(include_srcs=False)
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
            print('Independent Contributions:')
            for target, df in self.getIndependentContributions().items():
                print(target)
                print(df)
            print('------ Explanations ------')
            print(self.printExplanations())
            print('--------Individual Query Reports -------')
            for q in self.independ_queries:
                q.getReport()
        else:
            print('No results found.')
"""
