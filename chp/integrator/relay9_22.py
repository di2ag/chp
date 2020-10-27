'''
Source code developed by DI2AG.
Thayer School of Engineering at Dartmouth College
	Authors:    Dr. Eugene Santos, Jr
            Mr. Chase Yakaboski,
            Mr. Gregory Hyde,
            Dr. Keum Joo Kim
'''

import copy
import csv
import sys
import uuid
import pickle

from chp.query import Query
from chp.reasoner import Reasoner

from chp_data.bkb_handler import BkbDataHandler

class Relay9_22:
    def __init__(self, query, hosts_filename=None, num_processes_per_host=0):
        # query graph components
        self.query = query
        self.qg = self.query['query_graph']
        if 'knowledge_graph' not in list(self.query.keys()):
            self.kg = { "edges": [],
                        "nodes": []
                      }
        else:
            self.kg = self.query['knowledge_graph']
        if 'results' not in list(self.query.keys()):
            self.results = { "node_bindings": [],
                             "edge_bindings": []
                           }
        else:
            self.results = self.query['results']

        # Instiatate Reasoner
        self.bkb_data_handler = BkbDataHandler(dataset_version='1.2')
        self.reasoner = Reasoner(bkb_data_handler=self.bkb_data_handler,
                                hosts_filename=hosts_filename,
                                num_processes_per_host=num_processes_per_host)

        # prepare curie gene dict
        self.true_gene_contrib = dict()
        self.false_gene_contrib = dict()
        self.gene_curie_dict = dict()
        self.gene_to_curie = dict()
        with open(self.bkb_data_handler.gene_curie_path, 'r') as gene_file:
            reader = csv.reader(gene_file)
            next(reader)
            for row in reader:
                self.gene_curie_dict[row[1]] = row[0]
                self.true_gene_contrib[row[0]] = 0
                self.false_gene_contrib[row[0]] = 0
                self.gene_to_curie[row[0]] = row[1]
        # prepare curie drug dict
        self.drug_curie_dict = dict()
        with open(self.bkb_data_handler.drug_curie_path, 'r') as drug_file:
            reader = csv.reader(drug_file)
            next(reader)
            for row in reader:
                self.drug_curie_dict[row[1]] = row[0]

        # default query specification
        self.target_strategy = 'explicit'
        self.interpolation = 'standard'

    ##########################################################
    # checkQuery
    # Input:
    # Output: boolean True/False
    #--------------------------------------------------------
    # Description: NOT IMPLEMENTED. Checks the query graph
    # to ensure it meets our specifications

    def checkQuery(self):
        return True

    ##########################################################
    # buildQueries
    # Input:
    # Output: built query
    #--------------------------------------------------------
    # Description: Parses over sent query graph to form a BKB
    # query. The BKB query is returned.

    def buildQueries(self):

        evidence = dict()
        targets = list()
        acceptable_target_curies = ['EFO:0000714']
        for node in self.qg['nodes']:
            if node['type'] == 'Drug':
                drug_curie = node['curie']
                try:
                    self.drug = self.drug_curie_dict[drug_curie]
                except:
                    sys.exit('Invalid CHEMBL Identifier. Must be in form CHEMBL:<ID>')
                evidence['drug_{}'.format(self.drug)] = 'True'
        for edge in self.qg['edges']:
            if edge['type'] == 'chemical_to_disease_or_phenotypic_feature_association':
                if 'value' in edge.keys():
                    self.days = edge['value']
                else:
                    self.days = 970

        targets.append(('Survival_Time', '>=', self.days))

        query = Query(evidence=evidence,
                      targets=[],
                      meta_evidence=None,
                      meta_targets=targets,
                      type='updating')
        self.chp_query = query
        return query

    ##########################################################
    # runQueries
    # Input:
    # Output:
    #--------------------------------------------------------
    # Description: Runs built BKB query to calculate probability
    # of survival. Outputs are checked for -1s. If BOTH are -1
    # they are left alone (should interpret as no inference and
    # should document in openAPI). Otherwise probabilities are
    # normalized

    def runQueries(self):
        query = self.reasoner.analyze_query(copy.deepcopy(self.chp_query),
                                            save_dir=None,
                                            target_strategy=self.target_strategy,
                                            interpolation=self.interpolation
                                            )
        self.target_info = []
        for update, prob in query.result.updates.items():
            comp_idx, state_idx = update
            comp_name = query.bkb.getComponentName(comp_idx)
            state_name = query.bkb.getComponentINodeName(comp_idx, state_idx)
            #print(comp_name, state_name, prob)
            self.target_info.append([comp_name, state_name, prob])

        if self.target_info[0][1] == 'True':
            self.truth_assignment = self.target_info[0][2]
            self.false_assignment = self.target_info[1][2]
        else:
            self.truth_assignment = self.target_info[1][2]
            self.false_assignment = self.target_info[0][2]

        if self.truth_assignment != -1 and self.false_assignment != -1:
            prob_sum = self.truth_assignment + self.false_assignment
            self.truth_assignment /= prob_sum
            self.false_assignment /= prob_sum
            sensitivities = True
        elif self.truth_assignment == -1 and self.false_assignment != -1:
            self.truth_assignment = 0
            prob_sum = self.truth_assignment + self.false_assignment
            self.truth_assignment /= prob_sum
            self.false_assignment /= prob_sum
        elif self.truth_assignment != -1 and self.false_assignment == -1:
            self.false_assignment = 0
            prob_sum = self.truth_assignment + self.false_assignment
            self.truth_assignment /= prob_sum
            self.false_assignment /= prob_sum

        report = query.jsonExplanations(contributions_include_srcs=False,
                                        contributions_top_n_inodes=10,
                                        contributions_ignore_prefixes=['_'])
        self.report = {'Patient Analysis': report['Patient Analysis'],
                       'Contribution Analysis': report['Contributions Analysis']}
        # total contrib values
        true_contrib = self.report['Contribution Analysis']['Survival_Time >= {} = True'.format(self.days)]['drug_{} = True'.format(self.drug)]
        false_contrib = self.report['Contribution Analysis']['Survival_Time >= {} = False'.format(self.days)]['drug_{} = True'.format(self.drug)]
        # total patients in contrib cat
        true_pats = len(self.report['Patient Analysis']['All Involved Patients']['Survival_Time >= {} = True'.format(self.days)].keys())
        false_pats = len(self.report['Patient Analysis']['All Involved Patients']['Survival_Time >= {} = False'.format(self.days)].keys())
        # true_individual contrib
        true_ind_cont = float(true_contrib)/float(true_pats)
        false_ind_cont = float(false_contrib)/float(false_pats)

        self.gene_list_contrib = dict()


        patient_dict = pickle.load(open(self.bkb_data_handler.patient_data_pk_path, 'rb'))
        for key in patient_dict:
            pat = patient_dict[key]
            if pat['Survival_Time'] > self.days and self.drug in pat['Drug_Name(s)']:
                for gene in pat['Patient_Genes']:
                    self.true_gene_contrib[gene] += true_ind_cont
            elif pat['Survival_Time'] < self.days and self.drug in pat['Drug_Name(s)']:
                for gene in pat['Patient_Genes']:
                    self.false_gene_contrib[gene] += false_ind_cont
        self.report['Contribution Analysis']['Survival_Time >= {} = True'.format(self.days)] = {'{}-{}'.format(k,self.gene_to_curie[k]): v for k,v in sorted(self.true_gene_contrib.items(), key=lambda item: item[1], reverse=True)}
        self.report['Contribution Analysis']['Survival_Time >= {} = False'.format(self.days)] = {'{}-{}'.format(k,self.gene_to_curie[k]): v for k,v in sorted(self.false_gene_contrib.items(), key=lambda item: item[1], reverse=True)}


    ##########################################################
    # constructDecoratedKG
    # Input:
    # Output: endpoint query response
    #--------------------------------------------------------
    # Description: knowledge graph is a copy of query graph where
    # we add a edge property, 'has_confidence_level' annotated with
    # our determined value. The whole query response is formed and
    # returned.

    def constructDecoratedKG(self):
        self.kg = copy.deepcopy(self.qg)
        results = {'edge_bindings':[], 'node_bindings':[]}
        # update target node info and form edge pair combos for results graph
        edge_pairs = list()
        for edge in self.kg['edges']:
            if edge['type'] == 'chemical_to_disease_or_phenotypic_feature_association':
                edge['has_confidence_level'] = self.truth_assignment
                #they may want descriptions later
                edge['Description'] = self.report['Contribution Analysis']
            qg_id = edge['id']
            kg_id = str(uuid.uuid4())
            edge_pairs.append([qg_id,kg_id])
            edge['id'] = kg_id
        # form node pair combos for results graph
        node_pairs = list()
        for node in self.kg['nodes']:
            qg_id = node['id']
            kg_id = str(uuid.uuid4())
            node_pairs.append([qg_id,kg_id])

        # update edge/node bindings in results graph
        for edge_pair in edge_pairs:
            results['edge_bindings'].append({'qg_id':edge_pair[0], 'kg_id':edge_pair[1]})
        for node_pair in node_pairs:
            results['node_bindings'].append({'qg_id':node_pair[0], 'kg_id':node_pair[1]})

        # query response
        reasoner_std = {'query_graph': self.qg,
                        'knowledge_graph': self.kg,
                        'results': results}
        return reasoner_std
