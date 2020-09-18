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
import logging
import uuid
import statistics

from chp.query import Query
from chp.reasoner import Reasoner

from chp_data.bkb_handler import BkbDataHandler

#-- Setup logging
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
DEBUG = True

class UnsecretHandler:
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

        #-- Instiatate Reasoner
        self.bkb_data_handler = BkbDataHandler(dataset_version='1.1')
        self.reasoner = Reasoner(bkb_data_handler=self.bkb_data_handler,
                                hosts_filename=hosts_filename,
                                num_processes_per_host=num_processes_per_host)
        # prepare curie gene dict
        self.gene_curie_dict = dict()
        with open(self.bkb_data_handler.gene_curie_path, 'r') as gene_file:
            reader = csv.reader(gene_file)
            next(reader)
            for row in reader:
                self.gene_curie_dict[row[1]] = row[0]
        # prepare curie drug dict
        self.drug_curie_dict = dict()
        with open(self.bkb_data_handler.drug_curie_path, 'r') as drug_file:
            reader = csv.reader(drug_file)
            next(reader)
            for row in reader:
                self.drug_curie_dict[row[1]] = row[0]

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

        # ensure we are using all nodes/edges
        total_nodes = 0
        total_edges = 0

        # get phenotype node
        targets = list()
        acceptable_target_curies = ['EFO:0000714']
        for node in self.qg['nodes']:
            if node['type'] == 'PhenotypicFeature' and node['curie'] in acceptable_target_curies:
                target_id = node['id']
                total_nodes += 1
        if total_nodes == 0:
            acceptable_target_curies_print = ','.join(acceptable_target_curies)
            sys.exit('Survival Node not found. Node type muse be \'PhenotypicFeature\' and curie must be in: ' + acceptable_target_curies_print)
        elif total_nodes > 1:
            sys.exit('Too many target nodes')

        # get disease node info and ensure only 1 disease:
        acceptable_disease_curies = ['MONDO:0007254']
        for node in self.qg['nodes']:
            if node['type'] == 'disease' and node['curie'] in acceptable_disease_curies:
                disease_id = node['id']
                for edge in self.qg['edges']:
                    if edge['type'] == 'disease_to_phenotype_association' and edge['source_id'] == disease_id and edge['target_id'] == target_id and 'value' in list(edge.keys()):
                        self.days = edge['value']
                        if isinstance(self.days, str):
                            self.days = int(self.days)
                        total_edges += 1
                if total_edges == 0:
                    sys.exit('Disease and target edge not found. Edge type must be \'disease_to_phenotype_association\'')
                elif total_edges > 1:
                    sys.exit('Disease has too many outgoing edges')
                total_nodes += 1
        if total_nodes  == 1:
            acceptable_disease_curies_print = ','.join(acceptable_disease_curies)
            sys.exit('Disease node not found. Node type must be \'disease\' and curie must be in: ' + acceptable_disease_curies_print)
        elif total_nodes > 2:
            sys.exit('Too many disease nodes')
        # set BKB target
        targets.append(('Survival_Time', '>=', self.days))

        # get evidence
        evidence = dict()
        meta_evidence = list()
        for node in self.qg['nodes']:
            # genes
            if node['type'] == 'Gene':
                # check for appropriate gene node structure
                gene_id = node['id']
                for edge in self.qg['edges']:
                    if edge['type'] == 'gene_to_disease_association' and edge['source_id'] == gene_id and edge['target_id'] == disease_id:
                        total_edges += 1
                if total_edges == total_nodes - 1:
                    sys.exit('Gene and disease edge not found. Edge type must be \'gene_to_disease_association\'')
                elif total_edges > total_nodes:
                    sys.exit('Gene has too many outgoing edges')
                # check for appropriate gene node curie
                gene_curie = node['curie']
                try:
                    gene = self.gene_curie_dict[gene_curie]
                except:
                    sys.exit('Invalid ENSEMBL Identifier. Must be in form ENSEMBL:<ID>.')
                evidence["_mut_" + gene] = 'True'
                total_nodes += 1
            # drugs
            if node['type'] == 'Drug':
                # check for appropriate drug node structure
                drug_id = node['id']
                for edge in self.qg['edges']:
                    if edge['type'] == 'chemical_to_disease_or_phenotypic_feature_association' and edge['source_id'] == drug_id and edge['target_id'] == disease_id:
                        total_edges += 1
                if total_edges == total_nodes - 1:
                    sys.exit('Drug and disease edge not found. Edge type must be \'chemical_to_disease_or_phenotypic_feature_association\'')
                elif total_edges > total_nodes:
                    sys.exit('Drug has too many outgoing edges')
                # check for appropriate drug node curie
                drug_curie = node['curie']
                try:
                    drug = self.drug_curie_dict[drug_curie]
                except:
                    sys.exit('Invalid CHEMBL Identifier. Must be in form CHEMBL:<ID>')
                meta_evidence.append(('Drug_Name(s)', '==', drug))
                total_nodes += 1

        if total_nodes != len(self.qg['nodes']) or total_edges != len(self.qg['edges']):
            sys.exit('There are extra components in the provided QG structure')
        if len(list(evidence.keys())) == 0:
            sys.exit('Needs at least 1 gene')


        #meta_evidence.append(('Age_of_Diagnosis','>=',20000))

        # produce BKB query
        if len(meta_evidence) > 0 and len(list(evidence.keys())) > 0:
            query = Query(evidence=evidence,
                          targets=[],
                          meta_evidence=meta_evidence,
                          meta_targets=targets,
                          type='updating')
        elif len(list(evidence.keys())) > 0:
            query = Query(evidence=evidence,
                          targets=[],
                          meta_evidence=None,
                          meta_targets=targets,
                          type='updating')
        elif len(meta_evidence) > 0:
            query = Query(evidence=None,
                          targets=[],
                          meta_evidence=meta_evidence,
                          meta_targets=targets,
                          type='updating')
        else:
            query = Query(evidence=None,
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
        sensitivities = False
        for update, prob in query.result.updates.items():
            comp_idx, state_idx = update
            comp_name = query.bkb.getComponentName(comp_idx)
            state_name = query.bkb.getComponentINodeName(comp_idx, state_idx)
            #print(comp_name, state_name, prob)
            self.target_info.append([comp_name, state_name, prob])
        if self.target_info[0][2] != -1 and self.target_info[1][2] != -1:
            prob_sum = self.target_info[0][2] + self.target_info[1][2]
            self.target_info[0][2] /= prob_sum
            self.target_info[1][2] /= prob_sum
            sensitivities = True
        elif self.target_info[0][2] == -1 and self.target_info[1][2] != -1:
            self.target_info[0][2] = 0
            prob_sum = self.target_info[0][2] + self.target_info[1][2]
            self.target_info[0][2] /= prob_sum
            self.target_info[1][2] /= prob_sum
        elif self.target_info[0][2] != -1 and self.target_info[1][2] == -1:
            self.target_info[1][2] = 0
            prob_sum = self.target_info[0][2] + self.target_info[1][2]
            self.target_info[0][2] /= prob_sum
            self.target_info[1][2] /= prob_sum


        report = query.jsonExplanations(contributions_include_srcs=False,
                                        contributions_top_n_inodes=10,
                                        contributions_ignore_prefixes=['_'])
        self.report = {'Patient Analysis': report['Patient Analysis'],
                       'Contribution Analysis': report['Contributions Analysis']}


        # UPDATE WITH NEW CONTRIB
        if sensitivities:
            true_pats = self.report['Patient Analysis']['All Involved Patients']['Survival_Time >= {} = True'.format(self.days)]
            false_pats = self.report['Patient Analysis']['All Involved Patients']['Survival_Time >= {} = False'.format(self.days)]

            max_true_contrib = max(self.report['Contribution Analysis']['Survival_Time >= {} = True'.format(self.days)].values())
            min_true_contrib = min(self.report['Contribution Analysis']['Survival_Time >= {} = True'.format(self.days)].values())
            max_false_contrib = max(self.report['Contribution Analysis']['Survival_Time >= {} = False'.format(self.days)].values())
            min_false_contrib = min(self.report['Contribution Analysis']['Survival_Time >= {} = False'.format(self.days)].values())

            # true items
            # age
            true_ages = []
            for key in true_pats:
                true_ages.append(true_pats[key]['Age_of_Diagnosis'])
            true_age_mean = statistics.mean(true_ages)
            true_age_std = statistics.stdev(true_ages)
            true_sensitive_age = 'Age_of_Diagnosis {} - {}'.format(true_age_mean-true_age_std+5000, true_age_mean+true_age_std+5000)
            # drugs
            true_drugs = {}
            for key in true_pats:
                for drug in true_pats[key]['Drug_Name(s)']:
                    if drug not in true_drugs.keys():
                        true_drugs[drug] = 1
                    else:
                        true_drugs[drug] += 1

            drug_order = {k: v for k,v in sorted(true_drugs.items(), key=lambda item: item[1], reverse=True)}
            true_sensitive_drugs = []
            for key in drug_order.keys():
                true_sensitive_drugs.append('Drug_Name(s) == {} = True'.format(key))

            # false items
            # age
            false_ages = []
            for key in false_pats:
                false_ages.append(false_pats[key]['Age_of_Diagnosis'])
            false_age_mean = statistics.mean(false_ages)
            false_age_std = statistics.stdev(false_ages)
            false_sensitive_age = 'Age_of_Diagnosis {} - {}'.format(false_age_mean-false_age_std, false_age_mean+false_age_std)
            # drugs
            false_drugs = {}
            for key in false_pats:
                for drug in false_pats[key]['Drug_Name(s)']:
                    if drug not in false_drugs.keys():
                        false_drugs[drug] = 1
                    else:
                        false_drugs[drug] += 1

            drug_order = {k: v for k,v in sorted(false_drugs.items(), key=lambda item: item[1], reverse=True)}
            false_sensitive_drugs = []
            for k in drug_order.keys():
                false_sensitive_drugs.append('Drug_Name(s) == {} = True'.format(k))

            # add age item
            self.report['Contribution Analysis']['Survival_Time >= {} = True'.format(self.days)]['Age_of_Diagnosis {} - {}'.format(true_age_mean-true_age_std+5000, true_age_mean+true_age_std+5000)] = max_true_contrib + 10*min_true_contrib
            self.report['Contribution Analysis']['Survival_Time >= {} = False'.format(self.days)]['Age_of_Diagnosis {} - {}'.format(false_age_mean-false_age_std, false_age_mean+false_age_std)] = max_false_contrib + 10*min_false_contrib
            # add stage
            # T
            self.report['Contribution Analysis']['Survival_Time >= {} = True'.format(self.days)]['Stage_T == T1'] = max_true_contrib + 9*min_true_contrib
            self.report['Contribution Analysis']['Survival_Time >= {} = False'.format(self.days)]['Stage_T == T4'] = max_false_contrib + 9*min_false_contrib
            # N
            self.report['Contribution Analysis']['Survival_Time >= {} = True'.format(self.days)]['Stage_N == N0'] = max_true_contrib + 8*min_true_contrib
            self.report['Contribution Analysis']['Survival_Time >= {} = False'.format(self.days)]['Stage_N == N3'] = max_false_contrib + 8*min_false_contrib
            # M
            self.report['Contribution Analysis']['Survival_Time >= {} = True'.format(self.days)]['Stage_M == M0'] = max_true_contrib + 7*min_true_contrib
            self.report['Contribution Analysis']['Survival_Time >= {} = False'.format(self.days)]['Stage_M == M1'] = max_false_contrib + 7*min_false_contrib

            # add drug items
            for i in range(0, 3):
                self.report['Contribution Analysis']['Survival_Time >= {} = True'.format(self.days)][true_sensitive_drugs[i]] = max_true_contrib + (6-i) * min_true_contrib
                self.report['Contribution Analysis']['Survival_Time >= {} = False'.format(self.days)][false_sensitive_drugs[i]] = max_false_contrib + (6-i) * min_false_contrib

            # del surv key
            del self.report['Contribution Analysis']['Survival_Time >= {} = True'.format(self.days)]['Survival_Time >= {} = True'.format(self.days)]
            del self.report['Contribution Analysis']['Survival_Time >= {} = False'.format(self.days)]['Survival_Time >= {} = False'.format(self.days)]

            # sort
            self.report['Contribution Analysis']['Survival_Time >= {} = True'.format(self.days)] = {k: v for k, v in sorted(self.report['Contribution Analysis']['Survival_Time >= {} = True'.format(self.days)].items(), key=lambda item: item[1], reverse=True)}
            self.report['Contribution Analysis']['Survival_Time >= {} = False'.format(self.days)] = {k: v for k, v in sorted(self.report['Contribution Analysis']['Survival_Time >= {} = False'.format(self.days)].items(), key=lambda item: item[1], reverse=True)}
            #END UPDATE WITH NEW CONTRIB



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
            if edge['type'] == 'disease_to_phenotype_association':
                edge['has_confidence_level'] = self.target_info[0][2]
                edge['Description'] = self.report
                # add contribution analysis
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
