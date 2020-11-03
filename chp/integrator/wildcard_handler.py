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
from collections import defaultdict

from chp.query import Query
from chp.reasoner import Reasoner

from chp_data.bkb_handler import BkbDataHandler

class WildCardHandler:
    def __init__(self,
                 query,
                 hosts_filename=None,
                 num_processes_per_host=0,
                 max_results=10):
        # query graph components
        self.query = query
        self.max_results = max_results
        self.qg = self.query['query_graph']
        if 'knowledge_graph' not in list(self.query.keys()):
            self.kg = { "edges": {},
                        "nodes": {}
                      }
        else:
            self.kg = self.query['knowledge_graph']
        if 'results' not in list(self.query.keys()):
            self.results = [{ "node_bindings": {},
                             "edge_bindings": {}
                           }]
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
        self.contribution_target = None
        for node_id, node in self.qg['nodes'].items():
            if 'curie' not in node.keys():
                if self.contribution_target is None:
                    self.contribution_target = node['type']
                else:
                    sys.exit('You can only have one contribution target. Make sure to leave only one node with a black curie.')
            else:
                if node['type'] == 'chemical_substance':
                    drug_curie = node['curie']
                    try:
                        self.drug = self.drug_curie_dict[drug_curie]
                    except:
                        sys.exit('Invalid CHEMBL Identifier. Must be in form CHEMBL:<ID>')
                    evidence['drug_{}'.format(self.drug)] = 'True'
                elif node['type'] == 'gene':
                    gene_curie = node['curie']
                    try:
                        self.gene = self.gene_curie_dict[gene_curie]
                    except:
                        sys.exit('Invalid ENSEMBL Identifier. Must be in form ENSEMBL:<ID>')
                    evidence['mut_{}'.format(self.gene)] == 'True'
        for edge_id, edge in self.qg['edges'].items():
            if edge['type'] == 'disease_to_phenotypic_feature_association':
                if 'properties' in edge:
                    self.op = edge['properties']['qualifier']
                    self.value = edge['properties']['value']
                else:
                    self.op = '>='
                    self.value = 970

        targets.append(('Survival_Time', self.op, self.value))

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
        self.report = {'patient_analysis': report['Patient Analysis'],
                       'contribution_analysis': report['Contributions Analysis']}
        # total contrib values
        true_contrib = self.report['contribution_analysis']['Survival_Time {} {} = True'.format(self.op, self.value)]['drug_{} = True'.format(self.drug)]
        false_contrib = self.report['contribution_analysis']['Survival_Time {} {} = False'.format(self.op, self.value)]['drug_{} = True'.format(self.drug)]
        # total patients in contrib cat
        true_pats = len(self.report['patient_analysis']['All Involved Patients']['Survival_Time {} {} = True'.format(self.op, self.value)].keys())
        false_pats = len(self.report['patient_analysis']['All Involved Patients']['Survival_Time {} {} = False'.format(self.op, self.value)].keys())
        # true_individual contrib
        true_ind_cont = float(true_contrib)/float(true_pats)/self.truth_assignment
        false_ind_cont = float(false_contrib)/float(false_pats)/self.false_assignment

        self.gene_list_contrib = dict()


        patient_dict = pickle.load(open(self.bkb_data_handler.patient_data_pk_path, 'rb'))
        for key in patient_dict:
            pat = patient_dict[key]
            if pat['Survival_Time'] > self.value and self.drug in pat['Drug_Name(s)']:
                for gene in pat['Patient_Genes']:
                    self.true_gene_contrib[gene] += true_ind_cont
            elif pat['Survival_Time'] < self.value and self.drug in pat['Drug_Name(s)']:
                for gene in pat['Patient_Genes']:
                    self.false_gene_contrib[gene] += false_ind_cont
        self.report['contribution_analysis']['Survival_Time {} {} = True'.format(self.op, self.value)] = {'{}-{}'.format(k,self.gene_to_curie[k]): v for k,v in sorted(self.true_gene_contrib.items(), key=lambda item: item[1], reverse=True)}
        self.report['contribution_analysis']['Survival_Time {} {} = False'.format(self.op, self.value)] = {'{}-{}'.format(k,self.gene_to_curie[k]): v for k,v in sorted(self.false_gene_contrib.items(), key=lambda item: item[1], reverse=True)}


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
        #-- Return top 100 genes contributors for true target
        #-- Calculate Relative Contributions and sort
        rel_contrib_dict = {}
        for target, gene_contrib in self.report['contribution_analysis'].items():
            for gene, contrib in gene_contrib.items():
                if gene in rel_contrib_dict:
                    if 'True' in target:
                        rel_contrib_dict[gene] += contrib
                    else:
                        rel_contrib_dict[gene] -= contrib
                else:
                    if 'True' in target:
                        rel_contrib_dict[gene] = contrib
                    else:
                        rel_contrib_dict[gene] = -contrib
        rel_contrib = [(contrib,gene) for gene, contrib in sorted(rel_contrib_dict.items(), key=lambda x: abs(x[1]), reverse=True)]

        # Construct first result which is the result of the standard probablistic query. 
        self.kg = copy.deepcopy(self.qg)
        # Process Nodes
        node_pairs = defaultdict(None)
        contrib_qg_id = None
        for node_key in list(self.kg["nodes"].keys())[:]:
            qg_node_curie = self.kg['nodes'][node_key].pop('curie', None)
            if qg_node_curie is not None:
                self.kg['nodes'][qg_node_curie] = self.kg['nodes'].pop(node_key)
                if self.kg['nodes'][qg_node_curie]['type'] == 'gene':
                    self.kg['nodes'][qg_node_curie]['name'] = self.gene_curie_dict[qg_node_curie]
                elif self.kg['nodes'][qg_node_curie]['type'] == 'chemical_substance':
                    self.kg['nodes'][qg_node_curie]['name'] = self.drug_curie_dict[qg_node_curie]
                node_pairs[node_key] = qg_node_curie
            else:
                self.kg["nodes"].pop(node_key)

        # Process Edges
        edge_pairs = dict()
        knowledge_edges = 0
        for edge_key in list(self.kg['edges'].keys())[:]:
            if self.kg['edges'][edge_key]['type'] == 'gene_to_disease_association':
                self.kg['edges'].pop(edge_key)
            else:
                kg_id = 'kge{}'.format(knowledge_edges)
                knowledge_edges += 1
                self.kg['edges'][kg_id] = self.kg['edges'].pop(edge_key)
                self.kg['edges'][kg_id]['source_id'] = node_pairs[self.kg['edges'][kg_id]['source_id']]
                self.kg['edges'][kg_id]['target_id'] = node_pairs[self.kg['edges'][kg_id]['target_id']]
                edge_pairs[edge_key] = kg_id
                if self.kg['edges'][kg_id]['type'] == 'disease_to_phenotypic_feature_association':
                    self.kg['edges'][kg_id]['has_confidence_level'] = self.truth_assignment
                    if 'properties' in self.kg['edges'][kg_id].keys() and 'contributions' in self.kg['edges'][kg_id]['properties'].keys() and self.kg['edges'][kg_id]['properties']['contributions'] == True:
                        self.kg['edges'][kg_id]['Description'] = self.report

        # Put first result of standard prob query of only curie nodes (i.e. no wildcard nodes where used as evidence)
        for edge_pair_key in edge_pairs:
            self.results[0]['edge_bindings'][edge_pair_key] = { 'kg_id': str(edge_pairs[edge_pair_key])}
        for node_pair_key in node_pairs:
            self.results[0]['node_bindings'][node_pair_key] = { 'kg_id': str(node_pairs[node_pair_key])}

        # Build relative contribution results and added associated edges into knowledge graph
        for contrib, gene in rel_contrib[:self.max_results]:
            name, curie = gene.split('-')
            rg = copy.deepcopy(self.qg)
            _node_pairs = {}
            _edge_pairs = {}
            # Process node pairs
            for node_id, node in rg["nodes"].items():
                if node["type"] == self.contribution_target:
                    self.kg["nodes"][curie] = copy.deepcopy(node)
                    self.kg["nodes"][curie].update({"name": name})
                    _node_pairs[node_id] = curie
                else:
                    _node_pairs[node_id] = node_pairs[node_id]
            # Process edge pairs
            for edge_id, edge in rg["edges"].items():
                if self.contribution_target == 'gene' and edge["type"] == 'gene_to_disease_association':
                    knowledge_edges += 1
                    kg_edge_id = 'kge{}'.format(knowledge_edges)
                    self.kg["edges"][kg_edge_id] = copy.deepcopy(edge)
                    self.kg["edges"][kg_edge_id]["source_id"] = _node_pairs[self.kg["edges"][kg_edge_id]["source_id"]]
                    self.kg["edges"][kg_edge_id]["target_id"] = _node_pairs[self.kg["edges"][kg_edge_id]["target_id"]]
                    self.kg["edges"][kg_edge_id]["weight"] = contrib
                    _edge_pairs[edge_id] = kg_edge_id
                else:
                    _edge_pairs[edge_id] = edge_pairs[edge_id]
            # Process node and edge binding results
            _res = {"edge_bindings": {},
                    "node_bindings": {}}
            for edge_pair_key in _edge_pairs:
                _res["edge_bindings"][edge_pair_key] = { "kg_id": str(_edge_pairs[edge_pair_key])}
            for node_pair_key in _node_pairs:
                _res["node_bindings"][node_pair_key] = { "kg_id": str(_node_pairs[node_pair_key])}
            self.results.append(_res)

        # query response
        reasoner_std = {'query_graph': self.qg,
                        'knowledge_graph': self.kg,
                        'results': self.results}
        reasoner_std = {'message' : reasoner_std}
        return reasoner_std
