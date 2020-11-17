"""
    Source code developed by DI2AG.
    Thayer School of Engineering at Dartmouth College
    Authors:    Dr. Eugene Santos, Jr
                Mr. Chase Yakaboski,
                Mr. Gregory Hyde,
                Mr. Luke Veenhuis,
                Dr. Keum Joo Kim
"""

import copy
import csv
import sys
import pickle

from chp.query import Query
from chp.reasoner import Reasoner
from chp_data.bkb_handler import BkbDataHandler

class OneHopHandler:

    """ OneHopeHandler is the handler for 1-hop queries. That is
        query graphs (QGs) that consists of 2 nodes and a single edge.

        :param query: the query graph sent by the ARA.
        :type query: dict
        :param hosts_filename: a filename for a stored QG. Defaults to None
        :type hosts_filename: str
        :param num_processes_per_host: Not implemented thouroughly, but would be
            used for distributed reasoning.
        :type num_processes_per_host: int
        :param max_results: specific to 1-hop queries, specifies the number of
            wildcard genes to return.
        :type max_results: int
    """

    def __init__(self,
                 query,
                 hosts_filename=None,
                 num_processes_per_host=0,
                 max_results=50):
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
            self.results = []
        else:
            self.results = self.query['results']

        self.op = '>='
        self.value = 970

        # Instiatate Reasoner
        self.bkb_data_handler = BkbDataHandler(dataset_version='1.4')
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
                self.true_gene_contrib[row[1]] = 0
                self.false_gene_contrib[row[1]] = 0
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

    def checkQuery(self):
        """ Currently not implemented. Would check validity of query.
        """
        return True

    def buildQueries(self):
        """ Parses over the sent query graph to form a BKB query.

            :return: A  internal CHP query.
            :rtype: Query
        """

        evidence = dict()
        targets = list()
        self.contribution_target = None

        if len(self.qg['nodes']) > 2 or len(self.qg['edges']) > 1:
            sys.exit('1 hop quries can only have 2 nodes and 1 edge')

        # check edge for source and target
        edge_key = list(self.qg['edges'].keys())[0]
        edge = self.qg['edges'][edge_key]
        if 'subject' not in edge.keys() or 'object' not in edge.keys():
            sys.exit('Edge must have both a \'subject\' and and \'object\' key')
        subject = edge['subject']
        obj = edge['object']

        # get drug node
        if self.qg['nodes'][subject]['category'] != 'biolink:Drug':
            sys.exit('Subject node must be \'category\' biolink:Drug')
        elif 'id' not in self.qg['nodes'][subject].keys():
            sys.exit('Must have \'id\' key in drug node')
        self.drug_curie = self.qg['nodes'][subject]['id']
        if self.drug_curie not in self.drug_curie_dict.keys():
            sys.exit('Invalid CHEMBL Identifier. Must be CHEMBL:<ID>')
        evidence['demo_{}'.format(self.drug_curie)] = 'True'

        # ensure gene is wildcard
        if self.qg['nodes'][obj]['category'] != 'biolink:Gene':
            sys.exit('Object node must be \'category\' biolink:Gene')
        elif 'id' in self.qg['nodes'][obj].keys():
            sys.exit('Must NOT have \'id\' key in gene node')

        # default survival time
        targets.append(('survival_time', self.op, self.value))

        query = Query(evidence=evidence,
                      targets=[],
                      meta_evidence=None,
                      meta_targets=targets,
                      type='updating')
        self.chp_query = query
        return query

    def runQueries(self):

        """ Runs build BKB query to calculate probability of survival.
            Uniquely for the 1-hop query the probability is not going to
            be returned, rather the relative gene contributions are.
            Contributions for each gene are calculuated and classified under
            their true/false target assignments.
        """

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
                                        contributions_top_n_inodes=self.max_results,
                                        contributions_ignore_prefixes=['_'])
        self.report = {'patient_analysis': report['Patient Analysis'],
                       'contribution_analysis': report['Contributions Analysis']}

        if 'survival_time {} {} = True'.format(self.op, self.value) in self.report['contribution_analysis'].keys():
            true_contrib = self.report['contribution_analysis']['survival_time {} {} = True'.format(self.op, self.value)]['demo_{} = True'.format(self.drug_curie)]
            true_pats = len(self.report['patient_analysis']['All Involved Patients']['survival_time {} {} = True'.format(self.op, self.value)].keys())
            true_ind_cont = float(true_contrib)/float(true_pats)/self.truth_assignment
        else:
            true_contrib = 0
            true_pats = 0
            true_ind_cont = 0

        if 'survival_time {} {} = False'.format(self.op, self.value) in self.report['contribution_analysis'].keys():
            false_contrib = self.report['contribution_analysis']['survival_time {} {} = False'.format(self.op, self.value)]['demo_{} = True'.format(self.drug_curie)]
            false_pats = len(self.report['patient_analysis']['All Involved Patients']['survival_time {} {} = False'.format(self.op, self.value)].keys())
            false_ind_cont = float(false_contrib)/float(false_pats)/self.false_assignment
        else:
            false_contrib = 0
            false_pats = 0
            false_ind_cont = 0

        patient_dict = pickle.load(open(self.bkb_data_handler.patient_data_pk_path, 'rb'))
        for key in patient_dict:
            pat = patient_dict[key]
            if pat['survival_time'] > self.value and self.drug_curie in pat['drug_curies']:
                for gene in pat['gene_curies']:
                    self.true_gene_contrib[gene] += true_ind_cont
            elif pat['survival_time'] < self.value and self.drug_curie in pat['drug_curies']:
                for gene in pat['gene_curies']:
                    self.false_gene_contrib[gene] += false_ind_cont
        self.report['contribution_analysis']['survival_time {} {} = True'.format(self.op, self.value)] = {'{}-{}'.format(k,self.gene_curie_dict[k]): v for k,v in sorted(self.true_gene_contrib.items(), key=lambda item: item[1], reverse=True)}
        self.report['contribution_analysis']['survival_time {} {} = False'.format(self.op, self.value)] = {'{}-{}'.format(k,self.gene_curie_dict[k]): v for k,v in sorted(self.false_gene_contrib.items(), key=lambda item: item[1], reverse=True)}

    def constructDecoratedKG(self):

        """ Knowledge Graph (KG) is a copy of the Query graph. However we replace
            the wildcard gene node with all of the relative contributing genes.
            Edges are also updated to account for these added genes. Results map
            New gene nodes to the QG wildcard gene.

            :return: reasoner_std is our API response message combining KG and results.
            :rtype: dict
        """

        # Calculate Relative Contributions and sort
        rel_contrib_dict = {}
        true_contribs = self.report['contribution_analysis']['survival_time {} {} = True'.format(self.op, self.value)]
        false_contribs = self.report['contribution_analysis']['survival_time {} {} = False'.format(self.op, self.value)]
        contrib_keys = list(true_contribs.keys())
        for contrib_key in contrib_keys:
            rel_contrib_dict[contrib_key] = true_contribs[contrib_key] - false_contribs[contrib_key]
        rel_contrib = [(contrib,gene) for gene, contrib in sorted(rel_contrib_dict.items(), key=lambda x: abs(x[1]), reverse=True)]

        self.kg = copy.deepcopy(self.qg)

        self.edge_bindings = dict()
        self.node_bindings = dict()

        # get edge subject, object, edge label and pop edge
        edge_key = list(self.kg['edges'].keys())[0]
        edge = self.kg['edges'][edge_key]
        edge_label = edge['predicate']
        subject = edge['subject']
        obj = edge['object']
        self.kg['edges'].pop(edge_key)

        # move curie to key
        drug_curie = self.kg['nodes'][subject].pop('id')
        self.kg['nodes'][drug_curie] = self.kg['nodes'].pop(subject)
        self.kg['nodes'][drug_curie]['name'] = self.drug_curie_dict[drug_curie]
        self.node_bindings[subject] = drug_curie

        # remove wildcard gene node from kg
        self.kg['nodes'].pop(obj)

        # add kg gene nodes and edges
        edge_count = 0
        node_count = 1
        if len(rel_contrib) < self.max_results:
            self.max_results = len(rel_contrib)
        for contrib, gene in rel_contrib[:self.max_results]:
            # add node
            if 'DEFINE' in gene:
                continue
            gene_curie, gene_name = gene.split('-')
            self.kg['nodes'][gene_curie] = { 'name' : gene_name,
                                             'category' : 'biolink:Gene'}
            # add edge
            self.kg['edges']['kge{}'.format(edge_count)] = {'predicate' : 'biolink:ChemicalToGeneAssociation',
                                                            'subject' : drug_curie,
                                                            'object' : gene_curie,
                                                            'value' : contrib}
            # add to results
            node_binding = {subject : [{'id': drug_curie}],
                            'n{}'.format(node_count) : [{'id':gene_curie}]}
            edge_binding = {edge_key : [{'id':'kge{}'.format(edge_count)}]}
            self.results.append({'node_bindings': node_binding,
                            'edge_bindings': edge_binding})

            edge_count += 1
            node_count += 1

        # query response
        reasoner_std = {'query_graph': self.qg,
                        'knowledge_graph': self.kg,
                        'results': self.results}
        reasoner_std = {'message' : reasoner_std}
        return reasoner_std
