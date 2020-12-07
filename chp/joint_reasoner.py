"""
    Source code developed by DI2AG.
    Thayer School of Engineering at Dartmouth College
    Authors:    Dr. Eugene Santos, Jr
                Mr. Chase Yakaboski,
                Mr. Gregory Hyde,
                Mr. Luke Veenhuis,
                Dr. Keum Joo Kim
"""
import time
import tqdm
import os
import itertools
import math
import random
import pickle
import logging
from collections import defaultdict
import json
import numpy as np
import logging
import sys

#from chp.util import process_operator
#from pybkb.python_base.learning import BkbDataBuilder

from pybkb.python_base.reasoning.joint_reasoner import JointReasoner
from chp_data.bkb_handler import BkbDataHandler

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

def recursive_dd(max_depth, end_type, cur_depth=0):
    if cur_depth < max_depth:
        return defaultdict(recursive_dd(max_depth, end_type, cur_depth=cur_depth+1))
    else:
        return defaultdict(end_type)

class ChpJointReasoner:
    def __init__(
            self,
            bkb_handler=None,
            discretize=None,
            num_genes=None,
            num_drugs=None
            ):
        """ Joint Reasoner is responsible for taking in a patient_data dictionary and
            putting it in a nice format and just calculating joint probabilities explicitly.
            No BKBs required.

            :param patient_data: all patient data
            :type patient_data: dict
            :param bkb_handler: Handler to data file paths. Defaults to None.
            :type patient_data: chp_data.bkb_handler.BkbDataHandler
            :param discretize: Discretize continuous data to this many bins. Defaults to None. 
            :type patient_data: int
            :param num_genes: Only build the bkbs with the most frequent num_genes number of genes.
            :type num_genes: int
            :param num_drugs: Only build the bkbs with the most frequent num_drugs number of drugs.
            :type num_drugs: int
        """
        logger.info('Starting setup.')
        if bkb_handler is None:
            # Load default
            self.bkb_handler = BkbDataHandler()
        else:
            self.bkb_handler = bkb_handler
        # Load patient data
        with open(self.bkb_handler.patient_data_pk_path, 'rb') as f_:
            self.patient_data = pickle.load(f_)
        # Load curies
        with open(self.bkb_handler.curies_path, 'r') as f_:
            self.curies = json.load(f_)
        self.bins = discretize
        self.num_genes = num_genes
        self.num_drugs = num_drugs
        self.patient_data = self._process_patient_data(patient_data, discretize, num_genes, num_drugs)
        self.reasoner = JointReasoner(self.patient_data)
        self._calculate_feature_state_tables()
        self.num_patients = len(self.patient_data)
        logger.info('Completed setup.')

    def analyze_query(self, query):
        targets, evidence = self._process_query(query)

    def _process_query(self, query):
        # Process evidence
        evidence = query.evidence
        for meta_evidence in query.meta_evidence:

        
        targets = []

    def _process_meta_evidence(self, 




    def _calculate_feature_state_tables(self):
        # Used for just some functionalities
        bkb_builder = BkbDataBuilder(self.patient_data)
        self.unique_feature_set_counts = bkb_builder.unique_feature_set_counts
        self.unique_feature_set_examples = bkb_builder.unique_feature_set_examples
        self.features = bkb_builder.features_list
        self.states = bkb_builder.states
        self.feature_map = {idx: feature_name for idx, feature_name in enumerate(self.features)}

    def compute_joint_from_trapi_query_graph(self, query_graph):
        pass

    def compute_joint(self, evidence, targets, contribution_features=None):
        # Convert evidence to list indices for future feature set calculations
        processed_evidence = {}
        for feature, state in evidence.items():
            feature_idx = self.features.index(feature)
            processed_evidence[feature_idx] = state
        logger.info('Processed evidence.')
        processed_targets = []
        for target_feature in targets:
            processed_targets.append(self.features.index(target_feature))
        logger.info('Processed targets.')
        if contribution_features is not None:
            non_evidence_feature_indices = self._get_contribution_indices(
                contribution_features,
                processed_evidence,
                processed_targets,
            )
            logger.info('Processed contribution variables.')
        contributions_table = defaultdict(lambda: defaultdict(int))
        # Process updates
        results = defaultdict(int)
        evidence_inferences = 0
        logger.info('Starting reasoning.')
        start_time = time.time()
        traversed_patients = set()
        for (feature_set, count), (_, patients) in tqdm.tqdm(
            zip(self.unique_feature_set_counts.items(),
                self.unique_feature_set_examples.items()),
            desc='Finding supported inferences',
            total=len(self.unique_feature_set_counts),
            leave=False):
            if self._is_consistent_with_evidence(processed_evidence, feature_set):
                evidence_inferences += self._prob(count)
                for target in processed_targets:
                    for state in self.states[target]:
                        if self._is_consistent_with_target_state(target, state, feature_set):
                            # Don't double count a patient
                            num_patients_already_counted = len(set(patients) - traversed_patients)
                            if num_patients_already_counted != count:
                                count -= num_patients_already_counted
                                logger.info('Already encountered some patients.')
                            results[(self.feature_map[target], state)] += self._prob(count)
                            if contribution_features is not None:
                                contributions_table = self._update_contributions(
                                    contributions_table,
                                    target,
                                    state,
                                    feature_set,
                                    non_evidence_feature_indices,
                                    count,
                                )
        logger.info('Finished reasoning in {} seconds.'.format(time.time() - start_time))
        return dict(results), contributions_table

    def _get_contribution_indices(self, contribution_features, processed_evidence, processed_targets):
        if contribution_features == 'all':
            return list(set([feature_idx for feature_idx in range(len(self.features))]) - set(list(processed_evidence.keys()) + processed_targets))
        else:
            return [self.features.index(feature_name) for feature_name in contribution_features]

    def _prob(self, count):
        return float(count / self.num_patients)

    def _update_contributions(self, contributions_table, target, state, feature_set, non_evidence_feature_indices, count):
        for non_evidence_feature in non_evidence_feature_indices:
            target_name = self.feature_map[target]
            non_evidence_name = self.feature_map[non_evidence_feature]
            contributions_table[(target_name, state)][(non_evidence_name, feature_set[non_evidence_feature])] += self._prob(count)
        return contributions_table

    def _is_consistent_with_evidence(self, evidence, feature_set):
        for feature, state in evidence.items():
            if feature_set[feature] != state:
                return False
        return True

    def _is_consistent_with_target_state(self, target, state, feature_set):
        if feature_set[target] == state:
            return True
        return False

    def _process_feature_dict(self, feature_dict):
        processed_feature_dict = {}
        # Process phenotypic info
        for feature in feature_dict:
            for _curie, names in self.curies["phenotypic_feature"].items():
                if feature in names:
                    processed_feature_dict[_curie] = feature_dict[feature]
        # Process drugs
        for _curie in self.all_drug_counts:
            if _curie in feature_dict["drug_curies"]:
                processed_feature_dict[_curie] = 'True'
            else:
                processed_feature_dict[_curie] = 'False'
        # Process genes
        for _curie in self.all_gene_counts:
            if _curie in feature_dict["gene_curies"]:
                processed_feature_dict[_curie] = 'True'
            else:
                processed_feature_dict[_curie] = 'Unknown'
        return processed_feature_dict

    def _discretize(self, patient_dict, bins):
        continuous_feature_values = defaultdict(list)
        # Gather all continous variables
        for patient, feature_dict in patient_dict.items():
            for feature, state in feature_dict.items():
                if type(state) == int or type(state) == float:
                    #TODO: I am just going to convert to years here cause its easy this will need to be changed.
                    continuous_feature_values[feature].append(state / 365)
        # Get evenly spaced levels based on max and min values
        levels = {}
        for feature, values in continuous_feature_values.items():
            levels[feature] = list(np.round(np.linspace(min(values), max(values), num=bins-1), decimals=1))
        # Go back through patients and bin continuous values
        #print(levels)
        for patient, feature_dict in patient_dict.items():
            for feature in continuous_feature_values:
                #print(feature)
                #print(feature_dict[feature] / 365)
                _bin = np.digitize(feature_dict[feature] / 365, levels[feature])
                #print(_bin)
                if _bin == 0:
                    bin_name = '<{}'.format(levels[feature][0])
                elif _bin == bins-1:
                    bin_name = '>={}'.format(levels[feature][bins-1])
                else:
                    bin_name = '({}, {})'.format(levels[feature][_bin - 1], levels[feature][_bin])
                feature_dict[feature] = bin_name
                #print(bin_name)
                #input()
        return patient_dict, levels

    def _process_patient_data(self, patient_data, discretize, num_genes, num_drugs):
        """ Put's the patient data in the right dictionary format for pybkb's learning module.
        """
        processed_patient_dict = defaultdict(dict)
        # Collect all genes and drugs
        all_gene_counts = defaultdict(int)
        all_drug_counts = defaultdict(int)
        for patient, unprocessed_feature_dict in patient_data.items():
            for _gene in unprocessed_feature_dict["gene_curies"]:
                all_gene_counts[_gene] += 1
            for _drug in unprocessed_feature_dict["drug_curies"]:
                all_drug_counts[_drug] += 1
        # Only process num_genes and num_drugs
        if num_genes is not None:
            top_gene_counts = sorted([(_count, _gene) for _gene, _count in all_gene_counts.items()])[-num_genes:]
            self.all_gene_counts = {_gene: _count for _count, _gene in top_gene_counts}
        else:
            self.all_gene_counts = all_gene_counts
        if num_drugs is not None:
            top_drug_counts = sorted([(_count, _drug) for _drug, _count in all_drug_counts.items()])[-num_drugs:]
            self.all_drug_counts = {_drug: _count for _count, _drug in top_drug_counts}
        else:
            self.all_drug_counts = all_drug_counts
        # Process each feature_dict
        for patient, unprocessed_feature_dict in patient_data.items():
            processed_feature_dict = self._process_feature_dict(unprocessed_feature_dict)
            processed_patient_dict[str(patient)] = processed_feature_dict
        # Discretize continuous data
        if discretize is not None:
            processed_patient_dict, self.levels= self._discretize(processed_patient_dict, discretize)
        return processed_patient_dict


