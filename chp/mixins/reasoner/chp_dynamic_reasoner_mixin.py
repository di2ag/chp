import compress_pickle
import logging
import copy
import time

from pybkb.python_base.reasoning.reasoning import updating
from pybkb.python_base.learning.bkb_builder import LinkerBuilder

logger = logging.getLogger(__name__)

class ChpDynamicReasonerMixin:
    """ This reasoner class is responsible for running all dynamic BKB queries.
    """
    def _setup_reasoner(self):
        # Construct linker
        self.linker_builder = LinkerBuilder(self.patient_data)
        logger.info('Constructed Linker Builder from processed patient data.')
        # Load in gene prelinked bkb for bkb_data_handler or appropriate override
        if self.gene_prelinked_bkb_override is None:
            with open(self.bkb_handler.collapsed_gene_bkb_path, 'rb') as f_:
                self.gene_prelinked_bkb = compress_pickle.load(f_, compression='lz4')
            logger.info('Loaded in gene prelinked bkb from: {}'.format(self.bkb_handler.collapsed_gene_bkb_path))
        else:
            self.gene_prelinked_bkb = self.gene_prelinked_bkb_override
            logger.info('Loaded override gene prelinked bkb.')
        # Load in drug prelinked bkb for bkb_data_handler or appropriate override
        if self.drug_prelinked_bkb_override is None:
            with open(self.bkb_handler.collapsed_drug_bkb_path, 'rb') as f_:
                self.drug_prelinked_bkb = compress_pickle.load(f_, compression='lz4')
            logger.info('Loaded in drug prelinked bkb from: {}'.format(self.bkb_handler.collapsed_drug_bkb_path))
        else:
            self.drug_prelinked_bkb = self.drug_prelinked_bkb_override
            logger.info('Loaded override drug prelinked bkb.')

    def _pool_properties(self, query, bkb):
        """ Pools all the evidence and targets that aren't in the prelinked BKB and gets them
        ready to link.

        Args:
            :param query: The chp query that is to be run.
            :type query: chp.query.Query
            :param bkb: The prelinked bkb that is being used.
            :type bkb: pybkb.common.bayesianKnowledgeBase.BayesianKnowledgeBase

        Returns:
            :return: Tuple containing the feature properties dictionary that will be passed to the 
            PyBKB Linker module and the feature that should not be passed to the linker module.
            :rtype: tuple
        """
        features_not_to_format = []
        feature_properties = {}
        # Add dynamic evidence and targets
        feature_properties.update(query.dynamic_evidence)
        feature_properties.update(query.dynamic_targets)
        # Check normal evidence
        for feature, state in query.evidence.items():
            if not self.linker_builder.is_feature_in_bkb(feature, bkb):
                raise ValueError('Normal evidence should exist in the BKB.')
            else:
                features_not_to_format.append(feature)
        # Check meta evidence
        for feature, state in query.meta_evidence.items():
            meta_feature = '_' + feature
            if not self.linker_builder.is_feature_in_bkb(meta_feature, bkb):
                logger.info('Could not find interpolate feature: {} in bkb.'.format(meta_feature))
                features_not_to_format.append(feature)
            else:
                features_not_to_format.append(feature)
        return feature_properties, features_not_to_format
    
    def _check_evidence(self, evidence, bkb):
        """ Ensures all specified evidence, i.e. random variables and their respective states are actually in the linked BKB.
        """
        for feature, state in evidence.items():
            if not self.linker_builder.is_feature_in_bkb(feature, bkb):
                logger.info('Could not find feature: {} in bkb so we are removing it from the evidence'.format(feature))
                return False
            elif not self.linker_builder.is_feature_state_in_bkb(feature, state, bkb):
                logger.info('Could not find state: {} of feature: {} in bkb so we are removing it from the evidence'.format(state, feature))
                return False
        return True
    
    def _check_targets(self, targets, bkb):
        """ Ensures all specified target random variables are actually in the linked BKB.
        """
        for feature in targets:
            if not self.linker_builder.is_feature_in_bkb(feature, bkb):
                logger.info('Could not find feature: {} in bkb so we are removing it from the targets'.format(feature))
                return False
        return True

    def run_query(self, query, bkb_type='gene'):
        """ Method to calculate an update for a CHP Query.

        Args:
            :param query: The CHP query to run a update.
            :type query: chp.query.Query
            :param bkb_type: Used to specify the prelinked bkb to use. As of right now this can
            either by 'drug' or 'gene'.
            :type bkb_type: str

        :return: Augemented CHP Query with all the result attributes filled in according to the BKB update.
        :rtype: chp.query.Query
        """
        # Need to make a copy of prelinked bkb
        if bkb_type == 'gene':
            bkb = copy.deepcopy(self.gene_prelinked_bkb)
        elif bkb_type == 'drug':
            bkb = copy.deepcopy(self.drug_prelinked_bkb)
        else:
            raise ValueError('Unrecognized bkb type: {}'.format(bkb_type))
        # Pool any dynamic evidence and/or targets for linking
        feature_properties, features_not_to_format = self._pool_properties(query, bkb)
        # Link BKB based on dynamic evidence in query
        linked_bkb = self.linker_builder.link(feature_properties, bkb)
        # Compose evidence and targets
        evidence = query.compose_evidence()
        targets = query.compose_targets()
        # Run checks
        if not self._check_evidence(evidence, linked_bkb):
            logger.critical('Evidence check failed and pieces of evidence where removed. Check Log!')
            raise ValueError('Evidence failed. Check logs.')
        if not self._check_targets(targets, linked_bkb):
            logger.critical('Targets check failed and targets where removed. Check Log!')
            raise ValueError('Targets failed. Check logs.')
        # Run Updating
        start_time = time.time()
        res = updating(bkb,
                       evidence,
                       targets,
                       hosts_filename=self.hosts_filename,
                       num_processes_per_host=self.num_processes_per_host,
                       venv=self.venv,
                      )
        compute_time = time.time() - start_time
        logger.info('Ran update in {} seconds.'.format(compute_time))
        # Update query with results
        query.result = res
        query.compute_time = compute_time
        return query

""" Possibly deprecate on next version.
    
    def _format_evidence(self, query, features_not_to_format):
        new_format = {}
        # Process standard evidence
        for feature, state in query.evidence.items():
            if feature in features_not_to_format:
                new_format[feature] = '{}'.format(state)
            else:
                new_format[feature] = '== {}'.format(state)
        # Process dynamic evidence
        new_format.update(self._dynamic_to_standard_format(query.dynamic_evidence))
        return new_format

    def _dynamic_to_standard_format(self, dynamic_evidence):
        standard = {}
        if dynamic_evidence is not None:
            for feature, prop in dynamic_evidence.items():
                standard[feature] = '{} {}'.format(prop["op"], prop["value"])
        return standard
"""
