import compress_pickle
import logging
import copy
import time

from pybkb.python_base.reasoning.reasoning import updating
from pybkb.python_base.learning.bkb_builder import LinkerBuilder

logger = logging.getLogger(__name__)

class ChpDynamicReasonerMixin:
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
        features_not_to_format = []
        if query.dynamic_evidence is None:
            feature_properties = {}
        else:
            feature_properties = copy.deepcopy(query.dynamic_evidence)
        if query.dynamic_targets is not None:
            feature_properties.update(query.dynamic_targets)
        for feature, state in query.evidence.items():
            if not self.linker_builder.is_feature_in_bkb(feature, bkb):
                if feature[0] == '_':
                    logger.info('Could not find interpolate feature {} in bkb.'.format(feature))
                    features_not_to_format.append(feature)
                    continue
                feature_properties[feature] = {
                    "op": '==',
                    "value": state
                }
            else:
                features_not_to_format.append(feature)
        return feature_properties, features_not_to_format

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

    def _check_evidence(self, query, bkb):
        evidence = {}
        for feature, state in query.evidence.items():
            if not self.linker_builder.is_feature_in_bkb(feature, bkb):
                if feature[0] == '_':
                    logger.info('Could not find interpolated feature: {} in bkb so we are removing it from evidence.'.format(feature))
                else:
                    evidence[feature] = state
            else:
                evidence[feature] = state
        query.evidence = evidence
        return query

    def run_query(self, query, bkb_type='gene'):
        # Need to make a copy of prelinked bkb
        if bkb_type == 'gene':
            bkb = copy.deepcopy(self.gene_prelinked_bkb)
        elif bkb_type == 'drug':
            bkb = copy.deepcopy(self.drug_prelinked_bkb)
        else:
            raise ValueError('Unrecognized bkb type: {}'.format(bkb_type))
        query = self._check_evidence(query, bkb)
        # Pool any dynamic evidence and/or targets for linking
        feature_properties, features_not_to_format = self._pool_properties(query, bkb)
        # Link BKB based on dynamic evidence in query
        linked_bkb = self.linker_builder.link(feature_properties, bkb)
        # Update reasoning evidence with the dynamic evidence
        evidence = self._format_evidence(query, features_not_to_format)
        # Update reasoning targets with the dynamic targets
        if query.targets is None:
            targets = []
        else:
            targets = copy.copy(query.targets)
        if query.dynamic_targets is not None:
            targets += [target_feature for target_feature in query.dynamic_targets]
        # Run update
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
        query.bkb = linked_bkb
        return query
