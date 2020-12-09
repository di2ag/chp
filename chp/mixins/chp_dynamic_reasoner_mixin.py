import compress_pickle
import logging
import copy
import time

from pybkb.python_base.reasoning.reasoning import updating
from pybkb.python_base.learning.bkb_builder import LinkerBuilder

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

class ChpDynamicReasonerMixin:
    def _setup_reasoner(self):
        # Construct linker
        self.linker_builder = LinkerBuilder(self.patient_data)
        logger.info('Constructed Linker Builder from processed patient data.')
        # Load in prelinked bkb for bkb_data_handler
        with open(self.bkb_handler.collapsed_bkb_path, 'rb') as f_:
            self.prelinked_bkb = compress_pickle.load(f_, compression='lz4')
        logger.info('Loaded in prelinked bkb from: {}'.format(self.bkb_handler.collapsed_bkb_path))

    def _pool_properties(self, query, bkb):
        features_not_to_format = []
        if query.dynamic_evidence is None:
            feature_properties = {}
        else:
            feature_properties = copy.deepcopy(query.dynamic_evidence)
        feature_properties.update(query.dynamic_targets)
        for feature, state in query.evidence.items():
            if not self.linker_builder.is_feature_in_bkb(feature, bkb):
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

    def run_query(self, query):
        # Need to make a copy of prelinked bkb
        bkb = copy.deepcopy(self.prelinked_bkb)
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
        return query
