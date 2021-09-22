import logging
import pickle

from pybkb.python_base.reasoning.reasoning import updating
from pybkb.python_base.reasoning.joint_reasoner import JointReasoner

logger = logging.getLogger(__name__)
#logger.setLevel(logging.INFO)

class ChpJointReasonerMixin:
    def _setup_reasoner(self):
        if self.bkb_handler.gene_interpolations_path is not None:
            with open(self.bkb_handler.gene_interpolations_path, 'rb') as f:
                gene_interpolations = pickle.load(f)
        else:
            gene_interpolations_path = None

        if self.bkb_handler.drug_interpolations_path is not None:
            with open(self.bkb_handler.drug_interpolations_path, 'rb') as f:
                drug_interpolations = pickle.load(f)
        else:
            drug_interpolations_path = None

        self.joint_reasoner = JointReasoner(self.patient_data, 
                                            gene_interpolations=gene_interpolations,
                                            drug_interpolations=drug_interpolations)
        logger.info('Setup Joint Reasoner.')

    def run_query(self, query, interpolation_type=None, contribution_type=None):
        # Compute joint probability
        '''
        if bkb_type == 'gene':
            if return_contributions:
                contribution_feature_type = 'ENSEMBL'
        elif bkb_type == 'drug':
            if return_contributions:
                contribution_feature_type = 'CHEMBL.COMPOUND'
        '''
        if contribution_type == 'gene':
            contribution_feature_type = 'ENSEMBL'
        elif contribution_type == 'drug':
            contribution_feature_type = 'CHEMBL.COMPOUND'
        else:
            contribution_feature_type = None

        evidence = query.compose_evidence(with_dynamic=False, meta_tag=False)
        res, contrib = self.joint_reasoner.compute_joint(
            evidence,
            query.targets,
            continuous_evidence=query.dynamic_evidence,
            continuous_targets=query.dynamic_targets,
            contribution_features=contribution_feature_type,
            interpolation_type=interpolation_type
        )
        # Set query parameters
        query.result = res
        query.contributions = contrib
        query.from_joint_reasoner = True
        return query
