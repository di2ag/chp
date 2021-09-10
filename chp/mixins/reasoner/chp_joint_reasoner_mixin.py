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

    def run_query(self, query, bkb_type='gene'):
        # Compute joint probability
        evidence = query.compose_evidence(with_dynamic=False, meta_tag=False)
        res, contrib = self.joint_reasoner.compute_joint(
            evidence,
            query.targets,
            continuous_evidence=query.dynamic_evidence,
            continuous_targets=query.dynamic_targets,
            interp_type=bkb_type
        )
        # Set query parameters
        query.result = res
        query.contributions = contrib
        query.from_joint_reasoner = True
        return query
