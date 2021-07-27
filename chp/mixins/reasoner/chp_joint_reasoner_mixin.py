import logging

from pybkb.python_base.reasoning.reasoning import updating
from pybkb.python_base.reasoning.joint_reasoner import JointReasoner

logger = logging.getLogger(__name__)
#logger.setLevel(logging.INFO)

class ChpJointReasonerMixin:
    def _setup_reasoner(self):
        self.joint_reasoner = JointReasoner(self.patient_data)
        logger.info('Setup Joint Reasoner.')

    def run_query(self, query):
        # Compute joint probability
        evidence = query.compose_evidence(with_dynamic=False, meta_tag=False)
        res, contrib = self.joint_reasoner.compute_joint(
            evidence,
            query.targets,
            continuous_evidence=query.dynamic_evidence,
            continuous_targets=query.dynamic_targets,
        )
        # Set query parameters
        query.result = res
        query.contributions = contrib
        query.from_joint_reasoner = True
        return query
