import logging

from pybkb.python_base.reasoning.reasoning import updating
from pybkb.python_base.reasoning.joint_reasoner import JointReasoner

logger = logging.getLogger(__name__)
#logger.setLevel(logging.INFO)

class ChpJointReasonerMixin:
    def _setup_reasoner(self):
        self.joint_reasoner = JointReasoner(self.patient_data)
        logger.info('Setup Joint Reasoner.')

    def _process_evidence(self, evidence):
        """ Since no interpolation is going on remove the '_' from the gene evidence if
            it is passed.
        """
        new_format = {}
        for feature, state in evidence.items():
            if feature[0] == '_':
                new_format[feature[1:]] = state
            else:
                new_format[feature] = state
        return new_format

    def run_query(self, query):
        # Compute joint probability
        res, contrib = self.joint_reasoner.compute_joint(
            self._process_evidence(query.evidence),
            query.targets,
            continuous_evidence=query.dynamic_evidence,
            continuous_targets=query.dynamic_targets,
        )
        # Set query parameters
        query.result = res
        query.contributions = contrib
        query.from_joint_reasoner = True
        return query
