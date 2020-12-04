import json
import unittest
import pickle
import logging

from chp_data.bkb_handler import BkbDataHandler
from chp.joint_reasoner import JointReasoner

logging.basicConfig(level=logging.INFO)

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

class TestJointReasoner(unittest.TestCase):

    def setUp(self):
        self.bkb_handler = BkbDataHandler()
        with open(self.bkb_handler.patient_data_pk_path, 'rb') as f_:
            self.patient_data = pickle.load(f_)
        self.reasoner = JointReasoner(
            self.patient_data,
            self.bkb_handler,
            discretize=10,
            num_genes=None,
            num_drugs=None,
        )

    def test_joint_reason(self):
        #print(self.reasoner.features)
        #print(self.reasoner.states)

        # Specify evidence
        evidence = {'ENSEMBL:ENSG00000155657': 'True'}
        # Specify targets
        targets = ['EFO:0000714']
        # Specify contribution features
        contribution_features = [
            'EFO:0008007',
            'PATO:0000047',
            'EFO:0000714',
            'CHEMBL:CHEMBL83',
            'CHEMBL:CHEMBL53463',
            'CHEMBL:CHEMBL88'
        ]
        contribution_features = 'all'
        res, contrib =  self.reasoner.compute_joint(evidence, targets, contribution_features)
        print(res)
        #print(contrib)

