import json
import unittest
import pickle
import logging

from chp_data.bkb_handler import BkbDataHandler

from chp.reasoner import ChpDynamicReasoner, ChpJointReasoner
from chp.query import Query

logging.basicConfig(level=logging.INFO)

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

class TestJointReasoner(unittest.TestCase):

    def setUp(self):
        self.bkb_handler = BkbDataHandler()
        self.joint_reasoner = ChpJointReasoner(self.bkb_handler)

    def test_joint_reasoner_one_gene(self):
        # Specify evidence
        evidence = {'ENSEMBL:ENSG00000155657': 'True'}
        # Specify targets
        dynamic_targets = {
            "EFO:0000714": {
                "op": '>=',
                "value": 1000
            }
        }
        # Setup query
        query = Query(
            evidence=evidence,
            dynamic_targets=dynamic_targets
        )
        query =  self.joint_reasoner.run_query(query)

    def test_joint_reasoner_one_gene_one_drug(self):
        # Specify evidence
        evidence = {
            'ENSEMBL:ENSG00000155657': 'True',
            'CHEMBL:CHEMBL83': 'True',
        }
        # Specify targets
        dynamic_targets = {
            "EFO:0000714": {
                "op": '>=',
                "value": 1000
            }
        }
        # Setup query
        query = Query(
            evidence=evidence,
            dynamic_targets=dynamic_targets
        )
        query =  self.joint_reasoner.run_query(query)

class TestDynamicReasoner(unittest.TestCase):

    def setUp(self):
        self.bkb_handler = BkbDataHandler(
                 bkb_major_version='coulomb',
                 bkb_minor_version='1.0',
        )
        self.dynamic_reasoner = ChpDynamicReasoner(self.bkb_handler)

    def test_dynamic_reasoner_one_gene(self):
        # Specify evidence
        evidence = {'_ENSEMBL:ENSG00000155657': 'True'}
        # Specify targets
        dynamic_targets = {
            "EFO:0000714": {
                "op": '>=',
                "value": 1000
            }
        }
        # Setup query
        query = Query(
            evidence=evidence,
            dynamic_targets=dynamic_targets
        )
        query =  self.dynamic_reasoner.run_query(query)
        query.result.summary(include_contributions=False)

    def test_dynamic_reasoner_one_gene_one_drug(self):
        # Specify evidence
        evidence = {
            '_ENSEMBL:ENSG00000155657': 'True',
            'CHEMBL:CHEMBL83': 'True',
        }
        # Specify targets
        dynamic_targets = {
            "EFO:0000714": {
                "op": '>=',
                "value": 1000
            }
        }
        # Setup query
        query = Query(
            evidence=evidence,
            dynamic_targets=dynamic_targets
        )
        query =  self.dynamic_reasoner.run_query(query)
        query.result.summary(include_contributions=False)

    def test_dynamic_reasoner_two_gene_one_drug(self):
        # Specify evidence
        evidence = {
            '_ENSEMBL:ENSG00000155657': 'True',
            '_ENSEMBL:ENSG00000241973': 'True',
            'CHEMBL:CHEMBL83': 'True',
        }
        # Specify targets
        dynamic_targets = {
            "EFO:0000714": {
                "op": '>=',
                "value": 1000
            }
        }
        # Setup query
        query = Query(
            evidence=evidence,
            dynamic_targets=dynamic_targets
        )
        query =  self.dynamic_reasoner.run_query(query)
        query.result.summary(include_contributions=False)

    def test_dynamic_reasoner_one_drug_survival(self):
        # Specify evidence
        evidence = {
            '_CHEMBL:CHEMBL83': 'True',
        }
        # Specify targets
        dynamic_targets = {
            "EFO:0000714": {
                "op": '>=',
                "value": 1000
            }
        }
        # Setup query
        query = Query(
            evidence=evidence,
            dynamic_targets=dynamic_targets
        )
        query =  self.dynamic_reasoner.run_query(query, bkb_type='drug')
        query.result.summary(include_contributions=False)

    def test_dynamic_reasoner_two_drug_survival(self):
        # Specify evidence
        evidence = {
            '_CHEMBL:CHEMBL83': 'True',
            '_CHEMBL:CHEMBL1201247': 'True',
        }
        # Specify targets
        dynamic_targets = {
            "EFO:0000714": {
                "op": '>=',
                "value": 1000
            }
        }
        # Setup query
        query = Query(
            evidence=evidence,
            dynamic_targets=dynamic_targets
        )
        query =  self.dynamic_reasoner.run_query(query, bkb_type='drug')
        query.result.summary(include_contributions=False)
