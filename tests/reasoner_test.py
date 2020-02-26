import os
import sys
import pickle

sys.path.append('/home/cyakaboski/src/python/projects/bkb-pathway-provider/core')

from reasoner import Reasoner
from query import Query
from pybkb.core.common.bayesianKnowledgeBase import bayesianKnowledgeBase as BKB

#-- Load fused bkb
fused_bkb = BKB()
fused_bkb.load('10PatientFusion.bkb')

#-- Instantiate Reasoner
reasoner = Reasoner(fused_bkb, None)
reasoner.set_src_metadata('patient_data.pk')

#print(reasoner.metadata)

#-- Make query
query1 = Query(evidence=dict(),
               targets=list(),
               meta_evidence=[('Age_of_Diagnosis','>=', 20000)],
               meta_targets=[('Survival_Time', '>=', 500)],
               type='updating')

query1 = reasoner.analyze_query(query1)
query1.getReport()
