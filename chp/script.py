'''
Source code developed by DI2AG.
Thayer School of Engineering at Dartmouth College
Authors:    Dr. Eugene Santos, Jr
            Mr. Chase Yakaboski,
            Mr. Gregory Hyde,
            Dr. Keum Joo Kim
'''


#-- This test case took 35,000 seconds = 9.7 hours

#-- Import your version of BKB and reasoning appropriately. Below is how I do it with pybkb.
from pybkb.common.bayesianKnowledgeBase import bayesianKnowledgeBase as BKB
from pybkb.python_base.reasoning import updating

#-- Other imports
import pickle

#-- Load in evidence and targets
with open('query_evid_targs.pk', 'rb') as f_:
    data = pickle.load(f_)

#-- Pickled as dictionary.
evidence = data['evidence']
targets = data['targets']

#-- Load query BKB
bkb = BKB()
bkb.load('query_bkb_137.bkb')

#-- Run Updating
res = updating(bkb, evidence, targets)
