import os
import sys
import random

from pybkb import bayesianKnowledgeBase as BKB

#-- Change to your local data folder
NCATS_DIR = '/home/cyakaboski/src/python/projects/bkb-pathway-provider'

sys.path.append(os.path.join(NCATS_DIR, 'core'))

from query import Query
from reasoner import Reasoner

if __name__ == '__main__':
    bkb = BKB()
    bkb.load('/home/public/data/ncats/fullBKFs/patientBKFs/fusion.bkb')

    g = bkb.makeGraph(show=False, layout='neato')

    ##-- Select Random Evidence
    #num_evidence = 1
    #evidence = dict()
    #for _ in range(num_evidence):
    #    rand_rv = random.choice(bkb.components)
    #    while rand_rv.name in evidence:
    #        rand_rv = random.choice(bkb.components)
    #    #-- Get a random state
    #    rand_state = random.choice(rand_rv.states)
    #    evidence[rand_rv.name] = rand_state.name

    #print('Got Evidence')

    ##-- Select Random Targets
    #num_targets = 2
    #targets = list()
    #for _ in range(num_targets):
    #    rand_target = random.choice(bkb.components)
    #    while rand_target.name in targets or rand_target.name[:3] != 'mut':
    #        rand_target = random.choice(bkb.components)
    #    targets.append(rand_target.name)

    #print('Got Targets')

    reasoner = Reasoner(bkb, None)
    reasoner.set_src_metadata(os.path.join(NCATS_DIR,'core','src_dict.pik'))

    query0 = Query(evidence=dict(),
                   #targets=targets,
                   type='revision',
                   meta_evidence=[('Age_of_Diagnosis', '>=', 15000),
                                  ('Gender', '==', 'FEMALE')]
                  )

    query0 = reasoner.analyze_query(query0)
    query0.getReport()
