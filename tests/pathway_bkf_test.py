import os
import sys
import random

from pybkb import bayesianKnowledgeBase as BKB
from pybkb.core.cpp_base.fusion import fuse

#-- Change to your local data folder
NCATS_DIR = '/home/cyakaboski/src/python/projects/bkb-pathway-provider'

sys.path.append(os.path.join(NCATS_DIR, 'core'))

from query import Query
from reasoner import Reasoner
from reactomePathwayProcessor import ReactomePathwayProcessor

reactome_g2r_file = '/home/public/data/ncats/reactome/r2g_names_k1.csv'
reactome_p2p_file = '/home/public/data/ncats/reactome/p2p_k1.csv'

if __name__ == '__main__':
    processor = ReactomePathwayProcessor()
    processor.processGenePathways(reactome_g2r_file)
    processor.processHierarchyPathways(reactome_p2p_file)
    processor.processPathwayBKF()
    bkf_files, source_names = processor.BKFsToFile('/tmp/')

    bkfs = processor.bkfs

    fuse(bkfs,
         reliabilities=[random.random() for _ in range(len(bkfs))],
         source_names=source_names,
         file_prefix='react_only-')

    fused_bkb = BKB()
    fused_bkb.load('/tmp/react_only-fusion.bkb')

    #g = fused_bkb.makeGraph(show=True, layout='neato')

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

    #reasoner = Reasoner(bkb, None)
    #reasoner.set_src_metadata(os.path.join(NCATS_DIR,'core','src_dict.pik'))

    #query0 = Query(evidence=dict(),
    #               #targets=targets,
    #               type='revision',
    #               meta_evidence=[('Age_of_Diagnosis', '>=', 15000),
    #                              ('Gender', '==', 'FEMALE')]
    #              )

    #query0 = reasoner.analyze_query(query0)
    #query0.getReport()
