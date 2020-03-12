import os
import sys
import random
import itertools
from operator import ge, le, eq
import copy
random.seed(116)

from pybkb import bayesianKnowledgeBase as BKB
from pybkb.core.common.bayesianKnowledgeBase import BKB_I_node, BKB_component, BKB_S_node
from pybkb.core.python_base.fusion import fuse

sys.path.append('/home/cyakaboski/src/python/projects/bkb-pathway-provider/core')

from reasoner import Reasoner, _constructSNodesByHead
from query import Query

PATIENTS = ['Patient_{}'.format(i) for i in range(10)]
NUM_GENES = 5

bkfs = list()

#-- Make BKB frags
for j, _ in enumerate(PATIENTS):
    bkf = BKB()

    for i in range(NUM_GENES):
        #-- Setup Gene i component.
        comp_idx = bkf.addComponent('Gene{}_mutated'.format(i))
        stateTrue_idx = bkf.addComponentState(comp_idx, 'True')

        if random.choice([True, False]):
            bkf.addSNode(BKB_S_node(init_component_index=comp_idx,
                                    init_state_index=stateTrue_idx,
                                    init_probability=1))
    bkfs.append(bkf)

#-- Fuse patients together.
fused_bkb = fuse(bkfs,
                 [1 for _ in range(len(PATIENTS))],
                 [str(hash(patient)) for patient in PATIENTS],
                 working_dir=os.getcwd())

#-- Impose Structure while maintiaining mutual exclusivity.
def impose_structure(bkb, structure):
    bkb_ = copy.deepcopy(bkb)
    S_nodes_by_head = _constructSNodesByHead(bkb_)
    for gene1, dependecies in structure.items():
        gene1_comp_idx = bkb_.getComponentIndex(gene1)
        gene1_state_idx = bkb_.getComponentINodeIndex(gene1_comp_idx, 'True')
        #-- Gather all source Inodes
        src_nodes = list()
        for snode in S_nodes_by_head[gene1_comp_idx][gene1_state_idx]:
            for tail_idx in range(snode.getNumberTail()):
                tail = snode.getTail(tail_idx)
                src_nodes.append(tail)
            bkb_.removeSNode(snode)

        #-- Add in edge with scr_nodes attached
        for gene2, (add, prob) in dependecies.items():
            gene2_comp_idx = bkb_.getComponentIndex(gene2)
            gene2_state_idx = bkb_.getComponentINodeIndex(gene2_comp_idx, 'True')
            if add:
                for src_comp_idx, src_state_idx in src_nodes:
                    bkb_.addSNode(BKB_S_node(init_component_index=gene1_comp_idx,
                                                  init_state_index=gene1_state_idx,
                                                  init_probability=prob,
                                                  init_tail=[(src_comp_idx, src_state_idx), (gene2_comp_idx, gene2_state_idx)]))
    return bkb_


structure = {'Gene1_mutated':{'Gene2_mutated': (True, 1)}}
fused_bkb_wStruct = impose_structure(fused_bkb, structure)
#fused_bkb_wStruct.makeGraph()

#fused_bkb.makeGraph(layout='neato')

#-- Add demographic evidence.
GENDER_OPTIONS = ('Male', 'Female')
AGE_RANGE_YRS = (20, 90)
SURIVAL_YRS = (0, 10)

patient_data = {patient: {'Gender': random.choice(GENDER_OPTIONS),
                          'Age': random.randrange(AGE_RANGE_YRS[0], AGE_RANGE_YRS[1]),
                          'Survival': random.randrange(SURIVAL_YRS[0], SURIVAL_YRS[1])}
                for patient in PATIENTS}
patient_data_hash = dict()
for patient, dict_ in patient_data.items():
    patient_data_hash[hash(patient)] = dict_

print(patient_data)

#-- Denote Demographic evidence.
demo_ev = [('Age', '>=', 50),
           ('Gender', '==', 'Male')]

demo_tar = [('Survival', '>=', 2)]

reasoner1 = Reasoner(fused_bkb, None)
reasoner1.metadata = patient_data_hash
reasoner1.cpp_reasoning = True

reasoner2 = Reasoner(fused_bkb_wStruct, None)
reasoner2.metadata = patient_data_hash
reasoner2.cpp_reasoning = True

query0 = Query(evidence={'Gene{}_mutated'.format(i):'True' for i in range(NUM_GENES)},
               targets=list(),
               meta_evidence=demo_ev,
               meta_targets=demo_tar,
               type='updating')

query01 = reasoner1.analyze_query(copy.deepcopy(query0))
query02 = reasoner2.analyze_query(copy.deepcopy(query0))
query01.getReport()
query02.getReport()
#query01.bkb.makeGraph()
query02.bkb.makeGraph()

