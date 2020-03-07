import os
import sys
import random
import itertools
from operator import ge, le, eq
random.seed(116)

from pybkb import bayesianKnowledgeBase as BKB
from pybkb.core.common.bayesianKnowledgeBase import BKB_I_node, BKB_component, BKB_S_node
from pybkb.core.python_base.fusion import fuse
from pybkb.core.python_base.fusion_collapse import collapse_sources
from pybkb.core.cpp_base.reasoning import updating

sys.path.append('/home/cyakaboski/src/python/projects/bkb-pathway-provider/core')

from reasoner import Reasoner
from query import Query

PATIENTS = ['Patient{}'.format(i) for i in range(100)]
GENES = ['Gene{}_mut'.format(i) for i in range(10)]
GENE_VARIANTS = ['Variant{}'.format(i) for i in range(2)]

bkfs = list()

#-- Make BKB frags
for j, _ in enumerate(PATIENTS):
    bkf = BKB()

    for gene in GENES:
        #-- Setup Gene i component.
        comp_idx = bkf.addComponent(gene)
        stateTrue_idx = bkf.addComponentState(comp_idx, 'True')

        if random.choice([True, False]):
            bkf.addSNode(BKB_S_node(init_component_index=comp_idx,
                                    init_state_index=stateTrue_idx,
                                    init_probability=1))
            variant_comp_idx = bkf.addComponent(gene + '_Var')
            variant_state_idx = bkf.addComponentState(variant_comp_idx, random.choice(GENE_VARIANTS))
            bkf.addSNode(BKB_S_node(init_component_index=variant_comp_idx,
                                    init_state_index=variant_state_idx,
                                    init_probability=1,
                                    init_tail=[(comp_idx, stateTrue_idx)]))
#    if j == 0:
#        bkf.makeGraph()
    bkfs.append(bkf)

#-- Fuse patients together.
fused_bkb = fuse(bkfs,
                 [1 for _ in range(len(PATIENTS))],
                 [str(hash(patient)) for patient in PATIENTS],
                 working_dir=os.getcwd())

#fused_bkb.makeGraph(layout='neato')
#col_bkb = collapse_sources(fused_bkb)
#col_bkb.makeGraph(layout='neato')

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

reasoner = Reasoner(fused_bkb, None)
reasoner.metadata = patient_data_hash

query0 = Query(evidence={'Gene1_mut':'True',
                         'Gene2_mut':'True'},
               targets=list(),
               meta_evidence=demo_ev,
               meta_targets=demo_tar,
               type='updating')

query0 = reasoner.analyze_query(query0)
query0.getReport()
#query0.bkb.makeGraph()
