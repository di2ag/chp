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
from pybkb.core.python_base.reasoning import checkMutex
from pybkb.core.cpp_base.reasoning import updating

sys.path.append('/home/cyakaboski/src/python/projects/bkb-pathway-provider/core')

from reasoner import Reasoner
from query import Query

PATIENTS = ['Patient{}'.format(i) for i in range(50)]
GENES = ['Gene{}'.format(i) for i in range(10)]
GENE_VARIANTS = ['Variant{}'.format(i) for i in range(2)]

#-- Add demographic evidence.
GENDER_OPTIONS = ('Male', 'Female')
AGE_RANGE_YRS = (20, 90)
SURIVAL_YRS = (0, 10)

patient_data = {patient: {'Gender': random.choice(GENDER_OPTIONS),
                          'Age': random.randrange(AGE_RANGE_YRS[0], AGE_RANGE_YRS[1]),
                          'Survival': random.randrange(SURIVAL_YRS[0], SURIVAL_YRS[1])}
                for patient in PATIENTS}
bkfs = list()

#-- Make BKB frags
for patient in PATIENTS:
    bkf = BKB()

    for gene in GENES:
        if random.choice([True, False]):
            #-- Setup Gene i component.
            comp_idx = bkf.addComponent('mut_{}'.format(gene))
            stateTrue_idx = bkf.addComponentState(comp_idx, 'True')

            bkf.addSNode(BKB_S_node(init_component_index=comp_idx,
                                    init_state_index=stateTrue_idx,
                                    init_probability=1))
            variant_comp_idx = bkf.addComponent('mut-var_{}'.format(gene))
            random_variant = random.choice(GENE_VARIANTS)
            variant_state_idx = bkf.addComponentState(variant_comp_idx, random_variant)
            bkf.addSNode(BKB_S_node(init_component_index=variant_comp_idx,
                                    init_state_index=variant_state_idx,
                                    init_probability=1,
                                    init_tail=[(comp_idx, stateTrue_idx)]))
            if 'Patient_Genes' not in patient_data[patient]:
                patient_data[patient]['Patient_Genes'] = [gene]
            else:
                patient_data[patient]['Patient_Genes'].append(gene)
            if 'Patient_Gene_Variants' not in patient_data[patient]:
                patient_data[patient]['Patient_Gene_Variants'] = ['{}-{}'.format(gene, random_variant)]
            else:
                patient_data[patient]['Patient_Gene_Variants'].append('{}-{}'.format(gene, random_variant))

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

patient_data_hash = dict()
for patient, dict_ in patient_data.items():
    patient_data_hash[hash(patient)] = dict_

print(patient_data)

#-- Denote Demographic evidence.
demo_ev = [('Age', '>=', 50),
           ('Gender', '==', 'Male')]

demo_tar = [('Survival', '>=', 2)]

reasoner = Reasoner(fused_bkb, patient_data_hash)
random_genes = {'mut_{}'.format(gene): 'True' for gene in GENES[:5]}
query0 = Query(evidence=random_genes,
               targets=list(),
               meta_evidence=demo_ev,
               meta_targets=demo_tar,
               type='updating')

query0 = reasoner.analyze_query(query0, target_strategy='topological', interpolation='independence')
query0.getReport()
#print(checkMutex(query0.bkb))
#query0.bkb.makeGraph()
