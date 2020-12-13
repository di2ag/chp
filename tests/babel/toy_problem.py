import copy
import os
import sys
import random
import itertools
from operator import ge, le, eq
random.seed(116)

sys.path.append('/home/cyakaboski/src/python/projects/bkb-pathway-provider/core')

from reasoner import Reasoner
from query import Query

from pybkb import bayesianKnowledgeBase as BKB
from pybkb.core.common.bayesianKnowledgeBase import BKB_I_node, BKB_component, BKB_S_node
from pybkb.core.python_base.fusion import fuse
from pybkb.core.python_base.reasoning import updating

PATIENTS = ['Patient{}'.format(i) for i in range(5)]
NUM_GENES = 3
GENES = ['Gene{}'.format(i) for i in range(NUM_GENES)]

#-- Create patient data
GENDER_OPTIONS = ('Male', 'Female')
AGE_RANGE_YRS = (20, 90)
SURIVAL_YRS = (0, 10)

patient_data = {patient: {'Gender': random.choice(GENDER_OPTIONS),
                          'Age': random.randrange(AGE_RANGE_YRS[0], AGE_RANGE_YRS[1]),
                          'Survival': random.randrange(SURIVAL_YRS[0], SURIVAL_YRS[1]),
                          'Genes': {gene: random.choice([True, False]) for gene in GENES}}
                for patient in PATIENTS}
patient_data_hash = dict()
for patient, dict_ in patient_data.items():
    patient_data_hash[hash(patient)] = dict_

print(patient_data)
patient_bkfs = list()

#--Priors Calculate Correlations
priors = dict()
for gene in GENES:
    count = 0
    for patient in PATIENTS:
        if patient_data[patient]['Genes'][gene]:
            count += 1
    priors[gene] = float(count / len(PATIENTS))


#--Pairs Calculate Correlations
counts = dict()
for gene1 in GENES:
    for gene2 in GENES:
        count = 0
        for patient in PATIENTS:
            if patient_data[patient]['Genes'][gene1] and patient_data[patient]['Genes'][gene2]:
                count += 1
    counts[(gene1, gene2)] = float(count / len(PATIENTS)) / priors[gene2]

#-- Make BKB Patient frags
for patient in PATIENTS:
    bkf = BKB()

    for gene in patient_data[patient]['Genes']:
        #-- Setup Gene i component.
        geneComp_idx = bkf.addComponent(gene)
        geneStateTrue_idx = bkf.addComponentState(geneComp_idx, 'True')

        s_node_prior = BKB_S_node(init_component_index=geneComp_idx,
                                    init_state_index=geneStateTrue_idx,
                                    init_probability=1)
        bkf.addSNode(s_node_prior)

    patient_bkfs.append(bkf)

#-- Define Sample Pathway
PATHWAYS = ['Pathway{}'.format(i) for i in range(1)]
pathway_data = {pathway: {gene: random.random()/10 for gene in set([random.choice(GENES) for _ in range(len(GENES))])}
                for pathway in PATHWAYS}

#-- Make BKB Pathway Frags
pathway_bkfs = list()
for pathway in PATHWAYS:
    bkf = BKB()
    for gene in pathway_data[pathway]:
        geneComp_idx = bkf.addComponent(gene)
        geneState_idx = bkf.addComponentState(geneComp_idx, 'True')

        bkf.addSNode(BKB_S_node(init_component_index=geneComp_idx,
                     init_state_index=geneState_idx,
                     init_probability=pathway_data[pathway][gene]))
    pathway_bkfs.append(bkf)

bkfs = patient_bkfs + pathway_bkfs
source_names = [str(hash(patient)) for patient in PATIENTS] + PATHWAYS

#-- Fuse patients together.
fused_bkb = fuse(bkfs,
                 [1 for _ in range(len(source_names))],
                 source_names,
                 working_dir=os.getcwd())

#fused_bkb.makeGraph(layout='neato')
print(counts)

for gene in GENES:
    print('='*20 + str(' No Connection:'))
    #-- Instantiate Reasoner
    fused_bkb_ = copy.deepcopy(fused_bkb)
    reasoner = Reasoner(fused_bkb_, None)
    reasoner.metadata = patient_data_hash

    #-- Make query
    query1 = Query(evidence={gene_:'True' for gene_ in [GENES[0], GENES[2]]},
                   targets=['_Source_[{}]_'.format(gene) for gene in GENES],
                   type='updating')
    query1 = reasoner.analyze_query(query1)

    query1.getReport()
    query1.bkb._name = 'No Connection'
    #query1.bkb.makeGraph(layout='neato')


    #-- Structure Learn
    fused_bkb_1 = copy.deepcopy(fused_bkb_)
    gene0Comp_idx = fused_bkb_1.getComponentIndex('Gene0')
    gene2Comp_idx = fused_bkb_1.getComponentIndex('Gene2')
    gene0State_idx = fused_bkb_1.getComponentINodeIndex(gene0Comp_idx, 'True')
    gene2State_idx = fused_bkb_1.getComponentINodeIndex(gene2Comp_idx, 'True')

    fused_bkb_1.addSNode(BKB_S_node(init_component_index=gene0Comp_idx,
                                     init_state_index=gene0State_idx,
                                     init_probability=counts[('Gene0','Gene2')],
                                     init_tail=[(gene2Comp_idx, gene2State_idx)]))

    print('='*20 + str(' With Connection:'))
    #-- Instantiate Reasoner
    reasoner = Reasoner(fused_bkb_1, None)
    reasoner.metadata = patient_data_hash

    #-- Make query
    query1 = Query(evidence={gene_:'True' for gene_ in [GENES[0], GENES[2]]},
                   targets=['_Source_[{}]_'.format(gene) for gene in GENES],
                   type='updating')
    query1 = reasoner.analyze_query(query1)

    query1.getReport()
    query1.bkb._name = 'With Connection'
    #query1.bkb.makeGraph(layout='neato')
    break
'''
#-- Denote Demographic evidence.
demo_ev = [('Age', '>=', 50),
           ('Gender', '==', 'Male'),
           ('Survival', '>=', 2)]

#-- Create Source component map
src_components = fused_bkb.getSrcComponents()
patient_name_src_map = dict()
for comp in src_components:
    for i in range(comp.getNumberStates()):
        inode = comp.getState(i)
        for patient in PATIENTS:
            if patient in inode.name:
                patient_name_src_map[patient] = inode.name

print(patient_name_src_map)


def process_operator(op):
    if op == '>=':
        return ge
    elif op == '<=':
        return le
    elif op == '==':
        return eq
    else:
        raise ValueError('Unknown Operator')

def processDependencyHead(head, bkb):
    head_ev, state = head
    prop_head, op_str_head, val_head = head_ev
    op_head = process_operator(op_str_head)

    comp_head = bkb.findComponent('{} {} {}'.format(prop_head, op_str_head, val_head))
    if state:
        i_node_head = bkb.findINode(comp_head, 'True')
    else:
        i_node_head = bkb.findINode(comp_head, 'False')

    return comp_head, i_node_head

def processDependecyTail(tail, bkb):
    processed_tail = list()
    for tail_ in tail:
        tail_ev, state = tail_
        prop_tail, op_str_tail, val_tail = tail_ev
        op_tail = process_operator(op_str_tail)

        comp_tail = bkb.findComponent('{} {} {}'.format(prop_tail, op_str_tail, val_tail))
        if state:
            i_node_tail = bkb.findINode(comp_tail, 'True')
        else:
            i_node_tail = bkb.findINode(comp_tail, 'False')
        processed_tail.append((comp_tail, i_node_tail))
    return processed_tail

def processOptionDependency(option, option_dependencies, bkb, src_population, src_population_data):
    #-- Get consistent option combinations
    tail_product = [combo for combo in itertools.product(option_dependencies, [True, False])]
    tail_combos = list()
    for combo in itertools.combinations(tail_product, r=len(option_dependencies)):
        combo = list(combo)
        prop_set = set()
        for ev_state in combo:
            ev_, state = ev_state
            prop_, op_, val_ = ev_
            prop_set.add(prop_)
        if len(prop_set) == len(combo):
            tail_combos.append(combo)
    combos = [combo for combo in itertools.product([(option, True), (option, False)], tail_combos)]

    #-- Calculate joint probabilities from data
    counts = list()
    for combo in combos:
        head, tail = combo
        #-- Put head and tail in one list
        combo = [head] + tail
        count = 0
        for entity in src_population:
            truth = list()
            for ev_state in combo:
                ev_, state = ev_state
                prop_, op_str_, val_ = ev_
                op_ = process_operator(op_str_)
                res = op_(src_population_data[entity][prop_], val_)
                truth.append(res == state)
            if all(truth):
                count += 1
        counts.append(count)
    probs = [float(count / len(src_population)) for count in counts]

    #-- Setup each S-node
    for j, combo in enumerate(combos):
        head, tail = combo
        #-- Process head
        comp_head, i_node_head = processDependencyHead(head, bkb)
        #-- Process Tail
        processed_tail = processDependecyTail(tail, bkb)
        #-- Add Snode
        if probs[j] > 0:
            bkb.addSNode(BKB_S_node(init_component=comp_head, init_state=i_node_head, init_probability=probs[j], init_tail=processed_tail))

    return bkb

def addDemographicOption(option, bkb, src_population, src_population_data, option_dependencies=list()):
    prop, op_str, val = option
    op = process_operator(op_str)
    matched_srcs = set()
    pop_count_true = 0
    for entity in src_population:
        if op(src_population_data[entity][prop], val):
            matched_srcs.add(entity)
            pop_count_true += 1
    pop_stats[prop] = float(pop_count_true / len(src_population))

    comp = BKB_component('{} {} {}'.format(prop, op_str, val))
    inode_true = BKB_I_node(init_name='True', init_component=comp)
    inode_false = BKB_I_node(init_name='False', init_component=comp)
    bkb.addComponent(comp)
    bkb.addComponentState(comp, inode_true)
    bkb.addComponentState(comp, inode_false)

    #-- Create option dictionary
    options_dict = {comp.name: inode_true.name}

    #-- If no chain rule dependencies, just add prior s-nodes
    if len(option_dependencies) == 0:
        snode_1 = BKB_S_node(init_component=comp, init_state=inode_true, init_probability=pop_stats[prop])
        snode_2 = BKB_S_node(init_component=comp, init_state=inode_false, init_probability=1-pop_stats[prop])
        bkb.addSNode(snode_1)
        bkb.addSNode(snode_2)

    #-- Process Dependencies
    else:
        bkb = processOptionDependency(option, option_dependencies, bkb, src_population, src_population_data)

    return bkb, options_dict, matched_srcs

def addSrcConnections(comp, inode_true, inode_false, bkb, matched_srcs, src_map):
        #-- Attach to source nodes
        src_components = bkb.getSrcComponents()
        for entity_name, src_name in src_map.items():
            for src_comp in src_components:
                src_state = bkb.findINode(src_comp, src_name)

                #-- Find Prior snode, capture probability and remove.
                cidx = fused_bkb.getComponentIndex(src_comp)
                iidx = src_comp.getStateIndex(src_state)
                s_node = fused_bkb.S_nodes_by_head[cidx][iidx][0]
                prob = s_node.probability
                bkb.removeSNode(s_node)

                if entity_name in matched_srcs:
                    bkb.addSNode(BKB_S_node(init_component=src_comp,
                                                  init_state=src_state,
                                                  init_probability=prob,
                                                  init_tail = [(comp, inode_true)]))
                else:
                    bkb.addSNode(BKB_S_node(init_component=src_comp,
                                                  init_state=src_state,
                                                  init_probability=prob,
                                                  init_tail = [(comp, inode_false)]))
        return bkb

#-- Construct demographic evidence
pop_stats = dict()
evidence = dict()
demo_ev.reverse()
for i, ev in enumerate(demo_ev):
    fused_bkb, evidence_, matched_srcs = addDemographicOption(ev, fused_bkb, PATIENTS, patient_data, option_dependencies=demo_ev[:i])
    evidence.update(evidence_)

#-- Process Sources
#-- If first piece of evidence connect to all patients.
comp = fused_bkb.findComponent('{} {} {}'.format(ev[0], ev[1], ev[2]))
inode_true = fused_bkb.findINode(comp, 'True')
inode_false = fused_bkb.findINode(comp, 'False')
fused_bkb = addSrcConnections(comp, inode_true, inode_false, fused_bkb, matched_srcs, patient_name_src_map)


#-- Set targets
print('Evidence:')
print(evidence)
targets1 = ['Gene1_mutated']
targets2 = ['Age >= 50']

#-- Reason
res = updating(fused_bkb,
               evidence=evidence,
               targets=targets1,
               working_dir=os.getcwd())
res.summary()

res = updating(fused_bkb,
               evidence=dict(),
               targets=targets2,
               working_dir=os.getcwd())
res.summary()

fused_bkb.name = 'Three piece of demographic evidence'
fused_bkb.makeGraph()

'''
