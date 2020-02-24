import os
import sys
import random
import itertools
from operator import ge, le, eq
random.seed(112)

from pybkb import bayesianKnowledgeBase as BKB
from pybkb.core.common.bayesianKnowledgeBase import BKB_I_node, BKB_component, BKB_S_node
from pybkb.core.python_base.fusion import fuse
from pybkb.core.cpp_base.reasoning import updating

PATIENTS = ['Patient_{}'.format(i) for i in range(5)]
NUM_GENES = 2

bkfs = list()

#-- Make BKB frags
for j, _ in enumerate(PATIENTS):
    bkf = BKB()

    for i in range(NUM_GENES):
        #-- Setup Gene i component.
        comp = BKB_component('Gene{}_mutated'.format(i))
        stateTrue = BKB_I_node(init_name='True', init_component=comp)
        stateFalse = BKB_I_node(init_name='False', init_component=comp)
        bkf.addComponent(comp)
        bkf.addComponentState(comp, stateTrue)
        bkf.addComponentState(comp, stateFalse)

        prob = random.choice([0,1])
        s_node_prior_1 = BKB_S_node(init_component=comp,
                                    init_state=stateTrue,
                                    init_probability=prob)
        s_node_prior_2 = BKB_S_node(init_component=comp,
                                    init_state=stateFalse,
                                    init_probability=1 - prob)
        bkf.addSNode(s_node_prior_1)
        bkf.addSNode(s_node_prior_2)

    bkfs.append(bkf)

#-- Fuse patients together.
fused_bkb = fuse(bkfs,
                 [1 for _ in range(len(PATIENTS))],
                 PATIENTS,
                 working_dir=os.getcwd())

#fused_bkb.makeGraph(layout='neato')

#-- Add demographic evidence.
GENDER_OPTIONS = ('Male', 'Female')
AGE_RANGE_YRS = (20, 90)
SURIVAL_YRS = (0, 10)

patient_data = {patient: {'Gender': random.choice(GENDER_OPTIONS),
                          'Age': random.randrange(AGE_RANGE_YRS[0], AGE_RANGE_YRS[1]),
                          'Survival': random.randrange(SURIVAL_YRS[0], SURIVAL_YRS[1])}
                for patient in PATIENTS}
print(patient_data)

#-- Denote Demographic evidence.
demo_ev = [('Age', ge, 50),
           ('Gender', eq, 'Male'),
           ('Survival', ge, 2)]

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

#-- Construct demographic evidence
pop_stats = dict()
evidence = dict()
demo_ev.reverse()
for i, ev in enumerate(demo_ev):
    prop, op, val = ev
    matched_patients = set()
    pop_count_true = 0
    for patient in PATIENTS:
        if op(patient_data[patient][prop], val):
            matched_patients.add(patient)
            pop_count_true += 1
    pop_stats[prop] = float(pop_count_true / len(PATIENTS))

    if op == ge:
        op_str = '>='
    elif op == le:
        op_str = '<='
    else:
        op_str = '=='

    comp = BKB_component('{} {} {}'.format(prop, op_str, val))
    inode_true = BKB_I_node(init_name='True', init_component=comp)
    inode_false = BKB_I_node(init_name='False', init_component=comp)
    fused_bkb.addComponent(comp)
    fused_bkb.addComponentState(comp, inode_true)
    fused_bkb.addComponentState(comp, inode_false)

    #-- Set evidence
    evidence[comp.name] = inode_true.name

    #-- If this is last piece of evidence, make prior snodes. 
    if i == 0:
        #-- Add prior snode
        snode_1 = BKB_S_node(init_component=comp, init_state=inode_true, init_probability=pop_stats[prop])
        snode_2 = BKB_S_node(init_component=comp, init_state=inode_false, init_probability=1-pop_stats[prop])
        fused_bkb.addSNode(snode_1)
        fused_bkb.addSNode(snode_2)


    #-- Connect pieces of demographic evidence or targets together via chain rule.
    else:
        #-- Collect all previous components
        prev_comps = list()
        for j in range(i):
            ev_j = demo_ev[-j]
            prop_j, op_j, val_j = ev_j

            if op_j == ge:
                op_str_j = '>='
            elif op_j == le:
                op_str_j = '<='
            else:
                op_str_j = '=='

            prev_comps.append(fused_bkb.findComponent('{} {} {}'.format(prop_j, op_str_j, val_j)))


        #-- Calculate joint probabilities from data
        tail_product = [combo for combo in itertools.product(demo_ev[:i], [True, False])]
        #print(tail_product)
        tail_combos = list()
        for combo in itertools.combinations(tail_product, r=i):
            combo = list(combo)
            #print(combo)
            prop_set = set()
            for ev_state in combo:
                #print(ev_state)
                ev_, state = ev_state
                prop_, op_, val_ = ev_
                prop_set.add(prop_)
            if len(prop_set) == len(combo):
                tail_combos.append(combo)
        #print(tail_combos)
        combos = [combo for combo in itertools.product([(ev, True), (ev, False)], tail_combos)]

        print(combos)
        counts = list()
        for combo in combos:
            head, tail = combo
            #-- Put head and tail in one list
            combo = [head] + tail
            count = 0
            for patient in PATIENTS:
                truth = list()
                for ev_state in combo:
                    #print(ev_state)
                    ev_, state = ev_state
                    prop_, op_, val_ = ev_
                    res = op_(patient_data[patient][prop_], val_)
                    truth.append(res == state)
                if all(truth):
                    count += 1
            counts.append(count)

        probs = [float(count / len(PATIENTS)) for count in counts]
        print(probs)
        #-- Setup each S-node
        for j, combo in enumerate(combos):
            #print(combo)
            head, tail = combo

            #-- Process head
            head_ev, state = head
            prop_head, op_head, val_head = head_ev

            if op_head == ge:
                op_str_head = '>='
            elif op_head == le:
                op_str_head = '<='
            else:
                op_str_head = '=='

            comp_head = fused_bkb.findComponent('{} {} {}'.format(prop_head, op_str_head, val_head))
            if state:
                i_node_head = fused_bkb.findINode(comp_head, 'True')
            else:
                i_node_head = fused_bkb.findINode(comp_head, 'False')

            #print(comp_head.name)
            #print(i_node_head.name)

            #-- Process Tail
            processed_tail = list()
            for tail_ in tail:
                tail_ev, state = tail_
                prop_tail, op_tail, val_tail = tail_ev

                if op_tail == ge:
                    op_str_tail = '>='
                elif op_tail == le:
                    op_str_tail = '<='
                else:
                    op_str_tail = '=='

                comp_tail = fused_bkb.findComponent('{} {} {}'.format(prop_tail, op_str_tail, val_tail))
                if state:
                    i_node_tail = fused_bkb.findINode(comp_tail, 'True')
                else:
                    i_node_tail = fused_bkb.findINode(comp_tail, 'False')
                processed_tail.append((comp_tail, i_node_tail))
            #print(processed_tail)
            #-- Add Snode
            if probs[j] > 0:
                fused_bkb.addSNode(BKB_S_node(init_component=comp_head, init_state=i_node_head, init_probability=probs[j], init_tail=processed_tail))

        #-- Process Sources
        #-- If first piece of evidence connect to all patients.
        if i == len(demo_ev) - 1:
            #-- Attach to source nodes
            for patient, src_name in patient_name_src_map.items():
                for src_comp in src_components:
                    src_state = fused_bkb.findINode(src_comp, src_name)

                    #-- Find Prior snode, capture probability and remove.
                    cidx = fused_bkb.getComponentIndex(src_comp)
                    iidx = src_comp.getStateIndex(src_state)
                    s_node = fused_bkb.S_nodes_by_head[cidx][iidx][0]
                    fuse_prob = s_node.probability
                    fused_bkb.removeSNode(s_node)

                    if patient in matched_patients:
                        fused_bkb.addSNode(BKB_S_node(init_component=src_comp,
                                                      init_state=src_state,
                                                      init_probability=fuse_prob,
                                                      init_tail = [(comp, inode_true)]))
                        #fused_bkb.addSNode(BKB_S_node(init_component=src_comp,
                        #                              init_state=src_state,
                        #                              init_probability=0,
                        #                              init_tail = [(comp, inode_false)]))
                    else:
                        #fused_bkb.addSNode(BKB_S_node(init_component=src_comp,
                        #                              init_state=src_state,
                        #                              init_probability=0,
                        #                              init_tail = [(comp, inode_true)]))
                        fused_bkb.addSNode(BKB_S_node(init_component=src_comp,
                                                      init_state=src_state,
                                                      init_probability=fuse_prob,
                                                      init_tail = [(comp, inode_false)]))




        '''
            #-- Treating as non-independent
            count_i_j = 0
            count_i_notj = 0
            count_noti_j = 0
            count_noti_notj = 0
            count_j = 0
            for patient in PATIENTS:
                if op(patient_data[patient][prop], val) and op_j(patient_data[patient][prop_j], val_j):
                    count_i_j += 1
                    count_j += 1
                elif op(patient_data[patient][prop], val) and not op_j(patient_data[patient][prop_j], val_j):
                    count_i_notj += 1
                elif not op(patient_data[patient][prop], val) and op_j(patient_data[patient][prop_j], val_j):
                    count_noti_j += 1
                    count_j += 1
                else:
                    count_noti_notj += 1
            pop_stats[(prop, prop_j)] = (float(count_i_j / count_j),
                                         float(count_i_notj / (len(PATIENTS) - count_j)),
                                         float(count_noti_j / count_j),
                                         float(count_noti_notj / (len(PATIENTS) - count_j)))

            comp_j = fused_bkb.findComponent('{} {} {}'.format(prop_j, op_str_j, val_j))
            #print([comp_name for comp_name in fused_bkb.components_index])
            #print('{} {} {}'.format(prop_j, op_str_j, val_j))
            inode_true_j = fused_bkb.findINode(comp_j, 'True')
            inode_false_j = fused_bkb.findINode(comp_j, 'False')

            #-- Add S-nodes
            fused_bkb.addSNode(BKB_S_node(init_component=comp_j, init_state=inode_true_j, init_probability=pop_stats[(prop, prop_j)][0]/(2*(j+1)), init_tail=[(comp, inode_true)]))
            fused_bkb.addSNode(BKB_S_node(init_component=comp_j, init_state=inode_false_j, init_probability=pop_stats[(prop, prop_j)][1]/(2*(j+1)), init_tail=[(comp, inode_true)]))
            fused_bkb.addSNode(BKB_S_node(init_component=comp_j, init_state=inode_true_j, init_probability=pop_stats[(prop, prop_j)][2]/(2*(j+1)), init_tail=[(comp, inode_false)]))
            fused_bkb.addSNode(BKB_S_node(init_component=comp_j, init_state=inode_false_j, init_probability=pop_stats[(prop, prop_j)][3]/(2*(j+1)), init_tail=[(comp, inode_false)]))

        #-- Delete the last source_j's prior snodes
        cidx = fused_bkb.getComponentIndex(comp_j)
        for iidx in range(comp_j.getNumberStates()):
            s_nodes = fused_bkb.S_nodes_by_head[cidx][iidx]
            for s_node_ in s_nodes:
                if s_node_.getNumberTail() == 0:
                    fused_bkb.removeSNode(s_node_)
'''
#-- Set targets
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

