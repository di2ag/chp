'''
Source code developed by DI2AG.
Thayer School of Engineering at Dartmouth College
Authors:    Dr. Eugene Santos, Jr
            Mr. Chase Yakaboski,
            Mr. Gregory Hyde,
            Dr. Keum Joo Kim
'''


import csv
import numpy as np
import random
import pandas as pd
import tqdm
from operator import ge, le, eq

from pybkb.common.bayesianKnowledgeBase import bayesianKnowledgeBase as BKB
from pybkb.common.bayesianKnowledgeBase import BKB_S_node, BKB_component, BKB_I_node


def convert_patient_dict_to_dataframe(patient_dict,
                                      include_genes=False,
                                      include_gene_variants=False,
                                      include_drugs=False,
                                      other_info=['Age_of_Diagnosis', 'Gender', 'Cancer_Type', 'PathN', 'PathT', 'PathM', 'Survival_Time'],
                                      encode_strings=False):
    genes = set()
    gene_variants = set()
    drugs = set()
    for patient, data_dict in patient_dict.items():
        genes.update(set(data_dict['Patient_Genes']))
        gene_variants.update(set(data_dict['Patient_Gene_Variants']))
        drugs.update(set(data_dict['Drug_Name(s)']))

    def helper_fn(collection, name):
        _data = dict()
        for item in tqdm.tqdm(collection, desc='Building out columns for {}'.format(name), leave=False):
            info = list()
            for _, data_dict in patient_dict.items():
                if item in data_dict[name]:
                    info.append(1)
                else:
                    info.append(0)
            _data[item] = info
        return _data

    data = dict()
    if include_genes:
        data.update(helper_fn(genes, 'Patient_Genes'))
    if include_gene_variants:
        data.update(helper_fn(gene_variants, 'Patient_Gene_Variants'))
    if include_drugs:
        data.update(helper_fn(drugs, 'Drug_Name(s)'))

    encodings = dict()
    for col in tqdm.tqdm(other_info, desc='Building out other info', leave=False):
        info = list()
        for _, data_dict in patient_dict.items():
            val = data_dict[col]
            if encode_strings and type(val) == str:
                if val not in encodings:
                    encodings[val] = len(encodings)
                info.append(encodings[val])
            else:
                info.append(data_dict[col])
        data[col] = info

    df = pd.DataFrame(data=data)
    return df

def process_operator(op):
    if op == '>=':
        return ge
    elif op == '<=':
        return le
    elif op == '==':
        return eq
    else:
        raise ValueError('Unknown Operator')

class NcatsProcessor:
    def __init__(self):
        self.bkfs = list()

    def read_expression_snodes(self, expression_levels_dict, tcga_expression_file, react_gene_file, pathways_file):
        #-- Process all gene_reaction s_nodes and assume there is a header with gene names.
        with open(react_gene_file, 'r') as csv_file:
            reader = csv.reader(csv_file)
            rows = [row for row in reader]
        #-- Process Gene Names from header:
        print(rows[0])
        genes = rows[0][1:]
        genes = [genes[i][:-4] for i in range(0,len(genes),2)]
        print(genes)
        for row in rows[1:]:
            bkf = BKB(name=row[0] + ' BKF')
            reaction_head_name = row[0]
            react_component = BKB_component(reaction_head_name)
            i_node_true = BKB_I_node('True', react_component)
            i_node_false = BKB_I_node('False', react_component)
            react_component.addINode(i_node_true)
            react_component.addINode(i_node_false)
            bkf.addComponent(react_component)
            #-- Extract and Make I-nodes
            for gene_j, row_i in enumerate(range(1,len(row), 2)):
                gene_name = genes[gene_j]
                gene_low = row[row_i]
                gene_high = row[row_i + 1]
                if gene_low != '':
                    print(gene_low, gene_high)
                    num_levels = expression_levels_dict[gene_name]
                    levels = [rag for rag in np.linspace(int(gene_low), int(gene_high), num_levels)]

                    component_states = ['{} | {}'.format(levels[i], levels[i+1]) for i in range(len(levels)-1)]
                    component = BKB_component(gene_name + ' ExpressionLv')
                    for state in component_states:
                        i_node = BKB_I_node(state, component)
                        component.addINode(i_node)
                    bkf.addComponent(component)
                    for idx in range(bkf.findComponent(component.name).getNumberStates()):
                        i_node = component.getState(idx)
                        bkf.addSNode(BKB_S_node(component, i_node, random.random()))
                        bkf.addSNode(BKB_S_node(react_component, i_node_true, random.random(), [(component, i_node)]))
                        bkf.addSNode(BKB_S_node(react_component, i_node_false, random.random(), [(component, i_node)]))
            self.bkfs.append(bkf)
        print('Done')

        #-- Process Pathway BKFs
        with open(pathways_file, 'r') as csv_file:
            reader = csv.reader(csv_file)
            bkf = BKB()
            for row in reader:
                tail = row[0]
                head = row[1]
                states = ['True', 'False']
                if bkf.findComponent(tail) == -1:
                    component_tail = BKB_component(tail)
                    for state in states:
                        component_tail.addINode(BKB_I_node(state, component_tail))
                    bkf.addComponent(component_tail)
                else:
                    component_tail = bkf.findComponent(tail)

                if bkf.findComponent(head) == -1:
                    component_head = BKB_component(head)
                    for state in states:
                        component_head.addINode(BKB_I_node(state, component_head))
                    bkf.addComponent(component_head)
                else:
                    component_head = bkf.findComponent(head)

                for state1 in states:
                    for state2 in states:
                        bkf.addSNode(BKB_S_node(component_head,
                                                component_head.findState(state1),
                                                random.random(),
                                                [(component_tail, component_tail.findState(state2))])
                                    )
        self.bkfs.append(bkf)
        print('Complete')


if __name__ == '__main__':
    processor = NcatsProcessor()
    num_genes = 5
    expression_levels_dict = {'Gene_{}'.format(i): 3 for i in range(num_genes)}
    tcga_file = '../simulations/test.csv'
    react_gene_file='../simulations/test-react_gene.snodes'
    react_pathways_file = '../simulations/test-pathways.snodes'

    processor.read_expression_snodes(expression_levels_dict, tcga_file, react_gene_file, react_pathways_file)

    pathways_bkb = processor.bkfs[-1]
    pathways_bkb.name = 'Pathways BKB from Reactome'
    pathways_bkb.makeGraph(layout='neato')

    processor.bkfs[0].makeGraph()
