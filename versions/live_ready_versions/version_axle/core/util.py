import csv
import numpy as np
import random

from pybkb import bayesianKnowledgeBase as BKB
from pybkb import BKB_S_node, BKB_component, BKB_I_node


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
