'''
Source code developed by DI2AG.
Thayer School of Engineering at Dartmouth College
Authors:    Dr. Eugene Santos, Jr
            Mr. Chase Yakaboski,
            Mr. Gregory Hyde,
            Dr. Keum Joo Kim
'''


import csv
import tqdm

from pybkb.common.bayesianKnowledgeBase import bayesianKnowledgeBase as BKB
from pybkb.common.bayesianKnowledgeBase import BKB_S_node, BKB_component, BKB_I_node

class ReactomePathwayProcessor:
    def __init__(self):
        self.bkfs = list()
        self.pathways = list()

    class Pathway:
        def __init__(self, tails, head, ID):
            self.tails = tails
            self.head = head
            self.ID = ID

    def processHierarchyPathways(self, pathwaysDirectory):
        #print("Processing hierarchy pathways...")
        with open(pathwaysDirectory, 'r') as csv_file:
            reader = csv.reader(csv_file)
            rows = [row for row in reader]

        for row in tqdm.tqdm(rows, desc='Processing hierarchy pathways'):
            tails = [row[0]]
            for tail in tails:
                tail = tail.replace(" ", "_")
                tail = tail.replace("/", "-")
            head = row[1].replace(" ", "_")
            head = head.replace("/", "-")
            pathway = self.Pathway(tails,head,head+"_hier")
            self.pathways.append(pathway)
        #print("Hierarchy pathways processed")

    def processGenePathways(self, pathwaysDirectory):
        #print("Processing gene pathways...")
        with open(pathwaysDirectory, 'r') as csv_file:
            reader = csv.reader(csv_file)
            rows = [row for row in reader]

        for row in tqdm.tqdm(rows, desc='Processing gene pathways'):
            head = row[0].replace(" ", "_")
            head = head.replace("/", "-")
            tails = list()
            for cell in row[1:]:
                if cell == '':
                    continue
                elif "UniProt" not in cell:
                    continue
                if cell.split(" ")[1] not in tails:
                    tails.append(cell.split(" ")[1])
            for tail in tails:
                tail = tail.replace(" ", "_")
                tail = tail.replace("/", "-")
            pathway = self.Pathway(tails,head, head+"_reaction")
            self.pathways.append(pathway)
        #print("Gene pathways Processed")

    def processPathwayBKF(self, exhaustiveOr=False):
        assert len(self.pathways) > 0, "Have not processed pathways yet"

        for pathway in tqdm.tqdm(self.pathways, desc='Processing pathway BKFs'):
            bkf = BKB(name=pathway.ID)
            if "_hier" in pathway.ID:
                #head
                pathwayParentComp_idx = bkf.addComponent('{}_active='.format(pathway.head))
                pathwayParentTrue_idx = bkf.addComponentState(pathwayParentComp_idx, 'True')

                #tail - should only be 1 for _hier bkfs
                pathwayChildComp_idx = bkf.addComponent('{}_active='.format(pathway.tails[0]))
                pathwayChildTrue_idx = bkf.addComponentState(pathwayChildComp_idx, 'True')

                #S-node
                bkf.addSNode(BKB_S_node(pathwayChildComp_idx, pathwayChildTrue_idx, 1.0))
                bkf.addSNode(BKB_S_node(pathwayParentComp_idx, pathwayParentTrue_idx, 1.0, [(pathwayChildComp_idx, pathwayChildTrue_idx)]))
            else:
                '''
                #head
                pathwayHierComp_idx = bkf.addComponent('{}_active='.format(pathway.head))
                pathwayHierTrue_idx = bkf.addComponentState(pathwayHierComp_idx, 'True')

                #tails
                for tail in pathway.tails:
                    statConditionComp_idx = bkf.addComponent('mu-STD<={}<=mu+STD'.format(tail))
                    statConditionTrue_idx = bkf.addComponentState(statConditionComp_idx, 'True')

                    bkf.addSNode(BKB_S_node(statConditionComp_idx, statConditionTrue_idx, 1.0))
                    bkf.addSNode(BKB_S_node(pathwayHierComp_idx, pathwayHierTrue_idx, 1.0, [(statConditionComp_idx, statConditionTrue_idx)]))

            self.bkfs.append(bkf)
                '''
#=======

                if exhaustiveOr:
                    print("not implemented")
                else:
                    # pathway
                    pathwayReactionComp_idx = bkf.addComponent(pathway.head + "_active=")
                    pathwayReactionTrue_idx = bkf.addComponentState(pathwayReactionComp_idx, 'True')

                    # gene slector
                    geneSelectorComp_idx = bkf.addComponent("Gene_combo=")

                    #tails
                    for tail in pathway.tails:
                        geneCombo_idx = bkf.addComponentState(geneSelectorComp_idx, tail)

                        statConditionComp = BKB_component("mu-STD>=" + tail + "<=mu+STD=")
                        statConditionTrue = BKB_I_node('True', statConditionComp)
                        statConditionComp_idx = bkf.addComponent("mu-STD>=" + tail + "<=mu+STD=")
                        statConditionTrue_idx = bkf.addComponentState(statConditionComp_idx, 'True')

                        bkf.addSNode(BKB_S_node(statConditionComp_idx, statConditionTrue_idx, 1.0))

                        bkf.addSNode(BKB_S_node(geneSelectorComp_idx, geneCombo_idx, 1.0, [(statConditionComp_idx, statConditionTrue_idx)]))

                        bkf.addSNode(BKB_S_node(pathwayReactionComp_idx, pathwayReactionTrue_idx, 1.0, [(geneSelectorComp_idx, geneCombo_idx)]))

            self.bkfs.append(bkf)
#>>>>>>> 11e4db965d4e84060bc6f402b8980cf50af7e7d7

    def BKFsToFile(self, outDirect):
        bkf_files = list()
        source_names = list()
        for bkf in self.bkfs:
            file_name = outDirect + bkf.getName() + '.bkf'
            bkf.save(file_name)
            bkf_files.append(file_name)
            source_names.append(str(bkf.getName()))

        return bkf_files, source_names

if __name__ == '__main__':
    PP = ReactomePathwayProcessor()
    PP.processHierarchyPathways('data/reactomePaths/p2p_k1.csv')
    PP.processGenePathways('data/reactomePaths/r2g_names_k1.csv')
    PP.processPathwayBKF()
    PP.BKFsToFile('PatientPathwayBKFs/')
