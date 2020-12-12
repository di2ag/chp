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
import matplotlib.pyplot as plt

class PathwayProcessor:
    def __init__(self):
        self.bkfs = list()
        self.pathways = list()

    class Pathway:
        def __init__(self, pathwayID, genes):
            self.pathwayID = pathwayID
            self.genes = genes

    def processPathways(self, pathwaysFile):
        #print("reading pathways file...")
        with open(pathwaysFile, 'r') as csv_file:
            reader = csv.reader(csv_file)
            rows = [row for row in reader]

        headers = rows[0]
        for row in tqdm.tqdm(rows[1:], desc='Reading pathways file'):
            genes = list()
            pathwayName = row[0].replace("/","-")
            pathwayName = pathwayName.replace(" ","_")
            for i in range(0, len(row)):
                if row[i] != "" and 'Low' in headers[i]:
                    # splits gene_low/high to just get gene
                    header = headers[i].split('_')
                    # reading low and high (row[i+1])
                    genes.append((header[0],row[i],row[i+1]))
            self.pathways.append(self.Pathway(pathwayName, genes))


    def processPathwayBKF(self, exhaustiveOr=False):
        assert len(self.pathways) > 0, "Have not processed pathways yet"
        for pathway in tqdm.tqdm(self.pathways, desc='Processing pathway BKFs'):
            bkf = BKB(name=pathway.pathwayID)

            pathwayActiveComp_idx = bkf.addComponent('{}_active='.format(pathway.pathwayID))
            pathwayActiveTrue_idx = bkf.addComponentState(pathwayActiveComp_idx, 'True')

            if exhaustiveOr:
                print("not implemented")
            else:
                geneSelectorComp_idx = bkf.addComponent("Gene_combo=")
                for gene in pathway.genes: #gene[0] = geneName, gene[1] = low value gene[2] = high value
                    geneCombo_idx = bkf.addComponentState(geneSelectorComp_idx, gene[0])

                    statConditionComp_idx = bkf.addComponent('mu-STD<={}<=mu+STD='.format(gene[0]))
                    statConditionTrue_idx = bkf.addComponentState(statConditionComp_idx, 'True')

                    bkf.addSNode(BKB_S_node(statConditionComp_idx, statConditionTrue_idx, 1.0))

                    bkf.addSNode(BKB_S_node(geneSelectorComp_idx, geneCombo_idx, 1.0, [(statConditionComp_idx, statConditionTrue_idx)]))

                    bkf.addSNode(BKB_S_node(pathwayActiveComp_idx, pathwayActiveTrue_idx, 1.0, [(geneSelectorComp_idx, geneCombo_idx)]))

                self.bkfs.append(bkf)
        assert len(self.pathways) > 0, "Have not processed pathways yet"

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
    PP = PathwayProcessor()
    PP.processPathways('data/pathways_analysis.csv')
    PP.processPathwayBKF()
    PP.BKFsToFile('PatientPathwayBKFs/')
