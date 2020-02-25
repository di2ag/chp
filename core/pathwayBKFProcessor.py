import csv
import tqdm

from pybkb import bayesianKnowledgeBase as BKB
from pybkb import BKB_S_node, BKB_component, BKB_I_node
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

            pathwayActiveComp = BKB_component(pathway.pathwayID + "_active=")
            pathwayActiveTrue = BKB_I_node('True',pathwayActiveComp)
            bkf.addComponent(pathwayActiveComp)
            bkf.addComponentState(pathwayActiveComp, pathwayActiveTrue)

            if exhaustiveOr:
                print("not implemented")
            else:
                geneSelectorComp = BKB_component("Gene_combo=")
                bkf.addComponent(geneSelectorComp)
                for gene in pathway.genes: #gene[0] = geneName, gene[1] = low value gene[2] = high value
                    geneCombo = BKB_I_node(gene[0],geneSelectorComp)
                    bkf.addComponentState(geneSelectorComp, geneCombo)

                    statConditionComp = BKB_component("mu-STD>=" + gene[0] + "<=mu+STD=")
                    statConditionTrue = BKB_I_node('True', statConditionComp)
                    bkf.addComponent(statConditionComp)
                    bkf.addComponentState(statConditionComp, statConditionTrue)

                    bkf.addSNode(BKB_S_node(statConditionComp, statConditionTrue, 1.0))

                    bkf.addSNode(BKB_S_node(geneSelectorComp, geneCombo, 1.0, [(statConditionComp, statConditionTrue)]))

                    bkf.addSNode(BKB_S_node(pathwayActiveComp, pathwayActiveTrue, 1.0, [(geneSelectorComp, geneCombo)]))

                self.bkfs.append(bkf)
        assert len(self.pathways) > 0, "Have not processed pathways yet"

    def BKFsToFile(self, outDirect):
        bkf_files = list()
        source_names = list()
        for bkf in self.bkfs:
            print(bkf.name)
            file_name = outDirect + bkf.name + '.bkf'
            bkf.save(file_name)
            bkf_files.append(file_name)
            source_names.append(str(bkf.name))

        return bkf_files, source_names


if __name__ == '__main__':
    PP = PathwayProcessor()
    PP.processPathways('data/pathways_analysis.csv')
    PP.processPathwayBKF()
    PP.BKFsToFile('PatientPathwayBKFs/')
