import csv
from pybkb import bayesianKnowledgeBase as BKB
from pybkb import BKB_S_node, BKB_component, BKB_I_node
import matplotlib.pyplot as plt

class PatientProcessor:
    def __init__(self):
        # each patient BKF
        self.bkfs = list()
        self.patients = list()
        self.clinicalCollected = False

    class Patient:
        # will need to account for mutatedReads
        def __init__(self, patID, cancer, mutatedGenes):
            # patient gene data
            self.patientID = patID
            self.cancerType = cancer
            self.mutatedGenes = mutatedGenes
            self.mutatedGeneReads = list()
            self.geneReadsStats = list() #implement

            # clinical later
            self.ageDiagnos = -1
            self.gender = ""
            self.survivalTime = -1

            # add drug info later
            # add radiation info later

    # processes patient gene data
    def processPatientGeneData(self, patientMutGenes, patientMutReads, geneReadStats):
        print("Reading file 1 - patient genes")
        with open(patientMutGenes, 'r') as csv_file:
            reader = csv.reader(csv_file)
            patientGeneDict = dict()
            # row[0] = cancer type row[1] = patientID, row[6] = mutated gene
            rows = [(row[0],row[1],row[6]) for row in reader]

        print("Processing file 1...")
        # create all patients
        mutatedPatGenes = list()
        # skip header
        cancerType = rows[1][0]
        patID = rows[1][1]
        for row in rows[1:]:
            if patID != row[1]:
                p = self.Patient(patID,cancerType,mutatedPatGenes)
                self.patients.append(p)
                patID = row[1]
                cancerType = row[0]
                mutatedPatGenes = list()
            mutatedPatGenes.append(row[2])
        csv_file.close()

        print("Reading file 2 - gene read stats")
        #implement
        print("Processing file 2...")
        #implement

        print("Reading file 3 - patient gene reads")
        # get gene reads from tumor tissue

        with open(patientMutReads, 'r') as csv_file:
            reader = csv.reader(csv_file)
            geneReadsDict = dict()
            # row[1] = patientID, row[6] = gene, row[7] = reads
            for row in reader:
                geneReadsDict[str(row[1])+str(row[6])] = row[7]
        print("Processing file 3...")
        for p in self.patients:
            geneMismatch = []
            for gene in p.mutatedGenes[:]:
                if str(p.patientID)+str(gene) in geneReadsDict:
                    p.mutatedGeneReads.append(geneReadsDict[str(p.patientID)+str(gene)])
                #gene read information not found
                else:
                    geneMismatch.append(gene)
                    p.mutatedGenes.remove(gene)
            if len(geneMismatch) > 0:
                print("Warning - Patient", p.patientID, "has", len(geneMismatch), "missing genes")

    def processClinicalData(self, clinicalData):
        assert len(self.patients) > 0, "Patient and gene data have not been read in yet. Call processPatientGeneData(patientMutGenesFile, patientMutReadsFile) first, before secondary processing functions are called."
        print("Reading clinical data...")
        with open(clinicalData, 'r') as csv_file:
            reader = csv.reader(csv_file)
            patientClinicalDict = dict()
            # row[0] = cancerType, row[1] = patientID, row[2] = age of diagnosis, row[3] = gender, row[7] = survival time
            for row in reader:
                patientClinicalDict[str(row[0])+str(row[1])] = (str(row[2]),str(row[3]),str(row[7]))

        for p in self.patients:
            if str(p.cancerType)+str(p.patientID) in patientClinicalDict:
                clinical = patientClinicalDict[str(p.cancerType)+str(p.patientID)]
                p.ageDiagnos = clinical[0]
                p.gender = clinical[1]
                p.survivalTime = clinical[2]
            else:
                print("WARNING - Patient", p, "clinical data does not exist")
        self.clinicalCollected = True
        print("Clinical data processed")

    # should be called after all Patient objects have been cosntructed
    def processPatientBKF(self):
        print("Writing BKF files...")
        assert len(self.patients) > 0, "Have not processed Patient and gene data yet."

        for pat in self.patients:
            bkf = BKB(name=pat.patientID + '_BKF')

            for gene in pat.mutatedGenes:
                # gene
                mutGeneComp = BKB_component("mut_" + gene + "=")
                iNodeGeneMut = BKB_I_node('True',mutGeneComp)
                mutGeneComp.addINode(iNodeGeneMut)
                bkf.addComponent(mutGeneComp)

                # stat condition bin
                statConditionComp = BKB_component("mu-STD>=" + gene + "<=mu+STD=")
                statConditionTrue = BKB_I_node('True',statConditionComp)
                statConditionFalse = BKB_I_node('False',statConditionComp)
                statConditionComp.addINode(statConditionTrue)
                statConditionComp.addINode(statConditionFalse)
                bkf.addComponent(statConditionComp)

                # form SNode  o---->[mut_<genename>=True]
                bkf.addSNode(BKB_S_node(mutGeneComp, iNodeGeneMut, 1.0))

                if True: #add stat condition given data, i.e. compare patient gene reads vs the population mean and std
                    # form SNode  [mut_<genename>=True]----o---->[mu-STD>=<genename><=mu+STD=True
                    bkf.addSNode(BKB_S_node(statConditionComp, statConditionTrue, 1.0, [(mutGeneComp, iNodeGeneMut)]))
                    # form SNode  [mut_<genename>=True]----o---->[mu-STD>=<genename><=mu+STD=False
                    bkf.addSNode(BKB_S_node(statConditionComp, statConditionFalse, 0.0, [(mutGeneComp, iNodeGeneMut)]))
                else:
                    # form SNode  [mut_<genename>=True]----o---->[mu-STD>=<genename><=mu+STD=True
                    bkf.addSNode(BKB_S_node(statConditionComp, statConditionTrue, 0.0, [(mutGeneComp, iNodeGeneMut)]))
                    # form SNode  [mut_<genename>=True]----o---->[mu-STD>=<genename><=mu+STD=False
                    bkf.addSNode(BKB_S_node(statConditionComp, statConditionFalse, 1.0, [(mutGeneComp, iNodeGeneMut)]))
            self.bkfs.append(bkf)

    def BKFsToFile(self, outDirect):
        allBKFHashNames = dict()
        for i in range(0, len(self.bkfs)):
            #i matches the self.patients to self.bkfs.
            hashVal, hashItem = self.BKFHash(i)
            allBKFHashNames[hashVal] = hashItem
            self.bkfs[i].save(outDirect + str(self.bkfs[i].name))
        # write all patient BKF hashs to file
        w = csv.writer(open(outDirect + "patientHashDict.csv", "w"))
        for key, val in allBKFHashNames.items():
            w.writerow([key,val])


    def BKFHash(self, bkfPatientIndex):
        patientName = list()
        #patientID
        patientName.append(("Patient_ID",self.patients[bkfPatientIndex].patientID))
        #patient Cancer type
        patientName.append(("Cancer_Type",self.patients[bkfPatientIndex].cancerType))
        #patient mutated Genes
        patientName.append(("Patient_Genes",((gene) for gene in self.patients[bkfPatientIndex].mutatedGenes)))
        #patient mutated Gene Reads
        patientName.append(("Patient_Gene_Reads",((geneRead) for geneRead in self.patients[bkfPatientIndex].mutatedGeneReads)))
        #patient age of diagnosis
        patientName.append(("Age_of_Diagnosis",self.patients[bkfPatientIndex].ageDiagnos))
        #patient gender
        patientName.append(("Gender",self.patients[bkfPatientIndex].gender))
        #patient survival time
        patientName.append(("Survival_Time",self.patients[bkfPatientIndex].survivalTime))

        patientHashName = hash(tuple(patientName))
        self.bkfs[bkfPatientIndex].name = patientHashName
        return patientHashName, patientName

if __name__ == '__main__':
    PP = PatientProcessor()
    PP.processPatientGeneData('data/wxs.csv', 'data/rnaseq_fpkm_uq_primary_tumor.csv', 'data/geneReadsStats.csv')
    PP.processClinicalData('data/clinical.csv')
    PP.processPatientBKF()
    PP.BKFsToFile('patientBKFs/')
