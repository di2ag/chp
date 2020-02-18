import csv
from pybkb import bayesianKnowledgeBase as BKB
from pybkb import BKB_S_node, BKB_component, BKB_I_node
import matplotlib.pyplot as plt
import pickle

class PatientProcessor:
    def __init__(self):
        # self.bkfs[i] and self.patients[i] reference each other
        self.bkfs = list()
        self.patients = list()
        self.clinicalCollected = False

    class Patient:
        # will need to account for mutatedReads
        def __init__(self, patientID, cancerType, mutatedGenes):
            # patient gene data
            self.patientID = patientID
            self.cancerType = cancerType
            self.mutatedGenes = mutatedGenes
            self.mutatedGeneReads = list()

            # clinical later
            self.ageDiagnos = None
            self.gender = None
            self.survivalTime = None

            # add drug info later
            # add radiation info later

    # processes patient gene data
    def processPatientGeneData(self, patientMutGenes, patientMutReads, geneReadStats):
        print("Reading patient files...")
        with open(patientMutGenes, 'r') as csv_file:
            reader = csv.reader(csv_file)
            patientGeneDict = dict()
            # row[0] = cancer type row[1] = patientID, row[6] = mutated gene
            next(reader)
            rows = [(row[0],row[1],row[6]) for row in reader]

        # create all patients
        mutatedPatGenes = list()
        cancerType = rows[0][0]
        patID = rows[0][1]
        for row in rows:
            if patID != row[1]:
                p = self.Patient(patID,cancerType,mutatedPatGenes)
                self.patients.append(p)
                patID = row[1]
                cancerType = row[0]
                mutatedPatGenes = list()
            mutatedPatGenes.append(row[2])
        csv_file.close()
        print("Patient files read.")

        print("Reading patient gene fpkm reads...")
        # get gene reads from tumor tissue

        with open(patientMutReads, 'r') as csv_file:
            reader = csv.reader(csv_file)
            geneReadsDict = dict()
            # row[1] = patientID, row[6] = gene, row[7] = reads
            next(reader)
            for row in reader:
                geneReadsDict[str(row[1])+str(row[6])] = float(row[7])
        for p in self.patients:
            geneMismatch = []
            for gene in p.mutatedGenes[:]:
                if str(p.patientID)+str(gene) in geneReadsDict:
                    p.mutatedGeneReads.append(geneReadsDict[p.patientID+gene])
                #gene read information not found
                else:
                    geneMismatch.append(gene)
                    p.mutatedGenes.remove(gene)
            if len(geneMismatch) > 0:
                print("Warning - Patient", p.patientID, "has", len(geneMismatch), "missing genes")
        print("Patient gene fpkm reads read.")

    def processClinicalData(self, clinicalData):
        assert len(self.patients) > 0, "Patient and gene data have not been read in yet. Call processPatientGeneData(patientMutGenesFile, patientMutReadsFile) first, before secondary processing functions are called."
        print("Reading patient clinical data...")
        with open(clinicalData, 'r') as csv_file:
            reader = csv.reader(csv_file)
            patientClinicalDict = dict()
            # row[0] = cancerType, row[1] = patientID, row[2] = age of diagnosis, row[3] = gender, row[7] = survival time
            next(reader)
            for row in reader:
                ageDiag = int(float(row[2])) if 'DNP' not in row[2] else row[2]
                gender = row[3]
                survTime = int(float(row[7])) if 'DNP' not in row[7] else row[7]
                patientClinicalDict[str(row[0])+str(row[1])] = (ageDiag,gender,survTime)

        for p in self.patients:
            if str(p.cancerType)+str(p.patientID) in patientClinicalDict:
                clinical = patientClinicalDict[str(p.cancerType)+str(p.patientID)]
                p.ageDiagnos = clinical[0]
                p.gender = clinical[1]
                p.survivalTime = clinical[2]
            else:
                print("WARNING - Patient", p, "clinical data does not exist")
        self.clinicalCollected = True
        print("Patient clinical data read.")

    # should be called after all Patient objects have been cosntructed
    def processPatientBKF(self):
        assert len(self.patients) > 0, "Have not processed Patient and gene data yet."
        print("Forming Patient BKFs...")

        for pat in self.patients:
            bkf = BKB(name = pat.patientID)
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
        print("Patient BKFs formed.")

    def BKFsToFile(self, outDirect):
        print("Writing patient BKFs to file...")
        allBKFHashNames = dict()
        for i in range(0, len(self.bkfs)):
            #i matches the self.patients to self.bkfs.
            patientHashVal, patientDict = self.BKFHash(i)
            allBKFHashNames[patientHashVal] = patientDict
            self.bkfs[i].save(outDirect + str(self.bkfs[i].name))
        # write all patient BKF hashs to file
        with open(outDirect + 'patient_data.pk', 'wb') as f_:
            pickle.dump(file=f_, obj=allBKFHashNames)
        print("Patient BKFs written.")


    def BKFHash(self, bkfPatientIndex):
        patientDict = dict()
        #patientID
        patientDict["Patient_ID"] = self.patients[bkfPatientIndex].patientID
        #patient Cancer type
        patientDict["Cancer_Type"] = self.patients[bkfPatientIndex].cancerType
        #patient mutated Genes
        patientDict["Patient_Genes"] = tuple(self.patients[bkfPatientIndex].mutatedGenes)
        #patient mutated Gene Reads
        patientDict["Patient_Gene_Reads"] = tuple(self.patients[bkfPatientIndex].mutatedGeneReads)
        #patient age of diagnosis
        patientDict["Age_of_Diagnosis"] = self.patients[bkfPatientIndex].ageDiagnos
        #patient gender
        patientDict["Gender"] = self.patients[bkfPatientIndex].gender
        #patient survival time
        patientDict["Survival_Time"] = self.patients[bkfPatientIndex].survivalTime

        patientHashVal = hash(self.patients[bkfPatientIndex].patientID)
        self.bkfs[bkfPatientIndex].name = patientHashVal
        return patientHashVal, patientDict

if __name__ == '__main__':
    PP = PatientProcessor()
    PP.processPatientGeneData('data/wxs.csv', 'data/rnaseq_fpkm_uq_primary_tumor.csv', 'data/geneReadsStats.csv')
    PP.processClinicalData('data/clinical.csv')
    PP.processPatientBKF()
    PP.BKFsToFile('patientBKFs/')
