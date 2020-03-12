import csv
import tqdm
import os
import itertools
import operator
import functools
import math
from pybkb import bayesianKnowledgeBase as BKB
from pybkb import BKB_S_node, BKB_component, BKB_I_node
import matplotlib.pyplot as plt
import pickle
#from concurrent.futures import ThreadPool

class PatientProcessor:
    def __init__(self):
        # self.bkfs[i] and self.patients[i] reference each other
        self.bkfs = list()
        self.patientXBKF = None
        self.geneFrags = None
        self.geneReliabilities = None
        self.geneFragmentsSources = None
        self.patients = list()
        self.holdouts = list()
        self.clinicalCollected = False
        self.radiationCollected = False
        self.radNameOnly = True
        self.drugCollected = False
        self.drugNameOnly = True
        self.sort = True

    class Patient:
        # will need to account for mutatedReads
        def __init__(self, patientID, cancerType, mutatedGenes, mutatedGeneVariants, variants):
            # patient gene data
            self.patientID = patientID
            self.cancerType = cancerType
            self.mutatedGenes = mutatedGenes
            self.mutatedGeneVariants = mutatedGeneVariants
            self.variants = variants
            self.mutatedGeneReads = list()

            # clinical data
            self.ageDiagnos = None
            self.gender = None
            self.pathT = None
            self.pathN = None
            self.pathM = None
            self.survivalTime = None

            # radiation data
            self.therapyName = list()
            self.therapySite = list()
            self.radiationDose = list()
            self.daysToTherapyStart = list()
            self.daysToTherapyEnd = list()

            # add drug info later
            self.drugName = list()
            self.drugDose = list()
            self.daysToDrugStart = list()
            self.daysToDrugEnd = list()

    # processes patient gene data
    def processPatientGeneData(self, patientMutGenes, patientMutReads):
        # reading initial patient data
        with open(patientMutGenes, 'r') as csv_file:
            reader = csv.reader(csv_file)
            patientGeneDict = dict()
            # row[0] = cancer type row[1] = patientID, row[6] = mutated gene
            next(reader)
            rows = [(row[0],row[1],row[6],row[7],row[10],row[13],row[16]) for row in reader]

        # create all patients
        mutatedPatGenes = list()
        mutatedPatGeneVariants = list()
        variants = list()
        cancerType = rows[0][0]
        patID = rows[0][1]
       
        for row in tqdm.tqdm(rows, desc='Reading patient files'):
            variantClassifications = list()
            if patID != row[1]:
                p = self.Patient(patID,cancerType,mutatedPatGenes,mutatedPatGeneVariants,variants)
                self.patients.append(p)
                patID = row[1]
                cancerType = row[0]
                mutatedPatGenes = list()
                mutatedPatGeneVariants = list()
                variants = list()
            if 'DNP' not in row[3]:
                variantClassifications.append(row[3])
            if 'DNP' not in row[4]:
                variantClassifications.append(row[4])
            if 'DNP' not in row[5]:
                variantClassifications.append(row[5])
            if 'DNP' not in row[6]:
                variantClassifications.append(row[6])
            if all(x == variantClassifications[0] for x in variantClassifications):
                mutatedPatGenes.append(row[2])
                variants.append(variantClassifications[0])
                mutatedPatGeneVariants.append(row[2]+"-"+variantClassifications[0])

        csv_file.close()
        '''
        # reading gene reads
        with open(patientMutReads, 'r') as csv_file:
            reader = csv.reader(csv_file)
            geneReadsDict = dict()
            # row[0] = CancerType row[1] = patientID, row[6] = gene, row[7] = reads
            next(reader)
            pbar = tqdm.tqdm(total=os.path.getsize(patientMutReads), desc='Reading patient gene fpkm reads')
            for row in reader:
                geneReadsDict[str(row[0])+str(row[1])+str(row[6])] = float(row[7])
                pbar.update(len(''.join(row).encode('utf-8')))
            pbar.close()

        for p in tqdm.tqdm(self.patients[:], desc='Reading patient mutations'):
            geneMismatch = []
            for gene in p.mutatedGenes[:]:
                if p.cancerType + p.patientID + gene in geneReadsDict:
                    p.mutatedGeneReads.append(geneReadsDict[p.cancerType+p.patientID+gene])
                #gene read information not found
                else:
                    geneMismatch.append(gene)
                    p.mutatedGenes.remove(gene)
            # rare case where we've removed all genes for a patient
            if len(p.mutatedGenes) == 0:
                print("No Reads - Removing: {}".format(p.patientID))
                self.patients.remove(p)
        '''

    def processClinicalData(self, clinicalData):
        assert len(self.patients) > 0, "Patient and gene data have not been read in yet. Call processPatientGeneData(patientMutGenesFile, patientMutReadsFile) first, before secondary processing functions are called."
        with open(clinicalData, 'r') as csv_file:
            reader = csv.reader(csv_file)
            patientClinicalDict = dict()
            # row[0] = cancerType, row[1] = patientID, row[2] = age of diagnosis, row[3] = gender, row[4] = Pathologic_T, row[5] = Pathologic_N, row[6] = Pathologic_M , row[7] = survival time
            pbar = tqdm.tqdm(total=os.path.getsize(clinicalData), desc='Reading clinical data')
            next(reader)
            for row in reader:
                if 'DNP' in row[2]:
                    continue
                ageDiag = int(float(row[2]))
                if 'DNP' in row[3]:
                    continue
                gender = row[3]
                if 'DNP' in row[4]:
                    continue
                pathT = row[4]
                if 'DNP' in row[5]:
                    continue
                pathN = row[5]
                if 'DNP' in row[6]:
                    continue
                pathM = row[6]
                if 'DNP' in row[7]:
                    continue
                survTime = int(float(row[7]))
                patientClinicalDict[str(row[0])+str(row[1])] = (ageDiag,gender,pathT,pathN,pathM,survTime)
                pbar.update(len(''.join(row).encode('utf-8')))
            pbar.close()
        for p in self.patients[:]:
            if str(p.cancerType)+str(p.patientID) in patientClinicalDict:
                clinical = patientClinicalDict[str(p.cancerType)+str(p.patientID)]
                p.ageDiagnos = clinical[0]
                p.gender = clinical[1]
                p.pathT = clinical[2]
                p.pathN = clinical[3]
                p.pathM = clinical[4]
                p.survivalTime = clinical[5]
            else:
                print("No Clinical - Removing: {}".format(p.patientID))
                self.patients.remove(p)
        self.clinicalCollected = True

    def processRadiationData(self, radiationData, radNameOnly=True):
        assert len(self.patients) > 0, "Patient and gene data have not been read in yet. Call processPatientGeneData(patientMutGenesFile, patientMutReadsFile) first, before secondary processing functions are called."
        self.radNameOnly = radNameOnly
        with open(radiationData, 'r') as csv_file:
            reader = csv.reader(csv_file)
            patientRadiationDict = dict()
            # row[0] = cancerType, row[1] = patientID, row[2] = Therapy name, row[3] = therapy Site, row[4] = dose, row[5] = doseUnit, row[6] = days to start, row[7] = days to end
            pbar = tqdm.tqdm(total=os.path.getsize(radiationData), desc='Reading radiation data')
            next(reader)
            for row in reader:
                if radNameOnly:
                    cancerType = row[0]
                    patID = row[1]
                    if 'UNKNOWN' in row[2] or 'Discrepancy' in row[2] or 'DNP' in row[2]:
                        continue
                    therapyName = row[2]
                    if cancerType + patID in patientRadiationDict:
                        patientRadiationDict[cancerType + patID].append(therapyName)
                    else:
                        patientRadiationDict[cancerType + patID] = [therapyName]
                    pbar.update(len(''.join(row).encode('utf-8')))
                else:
                    cancerType = row[0]
                    patID = row[1]
                    if 'UNKNOWN' in row[2] or 'Discrepancy' in row[2] or 'DNP' in row[2]:
                        continue
                    therapyName = row[2]
                    if 'UNKNOWN' in row[3] or 'Discrepancy' in row[3] or 'DNP' in row[3]:
                        continue
                    therapySite = row[3]
                    if 'UNKNOWN' in row[4] or 'Discrepancy' in row[4] or 'DNP' in row[4]:
                        continue
                    dose = None
                    try:
                        if row[5] == 'Gy':
                            dose = int(float(row[4])) * 100 #conversion rate between cGy and Gy. Gy = 100 cGy 
                        elif row[5] == 'cGy':
                            dose = int(float(row[4]))
                        else:
                            continue
                    except:
                        continue
                    if 'UNKNOWN' in row[6] or 'Discrepancy' in row[6] or 'DNP' in row[6]:
                        continue
                    daysToStart = int(float(row[6]))
                    if 'UNKNOWN' in row[7] or 'Discrepancy' in row[7] or 'DNP' in row[7]:
                        continue
                    daysToEnd = int(float(row[7]))
                    if cancerType + patID in patientRadiationDict:
                        patientRadiationDict[cancerType + patID].append((therapyName, therapySite, dose, daysToStart, daysToEnd))
                    else:
                        patientRadiationDict[cancerType + patID] = [(therapyName, therapySite, dose, daysToStart, daysToEnd)]
                pbar.update(len(''.join(row).encode('utf-8')))
            pbar.close()
        for p in self.patients[:]:
            if p.cancerType + p.patientID in patientRadiationDict:
                radiationData = patientRadiationDict[p.cancerType+p.patientID]
                for rData in radiationData:
                    if radNameOnly:
                        p.therapyName.append(rData)
                    else:
                        p.therapyName.append(rData[0])
                        p.therapySite.append(rData[1])
                        p.radiationDose.append(rData[2])
                        p.daysToTherapyStart.append(rData[3])
                        p.daysToTherapyEnd.append(rData[4])
            else:
                print("No Radiation - Removing: {}".format(p.patientID))
                self.patients.remove(p) 
        self.radiationCollected = True

    def processDrugData(self, drugData, drugNameOnly=True):
        assert len(self.patients) > 0, "Patient and gene data have not been read in yet. Call processPatientGeneData(patientMutGenesFile, patientMutReadsFile) first, before secondary processing functions are called."
        self.drugNameOnly = drugNameOnly
        with open(drugData, 'r') as csv_file:
            reader = csv.reader(csv_file)
            patientDrugDict = dict()
            # row[0] = cancerType, row[1] = patientID, row[2] = Drug name, row[3] = dose, row[4] = doseUnit, row[5] = days to start, row[6] = days to end
            pbar = tqdm.tqdm(total=os.path.getsize(drugData), desc='Reading drug data')
            next(reader)
            for row in reader:
                if drugNameOnly:
                    cancerType = row[0]
                    patID = row[1]
                    if 'UNKNOWN' in row[2] or 'Discrepancy' in row[2] or 'DNP' in row[2]:
                        continue
                    drugName = row[2]
                    if cancerType + patID in patientDrugDict:
                        patientDrugDict[cancerType + patID].append(drugName)
                    else:
                        patientDrugDict[cancerType + patID] = [drugName]
                else:
                    cancerType = row[0]
                    patID = row[1]
                    if 'UNKNOWN' in row[2] or 'Discrepancy' in row[2] or 'DNP' in row[2]:
                        continue
                    drugName = row[2]
                    dose = None
                    try:
                        if row[4] == 'mg':
                            dose = float(row[3]) #conversion rate between cGy and Gy. Gy = 100 cGy
                        elif row[4] == 'g':
                            dose = float(row[3])/1000.0
                        elif row[4] == 'ug':
                            dose = float(row[3])*1000.0
                        else:
                            continue
                    except:
                        continue
                    if 'UNKNOWN' in row[5] or 'Discrepancy' in row[5] or 'DNP' in row[5]:
                        continue
                    daysToStart = int(float(row[5]))
                    if 'UNKNOWN' in row[6] or 'Discrepancy' in row[6] or 'DNP' in row[6]:
                        continue
                    daysToEnd = int(float(row[6]))
                    if cancerType + patID in patientDrugDict:
                        patientDrugDict[cancerType + patID].append((drugName, dose, daysToStart, daysToEnd))
                    else:
                        patientDrugDict[cancerType + patID] = [(drugName, dose, daysToStart, daysToEnd)]
                pbar.update(len(''.join(row).encode('utf-8')))
            pbar.close()
        for p in self.patients[:]:
            if p.cancerType + p.patientID in patientDrugDict:
                drugData = patientDrugDict[p.cancerType+p.patientID]
                for dData in drugData:
                    if drugNameOnly:
                        p.drugName.append(dData)
                    else:
                        p.drugName.append(dData[0])
                        p.drugDose.append(dData[1])
                        p.daysToDrugStart.append(dData[2])
                        p.daysToDrugEnd.append(dData[3])
            else:
                rint("No Drug - Removing: {}".format(p.patientID))
                self.patients.remove(p)
        self.drugCollected = True

    def setAsideHoldouts(self, patients, holdouts):
        #self.patients should = patients at end. These will be placed in newPatients
        #self.holdouts = holdouts at end. these will be placed in newHoldouts
        #everything else is removed
        newPatients = list()
        newHoldouts = list()
        holdoutIDs = [ho[1] for ho in holdouts]
        for pat in self.patients:
            if pat.patientID in holdoutIDs:
                newHoldouts.append(pat)
            elif pat in patients:
                newPatients.append(pat)
            else:
                print("removing:",pat.patientID)
        self.patients = newPatients
        self.holdouts = newHoldouts

    #-- Default interaction is pair-wise.
    def getGeneInteractions(self, interaction=2):
        #-- Construct Patient Data dict
        patient_data = dict()
        pop_genes = set()
        counts_B = dict()
        counts_AB = dict()
        for i in tqdm.tqdm(range(len(self.patients)), desc='Processing patients', leave=False):
            patientHash, patientDict = self.BKFHash(i,True)
            patient_data[patientHash] = patientDict
            pat_genes = patientDict['Patient_Genes']
            pop_genes.update(set(pat_genes))
            total_combos = math.factorial(len(pat_genes)) / (math.factorial(interaction-1) * math.factorial(len(pat_genes) - (interaction-1)))
            for A in tqdm.tqdm(pat_genes, desc='Processing interactions', leave=False):
                #if A in gene_counts:
                #    gene_counts[A] += 1
                #else:
                #    gene_counts[A] = 1
                for B in itertools.combinations(set(pat_genes)-{A}, interaction-1):
                    AB = tuple([A] + list(B))
                    if B in counts_B:
                        counts_B[B] += 1
                    else:
                        counts_B[B] = 1
                    if AB in counts_AB:
                        counts_AB[AB] += 1
                    else:
                        counts_AB[AB] = 1
        probs = dict()
        for AB, count_AB in tqdm.tqdm(counts_AB.items(), desc='Calculating probabilities', leave=False):
            A = AB[0]
            B = AB[1:]
            if A not in probs:
                probs[A] = {B: count_AB / counts_B[B]}
            else:
                probs[A][B] = count_AB / counts_B[B]
        return probs

        '''
        #-- Calculate gene conditional probabilities
        probs = {gene: dict() for gene in pop_genes}
        total_combos = math.factorial(len(pop_genes)) / (math.factorial(interaction-1) * math.factorial(len(pop_genes) - (interaction-1)))
        for A in tqdm.tqdm(pop_genes, desc='Calculating Interaction Probabilities', leave=False):
            for B in tqdm.tqdm(itertools.combinations(set(pop_genes) - {A}, interaction - 1), total=total_combos):
                AB_set = set(list(B) + [A])
                #-- Calculating P(A|B) = P(A,B) / P(B)
                count_AB = 0
                count_B = 0
                for _, data_dict in patient_data.items():
                    if len(AB_set - set(data_dict['Patient_Genes'])) == 0:
                        count_AB += 1
                    if len(set(B) - set(data_dict['Patient_Genes'])) == 0:
                        count_B += 1
                prob_AB = float(count_AB / len(patient_data))
                prob_B = float(count_B / len(patient_data))
                try:
                    prob_A_B = prob_AB / prob_B
                    probs[A][B] = prob_A_B
                except ZeroDivisionError:
                    continue
        return probs
        '''

    # should be called after all Patient objects have been cosntructed
    def processPatientBKF(self, patientFalses=True):
        assert len(self.patients) > 0, "Have not processed Patient and gene data yet."
        # collect list of all unqiue genes
        all_gene_mutations = list()
        for pat in self.patients:
            for gene in pat.mutatedGenes:
                if gene not in all_gene_mutations:
                    all_gene_mutations.append(gene)

        # get gene causeality dictionary
        geneCauses = self.getGeneInteractions()

        geneFreqNotIn = None
        geneFreqHashSource = None
        if not patientFalses:
            geneFreqNotIn = dict()
            geneFreqHashSource = dict()
            for gene in all_gene_mutations:
                geneFreqNotIn[gene] = 0
                geneFreqHashSource[gene] = []


        # make all patient BKFs
        for pat in tqdm.tqdm(self.patients, desc='Forming patient BKFs'):
            bkf = BKB(name = pat.patientID)
            for aGene in all_gene_mutations:
                mutGeneComp_idx = bkf.addComponent('mut_{}'.format(aGene))
                if aGene in pat.mutatedGenes:
                    iNodeGeneMut_idx = bkf.addComponentState(mutGeneComp_idx, 'True')
                    bkf.addSNode(BKB_S_node(mutGeneComp_idx, iNodeGeneMut_idx, 1.0))

                    _mutGeneComp_idx = bkf.addComponent('_mut_{}'.format(aGene))
                    _iNodeGeneMut_idx = bkf.addComponentState(_mutGeneComp_idx, 'True')

                    mutVarComp_idx = bkf.addComponent('mut-var_{}'.format(aGene))
                    iNodeMutVar_idx = bkf.addComponentState(mutVarComp_idx, pat.variants[pat.mutatedGenes.index(aGene)])

                    bkf.addSNode(BKB_S_node(mutVarComp_idx, iNodeMutVar_idx, 1.0, [(_mutGeneComp_idx, _iNodeGeneMut_idx)]))

                else:
                    if geneFreqNotIn is not None:
                        geneFreqNotIn[aGene] += 1
                        geneFreqHashSource[aGene].append(hash(pat.patientID))
                    else:
                        iNodeGeneMut_idx = bkf.addComponentState(mutGeneComp_idx, 'False')
                        bkf.addSNode(BKB_S_node(mutGeneComp_idx, iNodeGeneMut_idx, 1.0))
            self.bkfs.append(bkf)

        #add individual gene false fragments
        if not patientFalses:
            self.geneFrags = list()
            self.geneReliabilities = list()
            self.geneFragmentsSources = list()
            for aGene in all_gene_mutations:
                genebkf = BKB(name = aGene)
                mutGeneComp_idx = genebkf.addComponent('mut_{}'.format(aGene))
                iNodeGeneMut_idx = genebkf.addComponentState(mutGeneComp_idx, 'False')
                genebkf.addSNode(BKB_S_node(mutGeneComp_idx, iNodeGeneMut_idx, 1.0))
                self.geneFrags.append(genebkf)
                self.geneReliabilities.append(geneFreqNotIn[aGene])
                source = None
                if self.sort:
                    source = ','.join([str(hashName) for hashName in sorted(geneFreqHashSource[aGene])]) 
                else:
                    source = ','.join([str(hashName) for hashName in geneFreqHashSource[aGene]])
                self.geneFragmentsSources.append(source)

        self.patientXBKF = BKB(name = "patientX")
        for aGene in all_gene_mutations:
            mutGeneComp_idx = self.patientXBKF.addComponent('mut_{}'.format(aGene))
            iNodeGeneMutTrue_idx = self.patientXBKF.addComponentState(mutGeneComp_idx, 'True')
            iNodeGeneMutFalse_idx = self.patientXBKF.addComponentState(mutGeneComp_idx, 'False')

            #get causality gene
            causes = geneCauses[aGene].keys()
            vals = geneCauses[aGene].values()
            maxCauseGene = None
            maxCauseVal = 0
            for bGene, val in geneCauses[aGene].items():
                if val > maxCauseVal and bGene[0] != aGene:
                    maxCauseVal = val
                    maxCauseGene = bGene[0]

            cMutGeneComp_idx = self.patientXBKF.addComponent('mut_{}'.format(maxCauseGene))
            iNodeCGeneMut_idx = self.patientXBKF.addComponentState(cMutGeneComp_idx, 'True')

            _mutGeneComp_idx = self.patientXBKF.addComponent('_mut_{}'.format(aGene))
            _iNodeGeneMut_idx = self.patientXBKF.addComponentState(_mutGeneComp_idx, 'True')

            self.patientXBKF.addSNode(BKB_S_node(_mutGeneComp_idx, _iNodeGeneMut_idx, 1.0, [(mutGeneComp_idx, iNodeGeneMutTrue_idx)]))
            self.patientXBKF.addSNode(BKB_S_node(_mutGeneComp_idx, _iNodeGeneMut_idx, maxCauseVal, [(mutGeneComp_idx, iNodeGeneMutFalse_idx),
                                                                                                    (cMutGeneComp_idx, iNodeCGeneMut_idx)]))

    def SubsetBKFsToFile(self, outDirect, indices):
        #return bkf_files, source_names, patient_data_file
        allBKFHashNames = dict()
        bkf_files = list()
        source_names = list()
        for idx in indices:
            patientHashVal, patientDict = self.BKFHash(idx, True)
            allBKFHashNames[patientHashVal] = patientDict
            bkf_files.append(outDirect + str(patientHashVal) + '.bkf')
            self.bkfs[idx].save(outDirect + str(patientHashVal) + '.bkf')
            source_names.append(str(patientHashVal))
        holdout_source_names = list()
        for i in range(0, len(self.holdouts)):
            patientHashVal, patientDict = self.BKFHash(i,False)
            allBKFHashNames[patientHashVal] = patientDict
            holdout_source_names.append(str(patientHashVal))
        # write all patient BKF hashs to file
        patient_data_file = outDirect + 'patient_data.pk'
        with open(patient_data_file, 'wb') as f_:
            pickle.dump(file=f_, obj=allBKFHashNames)
        return bkf_files, source_names, holdout_source_names, patient_data_file


    def BKFsToFile(self, outDirect):
        print("Writing patient BKFs to file...")
        allBKFHashNames = dict()
        for i in range(0, len(self.bkfs)):
            #i matches the self.patients to self.bkfs.
            patientHashVal, patientDict = self.BKFHash(i,True)
            allBKFHashNames[patientHashVal] = patientDict
            self.bkfs[i].save(outDirect + str(self.bkfs[i].getName()) + '.bkf')
        for i in range(0, len(self.holdouts)):
            patientHashVal, patientDict = self.BKFHash(i,False)
            allBKFHashNames[patientHashVal] = patientDict
        # write all patient BKF hashs to file
        with open(outDirect + 'patient_data.pk', 'wb') as f_:
            pickle.dump(file=f_, obj=allBKFHashNames)
        print("Patient BKFs written.")


    def BKFHash(self,index, patient):
        pat = None
        if patient:
            pat = self.patients[index]
        else:
            pat = self.holdouts[index]
        patientDict = dict()
        #patientID
        patientDict["Patient_ID"] = pat.patientID
        #patient Cancer type
        patientDict["Cancer_Type"] = pat.cancerType
        #patient mutated Genes
        patientDict["Patient_Genes"] = tuple(pat.mutatedGenes)
        #patient mutated Gene Variants
        patientDict["Patient_Gene_Variants"] = tuple(pat.mutatedGeneVariants)
        #patient variants
        patientDict["Patient_Variants"] = tuple(pat.variants)
        #patient mutated Gene Reads
        patientDict["Patient_Gene_Reads"] = tuple(pat.mutatedGeneReads)
        if self.clinicalCollected:
            #patient age of diagnosis
            patientDict["Age_of_Diagnosis"] = pat.ageDiagnos
            #patient gender
            patientDict["Gender"] = pat.gender
            #pathT
            patientDict["PathT"] = pat.pathT
            #pathN
            patientDict["PathN"] = pat.pathN
            #pathM
            patientDict["PathM"] = pat.pathM
            #patient survival time
            patientDict["Survival_Time"] = pat.survivalTime
        #patients can have no radiation data associated with them --------------------v
        if self.radiationCollected and len(pat.therapyName) > 0:
            if self.radNameOnly:
                patientDict["Therapy_Name(s)"] = tuple(pat.therapyName)
            else:
                # patient therapy name(s)
                patientDict["Therapy_Name(s)"] = tuple(pat.therapyName)
                # patient therapy site(s)
                patientDict["Therapy_Site(s)"] = tuple(pat.therapySite)
                # patient radiation dose(s)
                patientDict["Radiation_Dose(s)"] = tuple(pat.radiationDose)
                # patient days to therapy start
                patientDict["Days_To_Therapy_Start(s)"] = tuple(pat.daysToTherapyStart)
                # patient days to therapy end
                patientDict["Days_To_Therapy_End(s)"] = tuple(pat.daysToTherapyEnd)
        if self.drugCollected and len(pat.drugName) > 0:
            if self.drugNameOnly:
                patientDict["Drug_Name(s)"] = tuple(pat.drugName)
            else:
                # patient drug name(s)
                patientDict["Drug_Name(s)"] = tuple(pat.drugName)
                # patient drug dose(s)
                patientDict["Drug_Dose(s)"] = tuple(pat.drugDose)
                # patient days to drug start
                patientDict["Days_To_Drug_Start(s)"] = tuple(pat.daysToDrugStart)
                # patient days to drug end
                patientDict["Days_To_Drug_End(s)"] = tuple(pat.daysToDrugEnd)

        patientHashVal = hash(pat.patientID)
        return patientHashVal, patientDict

if __name__ == '__main__':
    PP = PatientProcessor()
    PP.processPatientGeneData('/home/public/data/ncats/data_drop_02-11-2020/wxs.csv',
                              '/home/public/data/ncats/data_drop_02-11-2020/rnaseq_fpkm_uq_primary_tumor.csv')
    #PP.processClinicalData('data/clinical.csv')
    PP.processPatientBKF()
    #PP.BKFsToFile('BKFs/')
    probs = PP.getGeneInteractions()
    print(probs)
