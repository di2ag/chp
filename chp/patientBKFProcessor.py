import csv
import tqdm
import os
import itertools
import operator
import functools
import math
import random
from pybkb.pybkb.common.bayesianKnowledgeBase import bayesianKnowledgeBase as BKB
from pybkb.pybkb.common.bayesianKnowledgeBase import BKB_S_node, BKB_component, BKB_I_node
import matplotlib.pyplot as plt
import pickle
import pandas as pd
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
            self.processType = list()
            self.processActivity = list()
            self.biologicalObject = list()
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
                cancerType = row[0]
                patID = row[1]
                if 'UNKNOWN' in row[2] or 'Discrepancy' in row[2] or 'DNP' in row[2]:
                    continue
                drugName = row[2]
                #row7 = ProcessType
                processType = row[7]
                processActivity = row[8]
                biologicalObject = row[9]
                if cancerType + patID in patientDrugDict:
                    if (drugName, processType, processActivity, biologicalObject) not in patientDrugDict[cancerType + patID]:
                        patientDrugDict[cancerType + patID].append((drugName, processType, processActivity, biologicalObject))
                else:
                    patientDrugDict[cancerType + patID] = [(drugName, processType, processActivity, biologicalObject)]
                pbar.update(len(''.join(row).encode('utf-8')))
            pbar.close()
        for p in self.patients[:]:
            if p.cancerType + p.patientID in patientDrugDict:
                drugData = patientDrugDict[p.cancerType+p.patientID]
                for dData in drugData:
                    p.drugName.append(dData[0])
                    p.processType.append(dData[1])
                    p.processActivity.append(dData[2])
                    p.biologicalObject.append(dData[3])
            else:
                print("No Drug - Removing: {}".format(p.patientID))
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

    def getGeneTFProbs(self, allGenes):
        FT_Dict = dict()
        pop_genes = set(allGenes)
        for pat in tqdm.tqdm(self.patients, desc='Processing patient TF counts', leave=False):
            pat_genes = set(pat.mutatedGenes)
            not_pat_genes = pop_genes - pat_genes
            for npg in not_pat_genes:
                if str(npg) not in FT_Dict.keys():
                    FT_Dict[str(npg)] = dict()
                for pg in pat_genes:
                    if str(pg) not in FT_Dict[str(npg)].keys():
                        FT_Dict[str(npg)][str(pg)] = 1
                    else:
                        FT_Dict[str(npg)][str(pg)] += 1
        for geneA in tqdm.tqdm(allGenes, desc='Processing patient probabilities', leave=False):
            for geneB in allGenes:
                if str(geneA) not in FT_Dict.keys() or str(geneB) not in FT_Dict[str(geneA)].keys():
                    continue
                FT_Dict[str(geneA)][str(geneB)] = FT_Dict[str(geneA)][str(geneB)]/len(self.patients)
        return FT_Dict

    def getPatientDataCorrelations(self, outDirect=None):
        assert len(self.patients) > 0, "Patient data has not been collected yet!"
        # get all column names:
        column_names = list()
        column_names.append("Survival_Time:")
        column_names.append("Age_of_Diagnosis:")

        for pat in tqdm.tqdm(self.patients, desc = 'One-hot encoding feature column headers'):
            for x in pat.variants:
                if "Variants:" + str(x) not in column_names:
                    column_names.append("Variants:" + str(x))
            for x in pat.mutatedGenes:
                if "Genes:" + str(x) not in column_names:
                    column_names.append("Genes:" + str(x))
            for x in pat.mutatedGeneVariants:
                if "Gene_var:" + str(x) not in column_names:
                    column_names.append("Gene_var:" + str(x))

            # only if drugs were collected
            if self.drugCollected:
                for x in pat.drugName:
                    if "Drug_Name:" + str(x) not in column_names:
                        column_names.append("Drug_Name:" + str(x))
                for x in pat.processType:
                    if "Process_Type:" + str(x) not in column_names:
                        column_names.append("Process_Type:" + str(x))
                for x in pat.processActivity:
                    if "Process_Activity:" + str(x) not in column_names:
                        column_names.append("Process_Activity:" + str(x))
                for x in pat.biologicalObject:
                    if "Biological_Object:" + str(x) not in column_names:
                        column_names.append("Biological_Object:" + str(x))
            # only if clinical was collected
            if self.clinicalCollected:
                if "Path_T:" + str(pat.pathT) not in column_names:
                    column_names.append("Path_T:" + str(pat.pathT))
                if "Path_N:" + str(pat.pathN) not in column_names:
                    column_names.append("Path_N:" + str(pat.pathN))
                if "Path_M:" + str(pat.pathM) not in column_names:
                    column_names.append("Path_M:" + str(pat.pathM))

        # form data lists and populate for each patient
        data = list()

        for pat in tqdm.tqdm(self.patients, desc='Populating one-hot encoded feature table'):
            pat_data = list()
            for cn in column_names:
                header = cn.split(':')
                header_type = header[0]
                header_val = header[1]

                if "Survival_Time" == header_type:
                    pat_data.append(pat.survivalTime)
                elif "Age_of_Diagnosis" == header_type:
                    pat_data.append(pat.ageDiagnos)
                elif "Drug_Name" == header_type:
                    pat_data.append(1 if header_val in pat.drugName else 0)
                elif "Process_Type" == header_type:
                    pat_data.append(1 if header_val in pat.processType else 0)
                elif "Process_Activity" == header_type:
                    pat_data.append(1 if header_val in pat.processActivity else 0)
                elif "Biological_Object" == header_type:
                    pat_data.append(1 if header_val in pat.biologicalObject else 0)
                elif "Variants" == header_type:
                    pat_data.append(1 if header_val in pat.variants else 0)
                elif "Genes" == header_type:
                    pat_data.append(1 if header_val in pat.mutatedGenes else 0)
                elif "Gene_var" == header_type:
                    pat_data.append(1 if header_val in pat.mutatedGeneVariants else 0)
                elif "Path_T" == header_type:
                    pat_data.append(1 if header_val == pat.pathT else 0)
                elif "Path_N" == header_type:
                    pat_data.append(1 if header_val == pat.pathN else 0)
                elif "Path_M" == header_type:
                    pat_data.append(1 if header_val == pat.pathM else 0)
                else:
                    print(header)
            data.append(pat_data)

        # form dataframe and calculate correlations to determine good predictors
        print("Running pandas dataframe coorelations")
        df = pd.DataFrame(data, columns=column_names)
        corrs = df.corr(method='pearson')
        if outDirect is not None:
            df.to_csv(outDirect + 'data_table.csv')
            corrs.to_csv(outDirect + 'data_corrs.csv')
        else:
            print('Saving coorelation data to /tmp')
            df.to_csv('/tmp/data_table.csv')
            df.to_csv('/tmp/data_corrs.csv')

    def processAxlePatientBKF(self):
        assert len(self.patients) > 0, "Have not processed Patient and gene data yet."
        #print("Forming Patient BKFs...")

        for pat in tqdm.tqdm(self.patients, desc='Forming patient BKFs'):
            bkf = BKB(name = pat.patientID)
            for idx, gene in enumerate(pat.mutatedGenes):
                # gene
                mutGeneComp_idx = bkf.addComponent('mut_{}'.format(gene))
                iNodeGeneMut_idx = bkf.addComponentState(mutGeneComp_idx, 'True')
                # form SNode  o---->[mut_<genename>=True]
                bkf.addSNode(BKB_S_node(mutGeneComp_idx, iNodeGeneMut_idx, 1.0))

                # var
                varComp_idx = bkf.addComponent('var_{}'.format(pat.variants[idx]))
                iNodeVar_idx = bkf.addComponentState(varComp_idx, 'True')
                bkf.addSNode(BKB_S_node(varComp_idx, iNodeVar_idx, 1.0))

                #mut_var combo
                mutVarComp_idx = bkf.addComponent('mut-var_{}'.format(gene+"-"+pat.variants[idx]))
                iNodeMutVar_idx = bkf.addComponentState(mutVarComp_idx, pat.variants[idx])
                # form SNode [mut_<genename>=True]---->o---->[mut-var_<genename>=<varianttype>]
                bkf.addSNode(BKB_S_node(mutVarComp_idx, iNodeMutVar_idx, 1.0, [(mutGeneComp_idx, iNodeGeneMut_idx),
                                                                                (varComp_idx, iNodeVar_idx)]))
            self.bkfs.append(bkf)

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
        geneCauses = self.getGeneTFProbs(all_gene_mutations)

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

        #read in topological sorting freq file
        topological = dict()
        with open('/home/public/data/ncats/data_drop_03-04-2020/gene_freq_in_wxs.csv', 'r') as csv_file:
            reader = csv.reader(csv_file)
            next(reader)
            count = 0
            for i, row in enumerate(reader):
                topological[str(row[0])] = i
                topological_genes = str(row[0])
                count += 1

        self.patientXBKF = BKB(name = "patientX")
        geneOrderings = dict()
        orderFrom = list()
        for aGene in all_gene_mutations:
            mutGeneComp_idx = self.patientXBKF.addComponent('mut_{}'.format(aGene))
            iNodeGeneMutTrue_idx = self.patientXBKF.addComponentState(mutGeneComp_idx, 'True')
            iNodeGeneMutFalse_idx = self.patientXBKF.addComponentState(mutGeneComp_idx, 'False')

            #get causality gene
            maxCauseGene = None
            maxCauseVal = 0
            matches = list()
            for bGene, val in geneCauses[aGene].items():
                if val > maxCauseVal and bGene != aGene and topological[bGene] < topological[aGene]:
                    maxCauseVal = val
                    maxCauseGene = bGene
            if maxCauseGene is None:
                print("no gene link for",aGene)

                _mutGeneComp_idx = self.patientXBKF.addComponent('_mut_{}'.format(aGene))
                _iNodeGeneMut_idx = self.patientXBKF.addComponentState(_mutGeneComp_idx, 'True')

                self.patientXBKF.addSNode(BKB_S_node(_mutGeneComp_idx, _iNodeGeneMut_idx, 1.0, [(mutGeneComp_idx, iNodeGeneMutTrue_idx)]))
                #self.patientXBKF.addSNode(BKB_S_node(_mutGeneComp_idx, _iNodeGeneMut_idx, maxCauseVal, [(mutGeneComp_idx, iNodeGeneMutFalse_idx)]))

            else:

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
            # patient drug name(s)
            patientDict["Drug_Name(s)"] = tuple(pat.drugName)
            # patient drug dose(s)
            patientDict["Biological_Object(s)"] = tuple(pat.biologicalObject)
            # patient days to drug start
            patientDict["Process_Activity(s)"] = tuple(pat.processActivity)
            # patient days to drug end
            patientDict["Process_Type(s)"] = tuple(pat.processType)

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
