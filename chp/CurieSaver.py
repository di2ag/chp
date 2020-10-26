from biothings_explorer.hint import Hint
from biothings_explorer.user_query_dispatcher import FindConnection
import os
import pickle
from multiprocessing.dummy import Pool
import time
import json
from chp_data.bkb_handler import BkbDataHandler

class ImmutableHint:
    def __init__(self, dictionary):
        self.hint = dictionary
    
    def getHint(self):
        return self.hint
    
    def __hash__(self):
        return hash(json.dumps(sorted(self.hint.items())))

class PatientProcessor:
    #load in a patient file
    def readPatientFile(self, fileName):
        with open(fileName, 'rb') as f_:
            patient_data = pickle.load(f_)
        return patient_data

    def processGenes(self, patient_data):

        def getGenes():
            genes = set()
            for patient_id, data in patient_data.items():
                genes.update(data['Patient_Genes'])
            return genes

        def getAllGeneHints(genes):
            def getGeneHints(gene):
                hint = Hint()
                hints = hint.query(gene)
                hints = hints['Gene']
                for ht in hints:
                    ht['chp_name'] = gene
                return hints

            start_time = time.time()

            gene_hints = []

            threadSize = 10

            for i in range(0,len(genes),threadSize):
                start_time = time.time()
                pool = Pool(threadSize)
                async_resp = pool.map(getGeneHints, list(genes)[:i + threadSize])
                pool.close()
                pool.join()
                print('Pool time = {}'.format(time.time() - start_time))

                for resp in async_resp:
                    gene_hints + resp

                pool.close()
                pool.join()
            return gene_hints

        genes = getGenes()
        print(type(genes))
        geneHints = getAllGeneHints(genes[:10])
        return geneHints

    def processDrugs(self, patient_data):
        def getDrugs():
            drugs = set()
            for patient_id, data in patient_data.items():
                drugs.update(data['Drug_Name(s)'])
            return drugs

        def getAllDrugHints(drugs):
            def getDrugHints(drug):
                hint = Hint()
                ht = hint.query(drug)
                return {ImmutableHint(ht): drug}
                
            start_time = time.time()

            drug_hints = {}

            threadSize = 10

            for i in range(0,len(drugs),threadSize):
                start_time = time.time()
                pool = Pool(threadSize)
                async_resp = pool.map(getDrugHints, list(drugs)[:i + threadSize])
                pool.close()
                pool.join()
                print('Pool time = {}'.format(time.time() - start_time))
                
                for resp in async_resp:
                    drug_hints.update(resp)

                pool.close()
                pool.join()
            return drug_hints

        drugs = getDrugs()

        drugHints = getAllDrugHints(drugs)
        return drugHints

    def processCancers(self, patient_data):
        def getCancerTypes():
            cancerTypes = set()
            for patient_id, data in patient_data.items():
                cancerTypes.update(data['Cancer_Type'])
            return cancerTypes

        def getAllCancerHints(cancers):
            def getCancerHints(cancer):
                hint = Hint()
                ht = hint.query(cancer)
                return {ImmutableHint(ht): cancer}
                
            start_time = time.time()

            cancer_hints = {}

            threadSize = 10

            for i in range(0,len(cancers),threadSize):
                start_time = time.time()
                pool = Pool(threadSize)
                async_resp = pool.map(getCancerHints, list(cancers)[:i + threadSize])
                pool.close()
                pool.join()
                print('Pool time = {}'.format(time.time() - start_time))
                
                for resp in async_resp:
                    cancer_hints.update(resp)

                pool.close()
                pool.join()
            return cancer_hints

        cancers = getCancerTypes()

        drugHints = getAllCancerHints(cancers)
        return drugHints

    def writeDictionaryToFile(self, dictionary, filename):
        with open(filename, 'wb') as handle:
            pickle.dump(dictionary, handle, protocol=pickle.HIGHEST_PROTOCOL)

class CurrieLookUp:
    def readCurrieFile(self, filename):
        self.dictionary = pickle.load(open(filename, 'rb'))

    def lookUp(self, currie):
        currie = ImmutableHint(currie)
        gene = self.dictionary[currie]

if __name__ == '__main__':
    bkb_handler = BkbDataHandler(dataset_version='1.1')
    patient_processor = PatientProcessor()
    patient_data = patient_processor.readPatientFile(bkb_handler.patient_data_pk_path)

    gene_hints = patient_processor.processGenes(patient_data)
    print(gene_hints)
