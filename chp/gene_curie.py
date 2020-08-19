import csv
import tqdm
import requests
import threading
import queue
import json

##########################################################
# curie_match
# Input: gene symbol
# Output: string Ensemble ID or 'not_found' if not found
#--------------------------------------------------------
# Description: cascade of api calls to determine esemble
# id. Order is determined by confidence in ensemble returned.
# If an ID is pulled earlier in the cascade, the ensemble is
# returned.

def curie_match(patient_gene):
    gene = patient_gene

    # symbol pass
    symbolPass, passReturn = symbol_pass(gene)
    if symbolPass:
        return passReturn

    # prev_symbol pass
    if not symbolPass:
        prevSymbolPass, passReturn = prev_symbol_pass(gene)
        if prevSymbolPass:
            return passReturn

    # allias_symbol pass
    if not prevSymbolPass:
        aliasSymbolPass, passReturn = alias_symbol_pass(gene)
        if aliasSymbolPass:
            return passReturn

    # entrez_id pass
    if not aliasSymbolPass:
        entrezIDPass, passReturn = entrez_id_pass(gene)
        if entrezIDPass:
            return passReturn

    return [gene, 'not_found', 'NA']

##########################################################
# symbol_pass
# Input: gene symbol
# Output: string Ensemble ID or 'not_found' if not found
#--------------------------------------------------------
# Description: The first in the api request cascade. Uses
# genenames using the symbol api request. Uses new
# symbol to request genenames using symbol.

def symbol_pass(gene):
    url = 'http://rest.genenames.org/fetch/symbol/{}'.format(gene)
    header = {'Accept': 'application/json'}
    r = requests.get(url, headers=header)
    r_json = r.json()

    response = r_json['response']
    hits = response['numFound']
    # default symbol pass
    symbolPass = False
    if int(hits) > 0:
        docs = response['docs']
        if 'ensembl_gene_id' in list(docs[0].keys()):
            curie = docs[0]['ensembl_gene_id']
            return True, [gene, curie, 'symbol']
    return False, [gene, 'not found', 'symbol']

##########################################################
# prev_symbol_pass
# Input: gene symbol
# Output: string Ensemble ID or 'not_found' if not found
#--------------------------------------------------------
# Description: The second in the api request cascade. Uses
# genenames using the prev_symbol api request. Uses new
# symbol to request genenames using symbol.

def prev_symbol_pass(gene):
    url = 'http://rest.genenames.org/fetch/prev_symbol/{}'.format(gene)
    header = {'Accept': 'application/json'}
    r = requests.get(url, headers=header)
    r_json = r.json()

    response = r_json['response']
    hits = response['numFound']
    if int(hits) > 0:
        docs = response['docs']
        if 'symbol' in list(docs[0].keys()):
            new_gene = docs[0]['symbol']
            url = 'http://rest.genenames.org/fetch/symbol/{}'.format(new_gene)
            r = requests.get(url, headers=header)
            r_json = r.json()

            response = r_json['response']
            hits = response['numFound']
            if int(hits) > 0:
                docs = response['docs']
                if 'ensembl_gene_id' in list(docs[0].keys()):
                    curie = docs[0]['ensembl_gene_id']
                    return True, [gene, curie, 'prev_symbol']
    return False, [gene, 'not found', 'prev_symbol']

##########################################################
# alias_symbol_pass
# Input: gene symbol
# Output: string Ensemble ID or 'not_found' if not found
#--------------------------------------------------------
# Description: The third in the api request cascade. Uses
# genenames using the alias_symbol api request. Uses new
# symbol to request genenames using symbol.

def alias_symbol_pass(gene):
    url = 'http://rest.genenames.org/fetch/alias_symbol/{}'.format(gene)
    header = {'Accept': 'application/json'}
    r = requests.get(url, headers=header)
    r_json = r.json()

    response = r_json['response']
    hits = response['numFound']
    if int(hits) > 0:
        docs = response['docs']
        if 'symbol' in list(docs[0].keys()):
            new_gene = docs[0]['symbol']
            url = 'http://rest.genenames.org/fetch/symbol/{}'.format(new_gene)
            r = requests.get(url, headers=header)
            r_json = r.json()

            response = r_json['response']
            hits = response['numFound']
            if int(hits) > 0:
                docs = response['docs']
                if 'ensembl_gene_id' in list(docs[0].keys()):
                    curie = docs[0]['ensembl_gene_id']
                    return True, [gene, curie, 'alias_symbol']
    return False, [gene, 'not found', 'alias_symbol']

##########################################################
# entrez_id_pass
# Input: gene symbol
# Output: string Ensemble ID or 'not_found' if not found
#--------------------------------------------------------
# Description: The fourth in the api request cascade. Uses
# mygene.info to extract an extrez_id. request is specified
# for human species only and top return. Uses entrez to 
# request genenames using the entrez_id request

def entrez_id_pass(gene):
    url = "https://mygene.info/v3/query?q=RAF1&species=human&size=1"
    header = {'Accept': 'application/json'}
    r = requests.get(url, headers=header)
    r_json = r.json()

    hits = r_json['hits']
    if len(hits) == 1:
        docs = hits
        if '_entrezgene' in list(docs[0].keys()):
            entrez_id = docs[0]['_entrezgene']
            url = 'http://rest.genenames.org/fetch/entrez_id/{}'.format(entrez_id)
            r = requests.get(url, headers=header)
            r_json = r.json()

            response = r_json['response']
            hits = response['numFound']
            if int(hits) > 0:
                docs = response['docs']
                if 'ensembl_gene_id' in list(docs[0].keys()):
                    curie = docs[0]['ensembl_gene_id']
                    return True, [gene, curie, 'entrez_id']
    return False, [gene, 'not found', 'etrez_id']




##########################################################
# get_gene_list
# input: patient file location
# output: list of unique gene symbols
#--------------------------------------------------------
# description: searches through each patient and appends
# unique gene symbol names to list of genes

def get_gene_list(patient_dir):
    # reading initial patient data
    with open(patient_dir, 'r') as csv_file:
        reader = csv.reader(csv_file)
        # row[6] = mutated gene
        next(reader)
        rows = [row[6] for row in reader]

    patient_gene_list = list()
    # gather all relevant records
    for row in rows: #tqdm.tqdm(rows, desc='Reading patient files'):
        gene = str(row)
        if gene not in patient_gene_list:
            patient_gene_list.append(gene)
    csv_file.close()
    return patient_gene_list


##########################################################
# spinoff_requests
# input: list of genes, max threads
# output: annotated curies for each gene
#--------------------------------------------------------
# description: iterates gene list and spawns of threads
# for each gene. Working threads are held in a queue
# and can only contain up to max threads at a time before
# calling a join on each thread. 

def spinoff_requests(gene_list, max_threads):
    # thread safe queue
    que = queue.Queue()
    threads = list()
    results = list()
    counter = 0
    for pat_gene in tqdm.tqdm(gene_list, desc='sending gene curie requests'):
        t = threading.Thread(target = lambda q, arg1: q.put(curie_match(arg1)), args=(que, pat_gene))
        threads.append(t)
        if len(threads) == max_threads or gene_list[-1] == pat_gene:
            for thread in threads:
                thread.start()
            for thread in threads:
                thread.join()
                results.append(que.get())
            threads = list()
            que.empty()
    return results

if __name__ == "__main__":

    patientMutGenes = '/home/public/data/ncats/data_drop_03-04-2020/wxs.csv'
    max_threads = 64

    gene_list = get_gene_list(patientMutGenes)
    results = spinoff_requests(gene_list, max_threads)

    # save out request results
    f = open('gene_curie_map.csv', 'w')
    f.write('Gene,Curie,find_type\n')
    for r in results:
        row = str(r[0])+','+str(r[1])+','+str(r[2])+'\n'
        f.write(row)
    f.close()
