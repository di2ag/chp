import csv
import tqdm
import requests
import threading
import queue
import json

##########################################################
# curie_match
# Input: gene symbol
# Output: string Ensemble ID or 'NOT FOUND' if not found
#--------------------------------------------------------
# Description: Parsing is done via genecards.org. This
# does not access API code, but rather scrapes raw html
# and parses. Uses three indices to parse:
#    1.) alias_start - Aliases matching. While searches
#        are deterministic, it appears there is a 1 to N
#        relationship with aliases. Must check these
#        aliases to ensure there is a match
#    2.) enternal_id - location where curies begin
#    3.) external_ids_end - location of external id end

def curie_match(patient_gene):
    gene = patient_gene
    url = 'https://www.genecards.org/cgi-bin/carddisp.pl?gene={}'.format(gene)
    header = {'User-agent': 'Mozilla/5.0'}
    r = requests.get(url, headers=header)
    r_text = r.text.encode('ascii','ignore').decode('utf-8')

    # alias start
    alias_start = r_text.find('Aliases for')
    # external ids start
    external_ids = r_text.find('External Ids')

    # checks if we have 1 - N relationship on ensemble IDs
    ensembl_hits = 0

    curie = ''
    if alias_start != -1 and external_ids != -1:
        alias_interval = r_text[alias_start:external_ids]
        # direct match with genecard gene symbol as opposed to an alias
        if 'Aliases for {} Gene'.format(gene) in alias_interval:
            # external ids end
            external_ids_end = r_text.find('</div>', external_ids)
            id_interval = r_text[external_ids:external_ids_end]
            ids = id_interval.split('<li>')
            for id in ids:
                if 'Ensembl' in id:
                    curie = id.split('=')[3].split('\"')[0]
                    ensembl_hits += 1
        else:
            aliases = alias_interval.split('<li>')
            geneFound = False
            for i in range(1, len(aliases)):
                alias = aliases[i].split('<sup>')[0].strip()
                if gene == alias:
                    geneFound = True
            # direct match with genecard alias
            if geneFound:
                # external ids end
                external_ids_end = r_text.find('</div>', external_ids)
                id_interval = r_text[external_ids:external_ids_end]
                ids = id_interval.split('<li>')
                for id in ids:
                    if 'Ensembl' in id:
                        curie = id.split('=')[3].split('\"')[0]
                        ensembl_hits += 1
    if curie == '':
        return [patient_gene,'NOT FOUND', ensembl_hits]
    return [patient_gene,curie,ensembl_hits]


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
        patientGeneDict = dict()
        # row[6] = mutated gene
        next(reader)
        rows = [row[6] for row in reader]

    patient_gene_list = list()
    # gather all relevant records
    for row in tqdm.tqdm(rows, desc='Reading patient files'):
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
    for pat_gene in tqdm.tqdm(gene_list, desc='sending gene curie requests'):
        if len(threads) == max_threads:
            for thread in threads:
                thread.start()
            for thread in threads:
                thread.join()
                results.append(que.get())
            threads = list()
        else:
            t = threading.Thread(target = lambda q, arg1: q.put(curie_match(arg1)), args=(que, pat_gene))
            threads.append(t)
    return results

if __name__ == "__main__":

    patientMutGenes = '/home/public/data/ncats/data_drop_03-04-2020/wxs.csv'
    max_threads = 128

    gene_list = get_gene_list(patientMutGenes)
    results = spinoff_requests(gene_list, max_threads)

    # save out request results
    f = open('gene_curie_map.csv', 'w')
    f.write('Gene,Curie,Ensembl_hits\n')
    for r in results:
        row = str(r[0])+','+str(r[1])+','+str(r[2])+'\n'
        f.write(row)
    f.close()
