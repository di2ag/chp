import csv
import tqdm
import requests
import threading
import queue
from multiprocessing import Pool
import json

def curie_match(patient_gene):
    gene = patient_gene
    url = 'https://www.genecards.org/cgi-bin/carddisp.pl?gene={}'.format(gene)
    header = {'User-agent': 'Mozilla/5.0'}
    r = requests.get(url, headers=header)
    r_text = r.text.encode('ascii','ignore').decode('utf-8')

    alias_start = r_text.find('Aliases for')
    external_ids = r_text.find('External Ids')

    curie = ''
    if alias_start != -1 and external_ids != -1:
        alias_interval = r_text[alias_start:external_ids]
        # direct string match
        if 'Aliases for {} Gene'.format(gene) in alias_interval:
            external_ids_end = r_text.find('</div>', external_ids)
            id_interval = r_text[external_ids:external_ids_end]
            ids = id_interval.split('<li>')
            for id in ids:
                if 'Ensembl' in id:
                    curie = id.split('=')[3].split('\"')[0]
        else:
            aliases = alias_interval.split('<li>')
            geneFound = False
            for i in range(1, len(aliases)):
                alias = aliases[i].split('<sup>')[0].strip()
                if gene == alias:
                    geneFound = True
            # alias match
            if geneFound:
                external_ids_end = r_text.find('</div>', external_ids)
                id_interval = r_text[external_ids:external_ids_end]
                ids = id_interval.split('<li>')
                for id in ids:
                    if 'Ensembl' in id:
                        curie = id.split('=')[3].split('\"')[0]
    if curie == '':
        return [patient_gene,'HELP']
    return [patient_gene,curie]

patientMutGenes = '/home/public/data/ncats/data_drop_03-04-2020/wxs.csv'
patientCurieFile = '/home/public/data/ncats/data_drop_02-11-2020/rnaseq_fpkm_uq_primary_tumor.csv'

# reading initial patient data
with open(patientMutGenes, 'r') as csv_file:
    reader = csv.reader(csv_file)
    patientGeneDict = dict()
    # row[0] = cancer type row[1] = patientID, row[6] = mutated gene
    next(reader)
    rows = [(row[0],row[1],row[6]) for row in reader]

print(len(rows))
patient_gene_list = list()

# gather all relevant records
for row in tqdm.tqdm(rows, desc='Reading patient files'):
    gene = str(row[2])
    if gene not in patient_gene_list:
        patient_gene_list.append(gene)
csv_file.close()

print(len(patient_gene_list))

#pool = Pool(processes=80)
#data_outputs = pool.map(curie_match, patient_gene_list)

# thread safe queue
que = queue.Queue()
threads = list()
max_threads = 128
results = list()
for pat_gene in tqdm.tqdm(patient_gene_list, desc='sending gene curie requests'):
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

f = open('curie_map.csv', 'w')
f.write('Gene,Curie\n')

for r in results:
    row = str(r[0])+','+str(r[1])+'\n'
    f.write(row)
f.close()
