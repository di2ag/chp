import os
import copy
import argparse
import csv
import random
import numpy as np

REACTONS_DNA_MISMATCH = [
            'MLH1 variants-defective DNA mismatch repair',
            'MSH2 variant:MSH6-defective DNA mismatch repair',
            'MSH2 variant:MSH3-defective DNA mismatch repair',
            'MSH3 variant:MSH2-defective DNA mismatch repair',
            'MSH6 variant:MSH2-defective DNA mismatch repair',
            'PMS2 variants-defective DNA mismatch repair'
]

REACTIONS_TO_PATHWAYS = {
    'MLH1 variants-defective DNA mismatch repair': 'Defective Mismatch Repair Associated With MLH1',
    'MSH2 variant:MSH6-defective DNA mismatch repair': 'Defective Mismatch Repair Associated With MSH2',
    'MSH2 variant:MSH3-defective DNA mismatch repair': 'MSH3 variant:MSH2-defective DNA mismatch repair',
    'MSH3 variant:MSH2-defective DNA mismatch repair': 'Defective Mismatch Repair Associated With MSH3',
    'MSH6 variant:MSH2-defective DNA mismatch repair': 'Defective Mismatch Repair Associated With MSH6',
    'PMS2 variants-defective DNA mismatch repair': 'Defective Mismatch Repair Associated With PMS2',
}

PATHWAYS_TO_PATHWAYS = {
    'Defective Mismatch Repair Associated With MLH1': 'Diseases of Mismatch Repair (MMR)',
    'Defective Mismatch Repair Associated With MSH2':'Diseases of Mismatch Repair (MMR)',
    'MSH3 variant:MSH2-defective DNA mismatch repair':'Diseases of Mismatch Repair (MMR)',
    'Defective Mismatch Repair Associated With MSH3':'Diseases of Mismatch Repair (MMR)',
    'Defective Mismatch Repair Associated With MSH6':'Diseases of Mismatch Repair (MMR)',
    'Defective Mismatch Repair Associated With PMS2':'Diseases of Mismatch Repair (MMR)',
}

def simulate_gene_to_reaction_snodes(filename, tcga_expression_file, num_genes, scenario='mismatch'):
    #-- Pick some genes randomly that are involved in each reaction in scenario.
    if scenario == 'mismatch':
        reaction_to_genes = dict()
        genes = ['Gene_{}'.format(i) for i in range(num_genes)]
        for reaction in REACTONS_DNA_MISMATCH:
            rand_1 = random.randint(1, num_genes)
            reaction_to_genes[reaction] = random.sample(genes, rand_1)
    else:
        raise ValueError('Unknown scenario')

    #-- Read in the expression file to get counts
    with open(tcga_expression_file, 'r') as csv_file:
        reader = csv.reader(csv_file)
        gene_counts = {gene: [] for gene in genes}
        for row in reader:
            gene_ = row[-2]
            count_ = row[-1]
            gene_counts[gene_].append(count_)

    #-- Now make CSV s-node file
    header = ['reaction_name']
    for gene in genes:
        header.extend([gene+'_low', gene+'_high'])
    lines = [header]
    for reaction, react_genes in reaction_to_genes.items():
        line = [reaction]
        for gene_ in genes:
            if gene_ in react_genes:
                counts = np.asarray(gene_counts[gene_]).astype(np.int64)
                #print(counts)
                line.extend([min(counts), max(counts)])
            else:
                line.extend(['',''])
        lines.append(line)
    react_genes_snodes = filename + '-react_gene.snodes'
    react_pathway_snodes = filename + '-pathways.snodes'
    with open(react_genes_snodes, 'w') as f_:
        writer = csv.writer(f_)
        writer.writerows(lines)
    print('Wrote file to: {}'.format(react_genes_snodes))

    with open(react_pathway_snodes, 'w') as csv_file:
        writer = csv.writer(csv_file)
        writer.writerows(list(REACTIONS_TO_PATHWAYS.items()))
        writer.writerows(list(PATHWAYS_TO_PATHWAYS.items()))
    print('Wrote file to: {}'.format(react_pathway_snodes))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('file')
    parser.add_argument('tcga_file')
    parser.add_argument('num_genes', type=int)
    parser.add_argument('--scenario', default='mismatch')

    args = parser.parse_args()

    simulate_gene_to_reaction_snodes(args.file, args.tcga_file, args.num_genes, args.scenario)

