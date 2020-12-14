import os
import copy
import argparse
import csv
import random

MAX_EXPRESSION = 5000

def simulate_expression_data(filepath, num_genes, num_patients, num_projects):
    lines = []
    for i in range(num_projects):
        for j in range(num_patients):
            for k in range(num_genes):
                line = ['patient_{}'.format(j),'project_{}'.format(i),'Gene_{}'.format(k),random.randint(0, MAX_EXPRESSION)]
                lines.append(line)
    with open(filepath, 'w') as f_:
        writer = csv.writer(f_)
        writer.writerows(lines)
    print('Wrote file to: {}'.format(filepath))

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('file')
    parser.add_argument('num_genes', type=int)
    parser.add_argument('num_patients', type=int)
    parser.add_argument('num_projects', type=int)

    args = parser.parse_args()

    simulate_expression_data(args.file, args.num_genes, args.num_patients, args.num_projects)

