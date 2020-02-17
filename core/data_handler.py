import pandas as pd
import os
import pickle
import tqdm

class TcgaDataHandler:
    def __init__(self,
                 tumor_rna_expression_file=None,
                 normal_rna_expression_file=None,
                 clinical_data_file=None,
                 drug_data_file=None,
                 radiation_data_file=None,
                 mutation_data_file=None,
                 hash_files=False):

        self.tumor_rna_expression_file = tumor_rna_expression_file
        self.normal_rna_expression_file = normal_rna_expression_file
        self.clinical_data_file = clinical_data_file
        self.drug_data_file = drug_data_file
        self.radiation_data_file = radiation_data_file
        self.mutation_data_file = mutation_data_file
        self.hash_files = hash_files

    def __load_csv_file(self, data_file, type='rna_expression'):
        print('Loading Data...')
        if type == 'rna_expression':
            df = pd.read_csv(data_file)
            print('Loading Complete.')
        elif type == 'clinical':
            df = pd.read_csv(data_file)
            print('Loading Complete.')
        elif type == 'drug':
            df = pd.read_csv(data_file)
            print('Loading Complete.')
        elif type == 'radiation':
            df = pd.read_csv(data_file)
        elif type == 'mutation':
            df = pd.read_csv(data_file)
            print('Loading Complete.')
        else:
            raise ValueError('Data File type: {} not recognized'.format(type))

        return df

    def hashDataFiles(self, write_files=False, path='', file_prefix=''):
        #-- Hash Tumor RNA Expression
        if self.tumor_rna_expression_file is not None:
            #-- Helper Hasher Function
            def tumor_rna_hasher(master_dict, row):
                if row['Cancer Type'] not in master_dict:
                    master_dict[row['Cancer Type']] = dict()
                if row['Patient ID'] not in master_dict[row['Cancer Type']]:
                    master_dict[row['Cancer Type']][row['Patient ID']] = {'Gene': dict()}
                master_dict[row['Cancer Type']][row['Patient ID']]['Sample ID'] = row['Sample ID']
                master_dict[row['Cancer Type']][row['Patient ID']]['Portion ID'] = row['Portion ID']
                master_dict[row['Cancer Type']][row['Patient ID']]['Analyte ID'] = row['Analyte ID']
                master_dict[row['Cancer Type']][row['Patient ID']]['Aliquot ID'] = row['Aliquot ID']
                if row['Gene'] not in  master_dict[row['Cancer Type']][row['Patient ID']]['Gene']:
                    master_dict[row['Cancer Type']][row['Patient ID']]['Gene'][row['Gene']] = row['HTSeq FPKM-UQ']
                return master_dict

            tumor_rna_expression_df = self.__load_csv_file(self.tumor_rna_expression_file)
            self.tumor_rna_expression_dict = dict()
            print('Hashing...')
            for idx in tqdm.tqdm(tumor_rna_expression_df.index):
                row = tumor_rna_expression_df.loc[idx]
                self.tumor_rna_expression_dict = tumor_rna_hasher(self.tumor_rna_expression_dict, row)

            #self.tumor_rna_expression_dict = tumor_rna_expression_df.to_dict()
            if write_files:
                print('Writing Hashfile...')
                if path == '':
                    path = os.getcwd()
                with open(os.path.join(path, file_prefix + 'tumor_rna_expression.pk'), 'wb') as f_:
                    pickle.dump(file=f_, obj=self.tumor_rna_expression_dict)
                print('Finished writing.')
        return True

    def calculateGeneStatistics(self):
        if self.tumor_rna_expression_file is not None:
            rna_data = self.__load_csv_file(self.tumor_rna_expression_file)
        else:
            raise ValueError('No RNA file was specified.')

        #-- Calculate Gene Expression Level Categories over entire population through STD method.
        gene_express_data = rna_data[['Gene','HTSeq FPKM-UQ']].copy()
        pop_gene_data_means = gene_express_data.groupby('Gene').mean()
        pop_gene_data_means = pop_gene_data_means.rename(columns={'HTSeq FPKM-UQ': ('HTSeq FPKM-UQ', 'Mean')})
        pop_gene_data_stds = gene_express_data.groupby('Gene').std()
        pop_gene_data_stds = pop_gene_data_stds.rename(columns={'HTSeq FPKM-UQ': ('HTSeq FPKM-UQ','Std')})
        pop_gene_data = pd.concat([pop_gene_data_means, pop_gene_data_stds], axis=1)
        return pop_gene_data


class ReactomeDataHandler:
    def __init__(self,
                 r2g_file=None,
                 p2p_file=None):
        self.r2g_file = r2g_file
        self.p2p_file = p2p_file

    def buildPathwayBKB(self):
        #-- Process Pathway BKFs
        with open(pathways_file, 'r') as csv_file:
            reader = csv.reader(csv_file)
            bkf = BKB()
            for row in reader:
                tail = row[0]
                head = row[1]
                states = ['True', 'False']
                if bkf.findComponent(tail) == -1:
                    component_tail = BKB_component(tail)
                    for state in states:
                        component_tail.addINode(BKB_I_node(state, component_tail))
                    bkf.addComponent(component_tail)
                else:
                    component_tail = bkf.findComponent(tail)

                if bkf.findComponent(head) == -1:
                    component_head = BKB_component(head)
                    for state in states:
                        component_head.addINode(BKB_I_node(state, component_head))
                    bkf.addComponent(component_head)
                else:
                    component_head = bkf.findComponent(head)

                for state1 in states:
                    for state2 in states:
                        bkf.addSNode(BKB_S_node(component_head,
                                                component_head.findState(state1),
                                                random.random(),
                                                [(component_tail, component_tail.findState(state2))])
                                    )


if __name__ == '__main__':
    data_handler = TcgaDataHandler(
        tumor_rna_expression_file='/home/public/data/ncats/data_drop_02-11-2020/rnaseq_fpkm_uq_primary_tumor.csv'
    )
    data_handler.hashDataFiles(write_files=True)
    data_handler.calculateGeneStatistics()
