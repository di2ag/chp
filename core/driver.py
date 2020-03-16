import os
import sys
import argparse
import csv
import time
import pickle

from pybkb.core.common.bayesianKnowledgeBase import bayesianKnowledgeBase as BKB

sys.path.append('/home/cyakaboski/src/python/projects/bkb-pathway-provider/core')
#sys.path.append('/home/ncats/live/core')

from reasoner import Reasoner
from query import Query

#-- Reasoner Controls
INTERPOLATION = 'independence'
TARGET_STRATEGY = 'explicit'
MAX_NEW_EV = 5

AVALIABLE_DEMOGRAPHICS= set(['Gender', 'Process_Activity(s)', 'PathM', 'Process_Type(s)',\
                             'PathN', 'Survival_Time', 'Biological_Object(s)', 'PathT', 'Age_of_Diagnosis'])

class Driver:
    def __init__(self, config_file):
        #-- Read in configuration from config file.
        self.config = dict()
        with open(config_file, 'r') as csv_file:
            reader = csv.reader(csv_file)
            for row in reader:
                self.config[row[0]] = row[1]
        self.fused_bkb_path = self.config['fused_bkb_path']
        self.ncats_dir = self.config['ncats_dir']
        self.src_metadata_path = self.config['src_metadata_path']
        self.gene_var_direct = self.config['gene_var_direct']

        self.fused_bkb = BKB()
        self.fused_bkb.load(self.fused_bkb_path)
        self.reasoner = Reasoner(self.fused_bkb, gene_var_direct=self.gene_var_direct, max_new_ev=MAX_NEW_EV)
        self.reasoner.set_src_metadata(self.src_metadata_path)

    def run_query(self, query):
        result_query = self.reasoner.analyze_query(query, target_strategy=TARGET_STRATEGY, interpolation=INTERPOLATION)
        return result_query

    def run(self, i=0):
        if i == 0:
            print('Welcome to the BKB Pathways Driver.')
            input('Press any key to begin reasoning...')

        print("Let's build a query!")
        meta_evidence = self.chooseDemoEvidence()
        evidence = self.chooseStandardEvidence()
        meta_targets = self.chooseDemoTargets()
        targets = self.chooseStandardTargets()

        while True:
            print('Choose reasoning Type:\n1)\tRevision\n2)\tUpdating')
            reason_type = int(input('Your Choice: '))
            if reason_type == 1:
                reason_type = 'revision'
                break
            elif reason_type == 2:
                reason_type = 'updating'
                break
            else:
                print('Unrecognized reasoning type.')

        query = Query(evidence=evidence,
                      targets=targets,
                      meta_evidence=meta_evidence,
                      meta_targets=meta_targets,
                      type=reason_type)
        #print(meta_evidence)
        query = self.reasoner.analyze_query(query)
        query.getReport()
        again = input('Do you wish to reason again? ([y],n): ') or 'y'
        if again == 'y':
            self.run(i=i+1)
        else:
            return

    def chooseDemoEvidence(self):
        meta_evidence = list()
        while True:
            print('Choose Demographic Evidence:')
            for j, label in enumerate(self.reasoner.metadata_labels):
                print('{})\t{}'.format(j, label))
            print('{})\tNo More Demographic Evidence'.format(j+1))
            rv = int(input('Your Choice: '))
            if rv == j+1:
                return meta_evidence
            else:
                rvName = self.reasoner.metadata_labels[rv]
                while True:
                    op = input('Choose a operation ([==], >=, <=): ') or '=='
                    if op == '==' or op == '>=' or op == '<=':
                        print('Choose a value for demographic evidence.')
                        range_ = self.reasoner.metadata_ranges[rvName]
                        if type(range_) == set:
                            for item  in range_:
                                print('\t{}'.format(item))
                        else:
                            print('\tmin={}\n\tmax={}'.format(range_[0], range_[1]))
                        state = input('Your choice: ')
                        try:
                            state = float(state)
                        except:
                            pass
                        break
                    else:
                        print('Unrecognized operation!')
                while True:
                    print('This is your demographic choice:')
                    print('\t{} {} {}'.format(rvName, op, state))
                    correct = input('Is it correct? ([y], n): ') or 'y'
                    if correct == 'y':
                        meta_evidence.append((rvName, op, state))
                        break
                    elif correct == 'n':
                        break
                    else:
                        print('Unrecognized selection.')

    def chooseDemoTargets(self):
        meta_targets = list()
        while True:
            print('Choose Demographic Targets:')
            for j, label in enumerate(self.reasoner.metadata_labels):
                print('{})\t{}'.format(j, label))
            print('{})\tNo More Demographic Targets'.format(j+1))
            rv = int(input('Your Choice: '))
            if rv == j+1:
                return meta_targets
            else:
                rvName = self.reasoner.metadata_labels[rv]
                while True:
                    op = input('Choose a operation ([==], >=, <=): ') or '=='
                    if op == '==' or op == '>=' or op == '<=':
                        print('Choose a value for demographic target.')
                        range_ = self.reasoner.metadata_ranges[rvName]
                        if type(range_) == set:
                            for item  in range_:
                                print('\t{}'.format(item))
                        else:
                            print('\tmin={}\n\tmax={}'.format(range_[0], range_[1]))
                        state = input('Your choice: ')
                        try:
                            state = float(state)
                        except:
                            pass
                        break
                    else:
                        print('Unrecognized operation!')
                while True:
                    print('This is your demographic choice:')
                    print('\t{} {} {}'.format(rvName, op, state))
                    correct = input('Is it correct? ([y], n): ') or 'y'
                    if correct == 'y':
                        meta_targets.append((rvName, op, state))
                        break
                    elif correct == 'n':
                        break
                    else:
                        print('Unrecognized selection.')

    def chooseStandardEvidence(self):
        evidence = dict()
        while True:
            print('Choose Evidence:')
            time.sleep(1)
            for comp_idx in self.fused_bkb.getAllComponentIndices():
                comp_name = self.fused_bkb.getComponentName(comp_idx)
                if not 'Source' in comp_name and not 'source' in comp_name:
                    print('{})\t{}'.format(comp_idx, comp_name))
            print('{})\t{}'.format(comp_idx+1, 'Done selecting evidence.'))
            compIdx = int(input('Your Choice: '))
            if compIdx == comp_idx+1:
                return evidence
            print('Choose your state:')
            for state_idx in self.fused_bkb.getAllComponentINodeIndices(compIdx):
                state_name = self.fused_bkb.getComponentINodeName(compIdx, state_idx)
                print('{})\t{}'.format(state_idx, state_name))
            stateIdx = int(input('Your Choice: '))
            compName = self.fused_bkb.getComponentName(compIdx)
            stateName = self.fused_bkb.getComponentINodeName(compIdx, stateIdx)
            while True:
                print('This is your evidence choice:\n\t{} = {}'.format(compName, stateName))
                correct = input('Is that correct? ([y],n): ') or 'y'
                if correct == 'y':
                    evidence.update({compName: stateName})
                    break
                elif correct == 'n':
                    break
                else:
                    print('Unrecognized selection.')


    def chooseStandardTargets(self):
        targets = list()
        while True:
            print('Choose Targets:')
            time.sleep(1)
            for comp_idx in self.fused_bkb.getAllComponentIndices():
                comp_name = self.fused_bkb.getComponentName(comp_idx)
                if not 'Source' in comp_name and not 'source' in comp_name:
                    print('{})\t{}'.format(comp_idx, comp_name))
            print('{})\t{}'.format(comp_idx+1, 'Done selecting targets.'))
            compIdx = int(input('Your Choice: '))
            if compIdx == comp_idx+1:
                return targets
            compName = self.fused_bkb.getComponentName(compIdx)
            while True:
                print('This is your target choice:\n\t{}'.format(compName))
                correct = input('Is that correct? ([y],n): ') or 'y'
                if correct == 'y':
                    targets.append(compName)
                    break
                elif correct == 'n':
                    break
                else:
                    print('Unrecognized selection.')

    def collectVariables(self):
        #-- Get all inodes
        all_inode_names = self.fused_bkb.getINodeNames()
        #-- Filter out sources
        inode_names = [inode for inode in all_inode_names if 'Source' not in inode[0]]
        #-- Collect avaliable demographics
        demographics = {demo_name: list(demo_range)
                        for demo_name, demo_range in self.reasoner.metadata_ranges.items()
                        if demo_name in AVALIABLE_DEMOGRAPHICS}
        return {'genetic_info': inode_names,
                'demographic_info': demographics}

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='BKB Pathway Provider Minimal Driver')
    parser.add_argument('--config_file', default='driver.config', type=str)
    parser.add_argument('--headless', action='store_true')
    parser.add_argument('--query_file', type=str)
    parser.add_argument('--save_dir', type=str, default=os.getcwd())
    parser.add_argument('--get_variables', type=str, default=None)
    args = parser.parse_args()

    if args.headless:
        driver = Driver(args.config_file)
        #-- Load Query from query file
        query = Query().read(args.query_file)
        result_query = driver.run_query(query)

        #-- Save Query
        result_query.save(args.save_dir, only_json=True)
    elif args.get_variables is not None:
        driver = Driver(args.config_file)
        vars_ = driver.collectVariables()
        with open(args.get_variables, 'wb') as f_:
            pickle.dump(file=f_, obj=vars_)
    else:
        driver = Driver(args.config_file)
        driver.run()

