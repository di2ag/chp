"""
    Source code developed by DI2AG.
    Thayer School of Engineering at Dartmouth College
    Authors:    Dr. Eugene Santos, Jr
                Mr. Chase Yakaboski,
                Mr. Gregory Hyde,
                Mr. Luke Veenhuis,
                Dr. Keum Joo Kim
"""

import unittest
import copy
from chp.reasoner import Reasoner
from chp.query import Query
from chp_data.patient_bkb_builder import PatientBkbBuilder

class testSyntheticPatients(unittest.TestCase):

    """ Test suite for the synthetic data sets. The primary goal
        of these tests are to verify the integrity of the BKB structure.
        There are three primary tests:
            1.) test_synthetic_probs (Done)
            2.) test_synthetic_contributions (Done)
            3.) test_synthetic_sensitivities (TODO)
    """

    def test_synthetic_probs(self):
        """ Description: Builds a set of 11 patient sets all with a single
            gene mutation. In each sub set of patients the set is split
            in two. Half receive some drug treatment and half do not.
            Each preceting set of the 11 has a stronger drug treatment
            signal. That is, survival percentages (of patients >= 1000 days)
            for set 1 with drug treatment is 0% and without is 100%.
            For set two that is 10% and 90% respectively, until set 11 where
            they are flipped to be 100% and 0% respectively. We assert the
            following:
            for sets 1-5 - P(survival_time >= 1000 days | gene) >
                           P(survival_time >= 1000 days | gene & drug)
            for set 6 - P(survival_time >= 1000 days | gene) ==
                        P(survival_time >= 1000 days | gene & drug)
            for set 7-11 - P(survival_time >= 1000 days | gene) <
                           P(survival_time >= 1000 days | gene & drug)

            Note - P(survival_time >= 1000 days | gene) is always 50%
            Note - P(survival_time >= 1000 days | gene & drug) is
                   0%, 10%, ..., 100% for sets 1-11.
        """

        patient_dicts = self.build_synthetic_prob_dicts(num_pats=200)
        fused_bkbs = self.build_fused_axle_bkbs(patient_dicts)
        query_responses = self.run_synthetic_prob_queries(fused_bkbs, patient_dicts)
        for i in range(0, len(query_responses)):
            if i < 5:
                self.assertGreater(query_responses[i]['no_drug'], query_responses[i]['drug'])
            elif i == 5:
                # to 7 decimal places
                self.assertAlmostEqual(query_responses[i]['no_drug'], query_responses[i]['drug'], 7)
            else:
                self.assertGreater(query_responses[i]['drug'], query_responses[i]['no_drug'])

    def test_synthetic_contributions(self):

        """ Builds a set of 11 patient sets all with a single drug.
            In each sub set the set is split in two. The first half
            has a survival time of 999, and the second 1001. The
            halved sets are also further split down but some having
            some mutated gene. The gene split is designed so that
            the gene contribution beloning to the < 1000 survival
            group grows in its contribution value. That is, 0% of
            patients with survival_time < 1000 have the mutated gene
            in set 1 of 11 and 100% have the gene where survival_time
            > 1000. In set 2 this is 10% and 90%, respectively until
            set 11 where they are flipped to be 100% and 0%, respectively.
            We assert the following:
            for sets 1-5 - contribution(gene) in patients with survival_time
                           < 1000 is < contribution(gene) in patients with
                           survival_time > 1000.
            for set 6 - contribution(gene) in patients with survival_time
                           = 1000 is < contribution(gene) in patients with
                           survival_time > 1000.
            for set 7-11 - contribution(gene) in patients with survival_time
                           > 1000 is < contribution(gene) in patients with
                           survival_time > 1000.
        """

        patient_dicts = self.build_synthetic_contribution_dicts(num_pats=200)
        fused_bkbs = self.build_fused_wildcard_bkbs(patient_dicts)
        query_responses = self.run_synthetic_contribution_queries(fused_bkbs, patient_dicts)
        for i in range(0, len(query_responses)):
            if i < 5:
                self.assertGreater(query_responses[i][0]['RAF1'], query_responses[i][1]['RAF1'])
            elif i == 5:
                # to 7 decimal places
                self.assertAlmostEqual(query_responses[i][0]['RAF1'], query_responses[i][1]['RAF1'], 7)
            else:
                self.assertGreater(query_responses[i][1]['RAF1'], query_responses[i][0]['RAF1'])

    def run_synthetic_contribution_queries(self, fused_bkbs, patient_dicts):

        """ Runs the following query on a wildcard BKB:
                P(survival_time >= 1000 | drug)
            Calculates the wild card gene contributions using
            the patient analysis.

            :param fused_bkbs: A set of fused wildcard BKBs
            :type fused_bkbs: list
            :param patient_dicts: a set of patient_dicts
            :type fused_bkbs: list
            :return query_responses: a list of gene contributions in both
                the true/false survival instantiations for each query
            :rtype list:
        """

        query_responses = []
        for i in range(0,11):
            reasoner = Reasoner(fused_bkb=fused_bkbs[i], patient_data=patient_dicts[i])
            drug_evidence = {'demo_CYCLOPHOSPHAMIDE':'True'}
            query_drug = Query(evidence=drug_evidence,
                               targets=list(),
                               meta_evidence=None,
                               meta_targets=[('survival_time', '>=', 1000)],
                               type='updating')
            query_drug = reasoner.analyze_query(copy.deepcopy(query_drug),
                                                save_dir=None,
                                                target_strategy='explicit',
                                                interpolation='standard')
            report = query_drug.jsonExplanations(contributions_include_srcs=False,
                                        contributions_top_n_inodes=30,
                                        contributions_ignore_prefixes=['_'])

            t_norm, f_norm, t_unnorm, f_unnorm = self.get_synthetic_query_probs(query_drug)

            analysis = {'patient_analysis': report['Patient Analysis'],
                        'contribution_analysis': report['Contributions Analysis']}

            if 'survival_time >= 1000 = True' in analysis['contribution_analysis'].keys():
                true_contrib = analysis['contribution_analysis']['survival_time >= 1000 = True']['demo_CYCLOPHOSPHAMIDE = True']
                true_pats = len(analysis['patient_analysis']['All Involved Patients']['survival_time >= 1000 = True'].keys())
                true_ind_cont = float(true_contrib)/float(true_pats)/t_unnorm
            else:
                true_ind_cont = 0
            if 'survival_time >= 1000 = False' in analysis['contribution_analysis'].keys():
                false_contrib = analysis['contribution_analysis']['survival_time >= 1000 = False']['demo_CYCLOPHOSPHAMIDE = True']
                false_pats = len(analysis['patient_analysis']['All Involved Patients']['survival_time >= 1000 = False'].keys())
                false_ind_cont = float(false_contrib)/float(false_pats)/f_unnorm

            # gene contributions
            true_gene_contrib = {'RAF1':0}
            false_gene_contrib = {'RAF1':0}
            for pat_key, pat in patient_dicts[i].items():
                if pat['survival_time'] > 1000 and 'CYCLOPHOSPHAMIDE' in pat['drug_curies']:
                    if 'RAF1' in pat['gene_curies']:
                        true_gene_contrib['RAF1'] += true_ind_cont
                elif pat['survival_time'] < 1000 and 'CYCLOPHOSPHAMIDE' in pat['drug_curies']:
                    if 'RAF1' in pat['gene_curies']:
                        false_gene_contrib['RAF1'] += false_ind_cont
            query_responses.append((true_gene_contrib, false_gene_contrib))
        return query_responses


    def run_synthetic_prob_queries(self, fused_bkbs, patient_dicts):

        """ Each patient set two queries must be run:
                1.) P(survival_time >= 1000 | gene & drug)
                2.) P(survival_time >= 1000 | gene)
            The probabilities for each set are returned in a list of
            dictionaries indexed by 'drug' and 'no_drug'.

            :param fused_bkbs: a list of fused AXLE bkbs
            :type fused_bkbs: list
            :param patient_dicts: a list of patient dict
            :type patient_dicts: list
            :return query_responses: a list of truth probability
                assignments.
            :rtype query_responses: list
        """

        query_responses = []
        for i in range(0, 11):
            # with drug
            reasoner = Reasoner(fused_bkb=fused_bkbs[i], patient_data=patient_dicts[i])
            gene_evidence = {'mut_RAF1' : 'True'}
            query_drug = Query(evidence=gene_evidence,
                               targets=list(),
                               meta_evidence=[('drug_curies', '==','CYCLOPHOSPHAMIDE')],
                               meta_targets=[('survival_time', '>=', 1000)],
                               type='updating')
            query_drug1 = reasoner.analyze_query(copy.deepcopy(query_drug),
                                                save_dir=None,
                                                target_strategy='explicit',
                                                interpolation='standard')
            truth_assignment_drug, f_norm, t_unnorm, f_unnorm = self.get_synthetic_query_probs(query_drug1)

            # with no drug
            reasoner = Reasoner(fused_bkb=fused_bkbs[i], patient_data=patient_dicts[i])
            query_no_drug = Query(evidence=gene_evidence,
                                  targets=list(),
                                  meta_evidence=None,
                                  meta_targets=[('survival_time', '>=', 1000)],
                                  type='updating')
            query_no_drug = reasoner.analyze_query(query_no_drug,
                                                   save_dir=None,
                                                   target_strategy='explicit',
                                                   interpolation='standard')
            truth_assignment_no_drug, f_norm, t_unnorm, f_unnorm = self.get_synthetic_query_probs(query_no_drug)

            # answers
            query_responses.append({'drug':truth_assignment_drug,
                                    'no_drug':truth_assignment_no_drug})
        return query_responses

    def get_synthetic_query_probs(self, query):

        """ Indexes and a run query's update to extract
            end probabilities. These are then normalized. Both
            normalized and unnormalized True/False assignments
            are returned.

            :param query: A query that has been inferenced
            :type query: Query
            :return: truth_assignment, false_assignment, t_un, f_un
                with are the True/False assignments normalized and
                unnormalized.
            :rtype: int
        """

        target_info = []
        for update, prob in query.result.updates.items():
            comp_idx, state_idx = update
            comp_name = query.bkb.getComponentName(comp_idx)
            state_name = query.bkb.getComponentINodeName(comp_idx, state_idx)
            target_info.append([comp_name, state_name, prob])

        if target_info[0][1] == 'True':
            truth_assignment = target_info[0][2]
            false_assignment = target_info[1][2]
            t_un = target_info[0][2]
            f_un = target_info[1][2]
        else:
            truth_assignment = target_info[1][2]
            false_assignment = target_info[0][2]
            t_un = target_info[1][2]
            f_un = target_info[0][2]

        if truth_assignment != -1 and false_assignment != -1:
            prob_sum = truth_assignment + false_assignment
            truth_assignment /= prob_sum
            false_assignment /= prob_sum
            sensitivities = True
        elif truth_assignment == -1 and false_assignment != -1:
            truth_assignment = 0
            prob_sum = truth_assignment + false_assignment
            truth_assignment /= prob_sum
            false_assignment /= prob_sum
        elif truth_assignment != -1 and false_assignment == -1:
            false_assignment = 0
            prob_sum = truth_assignment + false_assignment
            truth_assignment /= prob_sum
            false_assignment /= prob_sum
        return truth_assignment, false_assignment, t_un, f_un

    def build_fused_wildcard_bkbs(self, patient_dicts):

        """ Builds wildcard BKBs from the set of
            patient_dicts.

            :param patient_dicts: A list of patient dicts
            :type patient_dicts: list
            :return BKBs: a list of fused wilcard BKBs.
            :rtype: list
        """

        BKBs = []
        for i in range(0, 11):
            PBB = PatientBkbBuilder(patient_dicts[i])
            PBB.processPatientWildcardBKF('drug_curies')
            fused_bkb = PBB.fuseBKB()
            BKBs.append(fused_bkb)
        return BKBs

    def build_fused_axle_bkbs(self, patient_dicts):

        """ Builds AXLE BKBs from the set of
            patient_dicts.

            :param patient_dicts: A list of patient dicts
            :type patient_dicts: list
            :return BKBs: a list of fused AXLE BKBs
            :rtype: list
        """

        BKBs = []
        for i in range(0, 11):
            PBB = PatientBkbBuilder(patient_dicts[i])
            PBB.processAxleBKF()
            fused_bkb = PBB.fuseBKB()
            BKBs.append(fused_bkb)
        return BKBs

    def build_synthetic_contribution_dicts(self, num_pats=2000):

        """ Builds a set of 11 patient dicts. Each subset is
            built of of num_pats number of patients. The first half
            has a survival_time < 1000 and the second half > 1000.
            These subsets are further divided where the group with
            survival < 1000 having 0% containing a gene mutation
            and 100% of patients with survival_time > 1000 having a
            mutated gene in set 1. In each subsequent set this proportion
            changes such that set 2 has 10% and 90%, respectively until
            set 11 where these are flipped to be 100% and 0%, respectively.

            :param num_pats: number of patients to build in each subset
            :type num_pats: int
            :return patient_sets: a list of the patient_dicts
            :rtype patient_sets: list
        """

        patient_sets = []
        for i in range(0, 11):
            patient_i_dict = dict()
            set_divide = i/10

            # how many in group 1 have gene
            group_1_divide = round((num_pats/2) * set_divide)
            # how many in group 2 have gene
            group_2_divide = round((num_pats/2) * (1-set_divide)) + (num_pats/2)

            for j in range(0, num_pats):
                if j < num_pats/2:
                    if j < group_1_divide:
                        patient_i_dict[j] = {'patient_id':j,
                                             'survival_time':999,
                                             'drug_curies':tuple(['CYCLOPHOSPHAMIDE']),
                                             'gene_curies':tuple(['RAF1'])}
                    else:
                        patient_i_dict[j] = {'patient_id':j,
                                             'survival_time':999,
                                             'drug_curies':tuple(['CYCLOPHOSPHAMIDE']),
                                             'gene_curies':tuple()}
                else:
                    if j < group_2_divide:
                        patient_i_dict[j] = {'patient_id':j,
                                             'survival_time':1001,
                                             'drug_curies':tuple(['CYCLOPHOSPHAMIDE']),
                                             'gene_curies':tuple(['RAF1'])}
                    else:
                        patient_i_dict[j] = {'patient_id':j,
                                             'survival_time':1001,
                                             'drug_curies':tuple(['CYCLOPHOSPHAMIDE']),
                                             'gene_curies':tuple()}
            patient_sets.append(patient_i_dict)
        return patient_sets

    def build_synthetic_prob_dicts(self, num_pats=2000):

        """ For each 11 slots we construct a set of
            <num_pats> patients. This group is split into 2 groups
            where the first half has taken a drug and the second half
            has not. In each iteration (1-11) we increase the survival
            percentage (w.r.t 1000 days) for the drug group and
            decrease the survival for the no drug group. Survivals for
            drug group are 0%, 10%, ... 100% percent. Inversely
            survivals for no drug group are 100%, 90%, ..., 0%). For
            pats that didn't survive 1000 days we set their survival
            arbitrarily to 999. Similarly, 1001 for the surviving pats.

            :param num_pats: number of patients to build in each subset
            :type num_pats: int
            :return patient_sets: a list of the patient_dicts
            :rtype patient_sets: list
        """

        patient_sets = []
        for i in range(0, 11):
            patient_i_dict = dict()
            set_divide = i/10

            # how many pats have surv >=
            group_1_divide = round((num_pats/2) * set_divide)
            # how many pats have surv <
            group_2_divide = round((num_pats/2) * (1-set_divide)) + (num_pats/2)

            if group_1_divide + group_2_divide != num_pats:
                raise Exception("pick num_pats where num_pats/2 is neatly partitions into 0%, 10%, ..., 100%")

            # First half is the drug group, second half is non-drug group. P(survival_time > X| gene) is 50%,
            # but the proportion shifts between the two groups based on i.
            for j in range(0, num_pats):
                # group 1 - drug group
                if j < int(num_pats/2):
                    if j >= group_1_divide:
                        survival_time = 999
                    else:
                        survival_time = 1001
                    patient_i_dict[j] = {'survival_time':survival_time,
                                         'gene_curies':tuple(['RAF1']),
                                         'drug_curies':tuple(['CYCLOPHOSPHAMIDE'])}
                # group 2 - non-drug group
                else:
                    if j >= group_2_divide:
                        survival_time = 999
                    else:
                        survival_time = 1001
                    patient_i_dict[j] = {'survival_time':survival_time,
                                      'gene_curies':tuple(['RAF1']),
                                      'drug_curies':tuple()}
            patient_sets.append(patient_i_dict)
        return patient_sets


if __name__ == '__main__':
    unittest.main()
