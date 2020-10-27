import unittest
import copy
from pybkb.common.bayesianKnowledgeBase import bayesianKnowledgeBase as BKB
from pybkb.common.bayesianKnowledgeBase import BKB_S_node, BKB_component, BKB_I_node
from pybkb.python_base.fusion import fuse as py_fuse
from chp.reasoner import Reasoner
from chp.query import Query

class testSyntheticPatients(unittest.TestCase):

    ##########################################################
    # test_synthetic
    # Input:
    # Output:
    #--------------------------------------------------------
    # Description: Builds a set of 11 patient sets all with a single
    # gene mutation. In each sub set of patients the set is split
    # in two. Half receive some drug treatment and half do not.
    # Each preceting set of the 11 has a stronger drug treatment
    # signal. That is, survival percentages (of patients >= 1000 days)
    # for set 1 with drug treatment is 0% and without is 100%.
    # For set two that is 10% and 90% respectively, until set 11 where
    # they are flipped to be 100% and 0% respectively. We assert the
    # following:
    # for sets 1-5 - P(survival_time >= 1000 days | gene) >
    #                P(survival_time >= 1000 days | gene & drug)
    # for set 6 - P(survival_time >= 1000 days | gene) ==
    #             P(survival_time >= 1000 days | gene & drug)
    # for set 7-11 - P(survival_time >= 1000 days | gene) <
    #                P(survival_time >= 1000 days | gene & drug)
    #
    # Note - P(survival_time >= 1000 days | gene) is always 50%
    # Note - P(survival_time >= 1000 days | gene & drug) is
    #        0%, 10%, ..., 100% for sets 1-11.

    def test_synthetic(self):
        patient_sets = self.build_synthetic_sets()
        fragments = self.build_synthetic_BKB_frags(patient_sets)
        dicts, source_names, reliabilities = self.build_synthetic_dicts(patient_sets)
        fused_bkbs = self.fuse_synthetic_frags(fragments, source_names, reliabilities)
        query_responses = self.run_synthetic_queries(fused_bkbs, dicts)
        for i in range(0, len(query_responses)):
            if i < 5:
                self.assertGreater(query_responses[i]['no_drug'], query_responses[i]['drug'])
            elif i == 5:
                # to 7 decimal places
                self.assertAlmostEqual(query_responses[i]['no_drug'], query_responses[i]['drug'], 7)
            else:
                self.assertGreater(query_responses[i]['drug'], query_responses[i]['no_drug'])

    ##########################################################
    # run_synthetic_queries
    # Input: fused_bkbs, dicts
    # Output: probability of survivals for each patient set
    #--------------------------------------------------------
    # Description: Each patient set two queries must be run:
    #    1.) P(survival_time >= 1000 | gene & drug)
    #    2.) P(survival_time >= 1000 | gene)
    # The probabilities for each set are returned in a list of
    # dictionaries indexed by 'drug' and 'no_drug'.

    def run_synthetic_queries(self, fused_bkbs, dicts):
        query_responses = []
        for i in range(0, 11):
            # with drug
            reasoner = Reasoner(fused_bkb=fused_bkbs[i], patient_data=dicts[i])
            gene_evidence = {'mut_RAF1' : 'True'}
            query_drug = Query(evidence=gene_evidence,
                               targets=list(),
                               meta_evidence=[('Drug_Name(s)', '==','CYCLOPHOSPHAMIDE')],
                               meta_targets=[('Survival_Time', '>=', 1000)],
                               type='updating')
            query_drug1 = reasoner.analyze_query(copy.deepcopy(query_drug),
                                                save_dir=None,
                                                target_strategy='explicit',
                                                interpolation='standard')
            truth_assignment_drug = self.get_synthetic_query_probs(query_drug1)

            # with no drug
            reasoner = Reasoner(fused_bkb=fused_bkbs[i], patient_data=dicts[i])
            query_no_drug = Query(evidence=gene_evidence,
                                  targets=list(),
                                  meta_evidence=None,
                                  meta_targets=[('Survival_Time', '>=', 1000)],
                                  type='updating')
            query_no_drug = reasoner.analyze_query(query_no_drug,
                                                   save_dir=None,
                                                   target_strategy='explicit',
                                                   interpolation='standard')
            truth_assignment_no_drug = self.get_synthetic_query_probs(query_no_drug)

            # answers
            query_responses.append({'drug':truth_assignment_drug,
                                    'no_drug':truth_assignment_no_drug})
        return query_responses

    ##########################################################
    # get_synthetic_query_probs
    # Input: finished query
    # Output: True probabilties from a finished query
    #--------------------------------------------------------
    # Description: indexes and a run query's update to extract
    # end probabilities. These are then normalized. Only the
    # True probability is returned

    def get_synthetic_query_probs(self, query):
        target_info = []
        for update, prob in query.result.updates.items():
            comp_idx, state_idx = update
            comp_name = query.bkb.getComponentName(comp_idx)
            state_name = query.bkb.getComponentINodeName(comp_idx, state_idx)
            target_info.append([comp_name, state_name, prob])

        if target_info[0][1] == 'True':
            truth_assignment = target_info[0][2]
            false_assignment = target_info[1][2]
        else:
            truth_assignment = target_info[1][2]
            false_assignment = target_info[0][2]

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
        return truth_assignment

    ##########################################################
    # fuse_synthetic_frags
    # Input: bkf fragments, source names and reliabilities
    # Output: set of 11 fused bkbs for each 11 patient sets
    #--------------------------------------------------------
    # Description: Runs bkb fusion on each of the 11 patient
    # sets. Returned is a fused BKB for each set.

    def fuse_synthetic_frags(self, fragments, source_names, reliabilities):
        fused_bkbs = []
        for i in range(0, 11):
            fused_bkb = py_fuse(fragments[i],
                                reliabilities = reliabilities[i],
                                source_names = source_names[i])
            fused_bkbs.append(fused_bkb)
        return fused_bkbs

    ##########################################################
    # build_synthetic_dicts
    # Input: each of the 11 patient sets
    # Output: set of 11 patient_dicts, set of 11 source sets
    #         and a set of 11 reliabilities sets.
    #--------------------------------------------------------
    # Description: for each set in the patient sets a dictionary
    # is made incorporating all patient drug (pat[3]) and
    # survival_time (pat[1]) info. Additionally we collect the
    # patient ids (pat[0]) and build out reliabilties for each
    # fragment (these are just 1).

    def build_synthetic_dicts(self, patient_sets):
        pat_dicts = []
        source_sets = []
        reliability_sets = []
        for i in range(0, 11):
            set_dict = dict()
            source_names = []
            reliabilities = []
            for pat in patient_sets[i]:
                # has drug
                if len(pat) == 4:
                    drug_tuple = tuple()
                    set_dict[pat[0]] = { 'Drug_Name(s)':tuple([pat[3]]),
                                         'Survival_Time':pat[1]}
                # has no drug
                else:
                    set_dict[pat[0]] = { 'Drug_Name(s)':tuple(),
                                         'Survival_Time':pat[1]}
                source_names.append(str(pat[0]))
                reliabilities.append(1)
            pat_dicts.append(set_dict)
            source_sets.append(source_names)
            reliability_sets.append(reliabilities)
        return pat_dicts, source_sets, reliability_sets

     ##########################################################
    # build_synthetic_BKB_frags
    # Input: set of 11 patient sets
    # Output: set of 11 patient frags
    #--------------------------------------------------------
    # Description: For each patient in each 11 patient sets we
    # construct a BKF. This is just a simple link of the form:
    #     1.0 o---->[mut_RAF1 = True]

    def build_synthetic_BKB_frags(self, patient_sets):
        bkf_sets = []
        for i in range(0, 11):
            bkf_set = []
            for pat in patient_sets[i]:
                bkf = BKB(name = pat[0])
                # gene
                mutGeneComp_idx = bkf.addComponent('mut_{}'.format(pat[2]))
                iNodeGeneMut_idx = bkf.addComponentState(mutGeneComp_idx, 'True')
                # form SNode  o---->[mut_<genename> = True]
                bkf.addSNode(BKB_S_node(mutGeneComp_idx, iNodeGeneMut_idx, 1.0))
                bkf_set.append(bkf)
            bkf_sets.append(bkf_set)
        return bkf_sets

    ##########################################################
    # build_synthetic_sets
    # Input: number of synthetic patients
    # Output: set of 11 patient sets
    #--------------------------------------------------------
    # Description: for each 11 slots we construct a set of
    # <num_pats> patients. This group is split into 2 groups
    # where the first half has taken a drug and the second half
    # has not. In each iteration (1-11) we increase the survival
    # percentage (w.r.t 1000 days) for the drug group and
    # decrease the survival for the no drug group. Survivals for
    # drug group are 0%, 10%, ... 100% percent. Inversely
    # survivals for no drug group are 100%, 90%, ..., 0%). For
    # pats that didn't survive 1000 days we set their survival
    # arbitrarily to 999. Similarly, 1001 for the surviving pats.

    def build_synthetic_sets(self, num_pats=2000):
        patient_sets = []
        for i in range(0, 11):
            patient_i_set = []
            set_divide = i/10

            # how many pats have surv >=
            group_1_divide = round((num_pats/2) * set_divide)
            # how many pats have surv <
            group_2_divide = round((num_pats/2) * (1-set_divide)) + (num_pats/2)

            if group_1_divide + group_2_divide != num_pats:
                raise Exception("pick num_pats where num_pats/2 is neatly partitions into 0%, 10%, ..., 100%")

            for j in range(0, num_pats):
                # group 1 - drug group
                if j < int(num_pats/2):
                    if j >= group_1_divide:
                        survival_time = 999
                    else:
                        survival_time = 1001
                    patient_i_set.append((j,survival_time,'RAF1', 'CYCLOPHOSPHAMIDE'))
                # group 2 - non-drug group
                else:
                    if j >= group_2_divide:
                        survival_time = 999
                    else:
                        survival_time = 1001
                    patient_i_set.append((j,survival_time,'RAF1'))
            patient_sets.append(patient_i_set)
        return patient_sets


if __name__ == '__main__':
    unittest.main()
