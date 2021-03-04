"""
    Source code developed by DI2AG.
    Thayer School of Engineering at Dartmouth College
    Authors:    Dr. Eugene Santos, Jr
                Mr. Chase Yakaboski,
                Mr. Gregory Hyde,
                Mr. Luke Veenhuis,
                Dr. Keum Joo Kim
"""

import requests
import json
import unittest

from chp_data.trapi_constants import *

from chp.trapi_interface import TrapiInterface
from reasoner_validator import validate_Message

class testWildCardHandler(unittest.TestCase):

    """ Unit tests for the wildcard_handler
        1. test_drug - queries Cyclophosphamide which has no sparsity issues
        2. test_sparse_drug - queries the drug Gemzar which has a sparsity issue.
               Should be handled by interpolation now.
        3. test_gene - queries TP53
        4. test with no gene/drug evidence but gene wildcard
        5. test with no gene/drug evidence but drug wildcard
        6. test with no gene/drug/survival evidence but with drug wildcard
        7. test illegal gene to disease edge
        8. test illegal chemical to disease edge
        9. test illegal disease to phenotype edge
        10. test illegal edge between gene and chemical
        11. test uknown edge
    """

    def test_drug(self):

        """ queries Cyclophosphamide which has no sparsity issue. Should return
            top 5 genes and a probability of survival
        """

        # empty response
        reasoner_std = { "query_graph": {}
                       }
        # empty query graph
        reasoner_std["query_graph"] = { "edges": {},
                                        "nodes": {}
                                      }

        # add in evidence drug
        drug = ('CYCLOPHOSPHAMIDE', 'CHEMBL:CHEMBL88')
        reasoner_std['query_graph']['nodes']['n0'] = {
                                                      'category':BIOLINK_DRUG,
                                                      'id':'{}'.format(drug[1])
                                                     }

        # add in gene node (to be filled by contribution analysis
        reasoner_std['query_graph']['nodes']['n1'] = {
                                                      'category':BIOLINK_GENE
                                                     }

        #add in disease node
        disease = ('Breast_Cancer', 'MONDO:0007254')
        reasoner_std['query_graph']['nodes']['n2'] = {
                                                      'category':BIOLINK_DISEASE,
                                                      'id':'{}'.format(disease[1])
                                                     }

        # add target survival node
        phenotype = ('Survival_Time', 'EFO:0000714')
        reasoner_std['query_graph']['nodes']['n3'] = {
                                                      'category': BIOLINK_PHENOTYPIC_FEATURE,
                                                      'id': '{}'.format(phenotype[1]),
                                                     }

        # link disease to target survival node
        reasoner_std['query_graph']['edges']['e0'] = {
                                                      'predicate':BIOLINK_DISEASE_TO_PHENOTYPIC_FEATURE_PREDICATE,
                                                      'subject':'n2',
                                                      'object':'n3',
                                                      'properties': {
                                                                     'qualifier': '>=',
                                                                     'days': 1000
                                                                    }
                                                     }

        # link genes/drugs to disease
        reasoner_std['query_graph']['edges']['e1'] = {
                                                      'predicate':BIOLINK_GENE_TO_DISEASE_PREDICATE,
                                                      'subject': 'n1',
                                                      'object': 'n2'
                                                     }
        reasoner_std['query_graph']['edges']['e2'] = {
                                                      'predicate':BIOLINK_CHEMICAL_TO_DISEASE_OR_PHENOTYPIC_FEATURE_PREDICATE,
                                                      'subject': 'n0',
                                                      'object': 'n2'
                                                     }

        # test input is TRAPI compliant
        validate_Message(reasoner_std) # doesn't return True/False for some reason... Will just present exception if not compliant

        handler = TrapiInterface(query=reasoner_std,
                                 max_results=2)
        queries = handler.build_chp_queries()
        queries = handler.run_chp_queries()
        reasoner_std_final = handler.construct_trapi_response()

        # test output is TRAPI compliant
        validate_Message(reasoner_std_final['message'])

        KG = reasoner_std_final["message"]['knowledge_graph']

        res_pretty = json.dumps(reasoner_std_final, indent=2)
        print(res_pretty)

        # extract probability
        for _, edge in KG['edges'].items():
            if edge['predicate'] == 'biolink:has_phenotype':
                p_survival = edge['attributes'][0]['value']
                break
        print("probability of survival:",p_survival)

    def test_sparse_drug(self):

        """ queries the drug Gemzar which has a sparsity issue. Used to throw an error.
            Should be handled.
        """

        # empty response
        reasoner_std = { "query_graph": {}
                       }

        # empty query graph
        reasoner_std["query_graph"] = { "edges": {},
                                        "nodes": {}
                                      }

        # add in evidence drug
        drug = ('GEMZAR', 'CHEMBL:CHEMBL888')
        reasoner_std['query_graph']['nodes']['n0'] = {
                                                      'category':BIOLINK_DRUG,
                                                      'id':'{}'.format(drug[1])
                                                     }

        # add in gene node (to be filled by contribution analysis
        reasoner_std['query_graph']['nodes']['n1'] = {
                                                      'category':BIOLINK_GENE
                                                     }

        #add in disease node
        disease = ('Breast_Cancer', 'MONDO:0007254')
        reasoner_std['query_graph']['nodes']['n2'] = {
                                                      'category':BIOLINK_DISEASE,
                                                      'id':'{}'.format(disease[1])
                                                     }

        # add target survival node
        phenotype = ('Survival_Time', 'EFO:0000714')
        reasoner_std['query_graph']['nodes']['n3'] = {
                                                      'category': BIOLINK_PHENOTYPIC_FEATURE,
                                                      'id': '{}'.format(phenotype[1]),
                                                     }

        # link disease to target survival node
        reasoner_std['query_graph']['edges']['e0'] = {
                                                      'predicate':BIOLINK_DISEASE_TO_PHENOTYPIC_FEATURE_PREDICATE,
                                                      'subject':'n2',
                                                      'object':'n3',
                                                      'properties': {
                                                                     'qualifier': '>=',
                                                                     'days': 940
                                                                    }
                                                     }

        # link genes/drugs to disease
        reasoner_std['query_graph']['edges']['e1'] = {
                                                      'predicate':BIOLINK_GENE_TO_DISEASE_PREDICATE,
                                                      'subject': 'n1',
                                                      'object': 'n2'
                                                     }
        reasoner_std['query_graph']['edges']['e2'] = {
                                                      'predicate':BIOLINK_CHEMICAL_TO_DISEASE_OR_PHENOTYPIC_FEATURE_PREDICATE,
                                                      'subject': 'n0',
                                                      'object': 'n2'
                                                     }

        # test input is TRAPI compliant
        validate_Message(reasoner_std) # doesn't return True/False for some reason... Will just present exception if not compliant

        handler = TrapiInterface(query=reasoner_std,
                                 max_results=2)
        queries = handler.build_chp_queries()
        queries = handler.run_chp_queries()
        reasoner_std_final = handler.construct_trapi_response()

        # test output is TRAPI compliant
        validate_Message(reasoner_std_final['message'])

        KG = reasoner_std_final["message"]['knowledge_graph']

        res_pretty = json.dumps(reasoner_std_final, indent=2)
        print(res_pretty)

        # extract probability
        for _, edge in KG['edges'].items():
            if edge['predicate'] == 'biolink:has_phenotype':
                p_survival = edge['attributes'][0]['value']
                break
        print("probability of survival:",p_survival)


    def test_gene(self):

        """ queries TP53. Should return 2 drugs.
        """

        # empty response
        reasoner_std = { "query_graph": dict()}
        # empty query graph
        reasoner_std["query_graph"] = { "edges": {},
                                        "nodes": {}
                                      }

        # add in evidence drug
        reasoner_std['query_graph']['nodes']['n0'] = {
                                                      'category':BIOLINK_DRUG
                                                     }

        # add in gene node (to be filled by contribution analysis
        reasoner_std['query_graph']['nodes']['n1'] = {
                                                      'category':BIOLINK_GENE,
                                                      'id': 'ENSEMBL:ENSG00000141510'
                                                     }

        #add in disease node
        disease = ('Breast_Cancer', 'MONDO:0007254')
        reasoner_std['query_graph']['nodes']['n2'] = {
                                                      'category':BIOLINK_DISEASE,
                                                      'id':'{}'.format(disease[1])
                                                     }

        # add target survival node
        phenotype = ('Survival_Time', 'EFO:0000714')
        reasoner_std['query_graph']['nodes']['n3'] = {
                                                      'category': BIOLINK_PHENOTYPIC_FEATURE,
                                                      'id': '{}'.format(phenotype[1]),
                                                     }

        # link disease to target survival node
        reasoner_std['query_graph']['edges']['e0'] = {
                                                      'predicate':BIOLINK_DISEASE_TO_PHENOTYPIC_FEATURE_PREDICATE,
                                                      'subject':'n2',
                                                      'object':'n3',
                                                      'properties': {
                                                                     'qualifier': '>=',
                                                                     'days': 940
                                                                    }
                                                     }

        # link genes/drugs to disease
        reasoner_std['query_graph']['edges']['e1'] = {
                                                      'predicate':BIOLINK_GENE_TO_DISEASE_PREDICATE,
                                                      'subject': 'n1',
                                                      'object': 'n2'
                                                     }
        reasoner_std['query_graph']['edges']['e2'] = {
                                                      'predicate':BIOLINK_CHEMICAL_TO_DISEASE_OR_PHENOTYPIC_FEATURE_PREDICATE,
                                                      'subject': 'n0',
                                                      'object': 'n2'
                                                     }

        # test input is TRAPI compliant
        validate_Message(reasoner_std) # doesn't return True/False for some reason... Will just present exception if not compliant

        handler = TrapiInterface(query=reasoner_std,
                                 max_results=2)
        queries = handler.build_chp_queries()
        queries = handler.run_chp_queries()
        reasoner_std_final = handler.construct_trapi_response()

        # test output is TRAPI compliant
        validate_Message(reasoner_std_final['message'])

        KG = reasoner_std_final["message"]['knowledge_graph']

        res_pretty = json.dumps(reasoner_std_final, indent=2)
        print(res_pretty)

        # extract probability
        for _, edge in KG['edges'].items():
            if edge['predicate'] == 'biolink:has_phenotype':
                p_survival = edge['attributes'][0]['value']
                break
        print("probability of survival:",p_survival)

    def test_no_gene_drug_evidence_with_gene_wildcard(self):

        """ queries with no gene/drug evidence, but with gene wildcard
        """

        # empty response
        reasoner_std = { "query_graph": {}
                       }

        # empty query graph
        reasoner_std["query_graph"] = { "edges": {},
                                        "nodes": {}
                                      }

        # add in gene node (to be filled by contribution analysis
        reasoner_std['query_graph']['nodes']['n0'] = {
                                                      'category':BIOLINK_GENE
                                                     }

        #add in disease node
        disease = ('Breast_Cancer', 'MONDO:0007254')
        reasoner_std['query_graph']['nodes']['n1'] = {
                                                      'category':BIOLINK_DISEASE,
                                                      'id':'{}'.format(disease[1])
                                                     }

        # add target survival node
        phenotype = ('Survival_Time', 'EFO:0000714')
        reasoner_std['query_graph']['nodes']['n2'] = {
                                                      'category': BIOLINK_PHENOTYPIC_FEATURE,
                                                      'id': '{}'.format(phenotype[1]),
                                                     }

        # link disease to target survival node
        reasoner_std['query_graph']['edges']['e0'] = {
                                                      'predicate':BIOLINK_DISEASE_TO_PHENOTYPIC_FEATURE_PREDICATE,
                                                      'subject':'n1',
                                                      'object':'n2',
                                                      'properties': {
                                                                     'qualifier': '>=',
                                                                     'days': 940
                                                                    }
                                                     }

        # link genes/drugs to disease
        reasoner_std['query_graph']['edges']['e1'] = {
                                                      'predicate':BIOLINK_GENE_TO_DISEASE_PREDICATE,
                                                      'subject': 'n0',
                                                      'object': 'n1'
                                                     }

        # test input is TRAPI compliant
        validate_Message(reasoner_std) # doesn't return True/False for some reason... Will just present exception if not compliant

        handler = TrapiInterface(query=reasoner_std,
                                 max_results=2)
        queries = handler.build_chp_queries()
        queries = handler.run_chp_queries()
        reasoner_std_final = handler.construct_trapi_response()

        # test output is TRAPI compliant
        validate_Message(reasoner_std_final['message'])

        KG = reasoner_std_final["message"]['knowledge_graph']

        res_pretty = json.dumps(reasoner_std_final, indent=2)
        print(res_pretty)

        # extract probability
        for _, edge in KG['edges'].items():
            if edge['predicate'] == 'biolink:has_phenotype':
                p_survival = edge['attributes'][0]['value']
                break
        print("probability of survival:",p_survival)

    def test_no_gene_drug_evidence_with_drug_wildcard(self):

        """ queries with no gene/drug evidence, but with drug wildcard
        """

        # empty response
        reasoner_std = { "query_graph": {}
                       }

        # empty query graph
        reasoner_std["query_graph"] = { "edges": {},
                                        "nodes": {}
                                      }

        # add in gene node (to be filled by contribution analysis
        reasoner_std['query_graph']['nodes']['n0'] = {
                                                      'category':BIOLINK_DRUG
                                                     }

        #add in disease node
        disease = ('Breast_Cancer', 'MONDO:0007254')
        reasoner_std['query_graph']['nodes']['n1'] = {
                                                      'category':BIOLINK_DISEASE,
                                                      'id':'{}'.format(disease[1])
                                                     }

        # add target survival node
        phenotype = ('Survival_Time', 'EFO:0000714')
        reasoner_std['query_graph']['nodes']['n2'] = {
                                                      'category': BIOLINK_PHENOTYPIC_FEATURE,
                                                      'id': '{}'.format(phenotype[1]),
                                                     }

        # link disease to target survival node
        reasoner_std['query_graph']['edges']['e0'] = {
                                                      'predicate':BIOLINK_DISEASE_TO_PHENOTYPIC_FEATURE_PREDICATE,
                                                      'subject':'n1',
                                                      'object':'n2',
                                                      'properties': {
                                                                     'qualifier': '>=',
                                                                     'days': 1000
                                                                    }
                                                     }

        # link genes/drugs to disease
        reasoner_std['query_graph']['edges']['e1'] = {
                                                      'predicate':BIOLINK_CHEMICAL_TO_DISEASE_OR_PHENOTYPIC_FEATURE_PREDICATE,
                                                      'subject': 'n0',
                                                      'object': 'n1'
                                                     }

        # test input is TRAPI compliant
        validate_Message(reasoner_std) # doesn't return True/False for some reason... Will just present exception if not compliant

        handler = TrapiInterface(query=reasoner_std,
                                 max_results=2)
        queries = handler.build_chp_queries()
        queries = handler.run_chp_queries()
        reasoner_std_final = handler.construct_trapi_response()

        # test output is TRAPI compliant
        validate_Message(reasoner_std_final['message'])

        KG = reasoner_std_final["message"]['knowledge_graph']

        res_pretty = json.dumps(reasoner_std_final, indent=2)
        print(res_pretty)

        # extract probability
        for _, edge in KG['edges'].items():
            if edge['predicate'] == 'biolink:has_phenotype':
                p_survival = edge['attributes'][0]['value']
                break
        print("probability of survival:",p_survival)

    def test_no_gene_drug_phenotypic_evidence_with_drug_wildcard(self):

        """ queries with no gene/drug/survival evidence, but with drug wildcard
        """

        # empty response
        reasoner_std = { "query_graph": {}
                       }

        # empty query graph
        reasoner_std["query_graph"] = { "edges": {},
                                        "nodes": {}
                                      }

        # add in gene node (to be filled by contribution analysis
        reasoner_std['query_graph']['nodes']['n0'] = {
                                                      'category':BIOLINK_DRUG
                                                     }

        #add in disease node
        disease = ('Breast_Cancer', 'MONDO:0007254')
        reasoner_std['query_graph']['nodes']['n1'] = {
                                                      'category':BIOLINK_DISEASE,
                                                      'id':'{}'.format(disease[1])
                                                     }

        # link genes/drugs to disease
        reasoner_std['query_graph']['edges']['e1'] = {
                                                      'predicate':BIOLINK_CHEMICAL_TO_DISEASE_OR_PHENOTYPIC_FEATURE_PREDICATE,
                                                      'subject': 'n0',
                                                      'object': 'n1'
                                                     }

        # test input is TRAPI compliant
        validate_Message(reasoner_std) # doesn't return True/False for some reason... Will just present exception if not compliant

        handler = TrapiInterface(query=reasoner_std,
                                 max_results=2)
        queries = handler.build_chp_queries()
        queries = handler.run_chp_queries()
        reasoner_std_final = handler.construct_trapi_response()

        # test output is TRAPI compliant
        validate_Message(reasoner_std_final['message'])

        KG = reasoner_std_final["message"]['knowledge_graph']

        res_pretty = json.dumps(reasoner_std_final, indent=2)
        print(res_pretty)

    def test_illegal_gene_to_disease(self):
        """ link between gene and disease node is backwards. Should return appropriate error message.
        """
        # empty response
        reasoner_std = { "query_graph": {}
                       }
        # empty query graph
        reasoner_std["query_graph"] = { "edges": {},
                                        "nodes": {}
                                      }

        # add in evidence drug
        drug = ('CYCLOPHOSPHAMIDE', 'CHEMBL:CHEMBL88')
        reasoner_std['query_graph']['nodes']['n0'] = {
                                                      'category':BIOLINK_DRUG,
                                                      'id':'{}'.format(drug[1])
                                                     }

        # add in gene node (to be filled by contribution analysis
        reasoner_std['query_graph']['nodes']['n1'] = {
                                                      'category':BIOLINK_GENE
                                                     }

        #add in disease node
        disease = ('Breast_Cancer', 'MONDO:0007254')
        reasoner_std['query_graph']['nodes']['n2'] = {
                                                      'category':BIOLINK_DISEASE,
                                                      'id':'{}'.format(disease[1])
                                                     }

        # add target survival node
        phenotype = ('Survival_Time', 'EFO:0000714')
        reasoner_std['query_graph']['nodes']['n3'] = {
                                                      'category': BIOLINK_PHENOTYPIC_FEATURE,
                                                      'id': '{}'.format(phenotype[1]),
                                                     }

        # link disease to target survival node
        reasoner_std['query_graph']['edges']['e0'] = {
                                                      'predicate':BIOLINK_DISEASE_TO_PHENOTYPIC_FEATURE_PREDICATE,
                                                      'subject':'n2',
                                                      'object':'n3',
                                                      'properties': {
                                                                     'qualifier': '>=',
                                                                     'days': 1000
                                                                    }
                                                     }

        # link genes/drugs to disease
        reasoner_std['query_graph']['edges']['e1'] = {
                                                      'predicate':BIOLINK_GENE_TO_DISEASE_PREDICATE,
                                                      'subject': 'n2',
                                                      'object': 'n1'
                                                     }
        reasoner_std['query_graph']['edges']['e2'] = {
                                                      'predicate':BIOLINK_CHEMICAL_TO_DISEASE_OR_PHENOTYPIC_FEATURE_PREDICATE,
                                                      'subject': 'n0',
                                                      'object': 'n2'
                                                     }

        # test input is TRAPI compliant
        validate_Message(reasoner_std) # doesn't return True/False for some reason... Will just present exception if not compliant

        handler = TrapiInterface(query=reasoner_std,
                                 max_results=2)
        queries = handler.build_chp_queries()
        self.assertEqual('Wildcard type detected. Gene to disease edge (edge id: e1) malformed.', handler.error_msg)

    def test_illegal_chemical_to_disease(self):
        """ link between drug and disease is backwards. Should return appropriate error message.
        """
        # empty response
        reasoner_std = { "query_graph": {}
                       }
        # empty query graph
        reasoner_std["query_graph"] = { "edges": {},
                                        "nodes": {}
                                      }

        # add in evidence drug
        drug = ('CYCLOPHOSPHAMIDE', 'CHEMBL:CHEMBL88')
        reasoner_std['query_graph']['nodes']['n0'] = {
                                                      'category':BIOLINK_DRUG,
                                                      'id':'{}'.format(drug[1])
                                                     }

        # add in gene node (to be filled by contribution analysis
        reasoner_std['query_graph']['nodes']['n1'] = {
                                                      'category':BIOLINK_GENE
                                                     }

        #add in disease node
        disease = ('Breast_Cancer', 'MONDO:0007254')
        reasoner_std['query_graph']['nodes']['n2'] = {
                                                      'category':BIOLINK_DISEASE,
                                                      'id':'{}'.format(disease[1])
                                                     }

        # add target survival node
        phenotype = ('Survival_Time', 'EFO:0000714')
        reasoner_std['query_graph']['nodes']['n3'] = {
                                                      'category': BIOLINK_PHENOTYPIC_FEATURE,
                                                      'id': '{}'.format(phenotype[1]),
                                                     }

        # link disease to target survival node
        reasoner_std['query_graph']['edges']['e0'] = {
                                                      'predicate':BIOLINK_DISEASE_TO_PHENOTYPIC_FEATURE_PREDICATE,
                                                      'subject':'n2',
                                                      'object':'n3',
                                                      'properties': {
                                                                     'qualifier': '>=',
                                                                     'days': 1000
                                                                    }
                                                     }
        # link genes/drugs to disease
        reasoner_std['query_graph']['edges']['e1'] = {
                                                      'predicate':BIOLINK_GENE_TO_DISEASE_PREDICATE,
                                                      'subject': 'n1',
                                                      'object': 'n2'
                                                     }
        reasoner_std['query_graph']['edges']['e2'] = {
                                                      'predicate':BIOLINK_CHEMICAL_TO_DISEASE_OR_PHENOTYPIC_FEATURE_PREDICATE,
                                                      'subject': 'n2',
                                                      'object': 'n0'
                                                     }

        # test input is TRAPI compliant
        validate_Message(reasoner_std) # doesn't return True/False for some reason... Will just present exception if not compliant

        handler = TrapiInterface(query=reasoner_std,
                                 max_results=2)
        queries = handler.build_chp_queries()
        self.assertEqual('Wildcard type detected. Chemical to disease edge (edge id: e2) malformed.', handler.error_msg)

    def test_illegal_disease_to_phenotype(self):
        """ link between disease to phenotype edge is backwards. Should return appropriate error message.
        """
        # empty response
        reasoner_std = { "query_graph": {}
                       }
        # empty query graph
        reasoner_std["query_graph"] = { "edges": {},
                                        "nodes": {}
                                      }

        # add in evidence drug
        drug = ('CYCLOPHOSPHAMIDE', 'CHEMBL:CHEMBL88')
        reasoner_std['query_graph']['nodes']['n0'] = {
                                                      'category':BIOLINK_DRUG,
                                                      'id':'{}'.format(drug[1])
                                                     }

        # add in gene node (to be filled by contribution analysis
        reasoner_std['query_graph']['nodes']['n1'] = {
                                                      'category':BIOLINK_GENE
                                                     }

        #add in disease node
        disease = ('Breast_Cancer', 'MONDO:0007254')
        reasoner_std['query_graph']['nodes']['n2'] = {
                                                      'category':BIOLINK_DISEASE,
                                                      'id':'{}'.format(disease[1])
                                                     }

        # add target survival node
        phenotype = ('Survival_Time', 'EFO:0000714')
        reasoner_std['query_graph']['nodes']['n3'] = {
                                                      'category': BIOLINK_PHENOTYPIC_FEATURE,
                                                      'id': '{}'.format(phenotype[1]),
                                                     }

        # link disease to target survival node
        reasoner_std['query_graph']['edges']['e0'] = {
                                                      'predicate':BIOLINK_DISEASE_TO_PHENOTYPIC_FEATURE_PREDICATE,
                                                      'subject':'n2',
                                                      'object':'n3',
                                                      'properties': {
                                                                     'qualifier': '>=',
                                                                     'days': 1000
                                                                    }
                                                     }

        # link genes/drugs to disease
        reasoner_std['query_graph']['edges']['e1'] = {
                                                      'predicate':BIOLINK_GENE_TO_DISEASE_PREDICATE,
                                                      'subject': 'n1',
                                                      'object': 'n2'
                                                     }
        reasoner_std['query_graph']['edges']['e2'] = {
                                                      'predicate':BIOLINK_CHEMICAL_TO_DISEASE_OR_PHENOTYPIC_FEATURE_PREDICATE,
                                                      'subject': 'n0',
                                                      'object': 'n2'
                                                     }
        reasoner_std['query_graph']['edges']['e3'] = {
                                                      'predicate':BIOLINK_DISEASE_TO_PHENOTYPIC_FEATURE_PREDICATE,
                                                      'subject': 'n3',
                                                      'object': 'n2'
                                                     }

        # test input is TRAPI compliant
        validate_Message(reasoner_std) # doesn't return True/False for some reason... Will just present exception if not compliant

        handler = TrapiInterface(query=reasoner_std,
                                 max_results=2)
        queries = handler.build_chp_queries()
        self.assertEqual('Wildcard type detected. Disease to phenotype edge (edge id: e3) malformed.', handler.error_msg)

    def test_illegal_edge(self):
        """ uses illegal gene-drug edge that is not compatible. Should return respective error message.
        """
        # empty response
        reasoner_std = { "query_graph": {}
                       }
        # empty query graph
        reasoner_std["query_graph"] = { "edges": {},
                                        "nodes": {}
                                      }

        # add in evidence drug
        drug = ('CYCLOPHOSPHAMIDE', 'CHEMBL:CHEMBL88')
        reasoner_std['query_graph']['nodes']['n0'] = {
                                                      'category':BIOLINK_DRUG,
                                                      'id':'{}'.format(drug[1])
                                                     }

        # add in gene node (to be filled by contribution analysis
        reasoner_std['query_graph']['nodes']['n1'] = {
                                                      'category':BIOLINK_GENE
                                                     }

        #add in disease node
        disease = ('Breast_Cancer', 'MONDO:0007254')
        reasoner_std['query_graph']['nodes']['n2'] = {
                                                      'category':BIOLINK_DISEASE,
                                                      'id':'{}'.format(disease[1])
                                                     }

        # add target survival node
        phenotype = ('Survival_Time', 'EFO:0000714')
        reasoner_std['query_graph']['nodes']['n3'] = {
                                                      'category': BIOLINK_PHENOTYPIC_FEATURE,
                                                      'id': '{}'.format(phenotype[1]),
                                                     }

        # link disease to target survival node
        reasoner_std['query_graph']['edges']['e0'] = {
                                                      'predicate':BIOLINK_DISEASE_TO_PHENOTYPIC_FEATURE_PREDICATE,
                                                      'subject':'n2',
                                                      'object':'n3',
                                                      'properties': {
                                                                     'qualifier': '>=',
                                                                     'days': 1000
                                                                    }
                                                     }

        # link genes/drugs to disease
        reasoner_std['query_graph']['edges']['e1'] = {
                                                      'predicate':BIOLINK_GENE_TO_DISEASE_PREDICATE,
                                                      'subject': 'n1',
                                                      'object': 'n2'
                                                     }
        reasoner_std['query_graph']['edges']['e2'] = {
                                                      'predicate':BIOLINK_CHEMICAL_TO_DISEASE_OR_PHENOTYPIC_FEATURE_PREDICATE,
                                                      'subject': 'n0',
                                                      'object': 'n2'
                                                     }
        reasoner_std['query_graph']['edges']['e3'] = {
                                                      'predicate':BIOLINK_DISEASE_TO_PHENOTYPIC_FEATURE_PREDICATE,
                                                      'subject': 'n2',
                                                      'object': 'n3'
                                                     }

        reasoner_std['query_graph']['edges']['e4'] = {
                                                      'predicate':BIOLINK_CHEMICAL_TO_GENE_PREDICATE,
                                                      'subject': 'n0',
                                                      'object': 'n1'
                                                     }

        # test input is TRAPI compliant
        validate_Message(reasoner_std) # doesn't return True/False for some reason... Will just present exception if not compliant

        handler = TrapiInterface(query=reasoner_std,
                                 max_results=2)
        queries = handler.build_chp_queries()
        self.assertEqual('Wildcard type detected. Received edge between gene and chemical (edge id: e4). This edge is incompatible with this wildcard type.', handler.error_msg)

    def test_unknown_edge(self):
        """ uses unknown edge. Should return respective error message.
        """
        # empty response
        reasoner_std = { "query_graph": {}
                       }
        # empty query graph
        reasoner_std["query_graph"] = { "edges": {},
                                        "nodes": {}
                                      }

        # add in evidence drug
        drug = ('CYCLOPHOSPHAMIDE', 'CHEMBL:CHEMBL88')
        reasoner_std['query_graph']['nodes']['n0'] = {
                                                      'category':BIOLINK_DRUG,
                                                      'id':'{}'.format(drug[1])
                                                     }

        # add in gene node (to be filled by contribution analysis
        reasoner_std['query_graph']['nodes']['n1'] = {
                                                      'category':BIOLINK_GENE
                                                     }

        #add in disease node
        disease = ('Breast_Cancer', 'MONDO:0007254')
        reasoner_std['query_graph']['nodes']['n2'] = {
                                                      'category':BIOLINK_DISEASE,
                                                      'id':'{}'.format(disease[1])
                                                     }

        # add target survival node
        phenotype = ('Survival_Time', 'EFO:0000714')
        reasoner_std['query_graph']['nodes']['n3'] = {
                                                      'category': BIOLINK_PHENOTYPIC_FEATURE,
                                                      'id': '{}'.format(phenotype[1]),
                                                     }

        # link disease to target survival node
        reasoner_std['query_graph']['edges']['e0'] = {
                                                      'predicate':BIOLINK_DISEASE_TO_PHENOTYPIC_FEATURE_PREDICATE,
                                                      'subject':'n2',
                                                      'object':'n3',
                                                      'properties': {
                                                                     'qualifier': '>=',
                                                                     'days': 1000
                                                                    }
                                                     }

        # link genes/drugs to disease
        reasoner_std['query_graph']['edges']['e1'] = {
                                                      'predicate':BIOLINK_GENE_TO_DISEASE_PREDICATE,
                                                      'subject': 'n1',
                                                      'object': 'n2'
                                                     }
        reasoner_std['query_graph']['edges']['e2'] = {
                                                      'predicate':BIOLINK_CHEMICAL_TO_DISEASE_OR_PHENOTYPIC_FEATURE_PREDICATE,
                                                      'subject': 'n0',
                                                      'object': 'n2'
                                                     }
        reasoner_std['query_graph']['edges']['e3'] = {
                                                      'predicate':BIOLINK_DISEASE_TO_PHENOTYPIC_FEATURE_PREDICATE,
                                                      'subject': 'n2',
                                                      'object': 'n3'
                                                     }

        reasoner_std['query_graph']['edges']['e4'] = {
                                                      'predicate':'biolink:unknown_edge',
                                                      'subject': 'n0',
                                                      'object': 'n1'
                                                     }

        # test input is TRAPI compliant
        validate_Message(reasoner_std) # doesn't return True/False for some reason... Will just present exception if not compliant

        handler = TrapiInterface(query=reasoner_std,
                                 max_results=2)
        queries = handler.build_chp_queries()
        self.assertEqual('Wildcard type detected. Unknown predicate type for edge (edge id: e4).', handler.error_msg)

if __name__ == '__main__':
    unittest.main()
