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

class testDefaultHandler(unittest.TestCase):

    """ Unit tests for default_handler. Uses 8 test cases:
        1. Testing negative result normalization
        2. Normal request with two genes and a drug
        3. Test with no evidence
        4. Test with no evidence, but omitting KG and Results (should be handled by handler)
        5. Test with one gene
        6. Test with one drug
        7. Test with no target (should crash)
        8. Test with no disease
        9. Test default survival
       10. Test illegal gene to disease
       11. Test illegal drug to disease
       12. Test illegal disease to phenotype
       13. Test illegal edge from gene to drug
       14. Test unknown edge
    """

    def test_negative(self):

        """ Testing negative result normalization
        """

        # empty response
        reasoner_std = { "query_graph": dict()
                       }

        # empty query graph
        reasoner_std["query_graph"] = { "edges": dict(),
                                        "nodes": dict()
                                      }
        # add in evidence gene
        gene1 = ('X', 'ENSEMBL:ENSG00000007174')
        reasoner_std['query_graph']['nodes']['n0'] = { 'category':BIOLINK_GENE,
                                                       'id':'{}'.format(gene1[1])
                                                     }
        # add in evidence drug
        drug = ('Y', 'CHEMBL:CHEMBL554')
        reasoner_std['query_graph']['nodes']['n1'] = { 'category':BIOLINK_DRUG,
                                                       'id':'{}'.format(drug[1])
                                                     }
        # add in disease node
        disease = ('Breast_Cancer', 'MONDO:0007254')
        reasoner_std['query_graph']['nodes']['n2'] = { 'category':BIOLINK_DISEASE,
                                                       'id':'{}'.format(disease[1])
                                                     }
        # add target survival node
        phenotype = ('Survival_Time', 'EFO:0000714')
        reasoner_std['query_graph']['nodes']['n3'] = { 'category': BIOLINK_PHENOTYPIC_FEATURE,
                                                       'id': '{}'.format(phenotype[1]),
                                                     }
        # link genes/drugs to disease
        reasoner_std['query_graph']['edges']['e0'] = { 'predicate':BIOLINK_GENE_TO_DISEASE_PREDICATE,
                                                       'subject': 'n0',
                                                       'object': 'n2'
                                                     }
        reasoner_std['query_graph']['edges']['e1'] = { 'predicate':BIOLINK_CHEMICAL_TO_DISEASE_OR_PHENOTYPIC_FEATURE_PREDICATE,
                                                       'subject': 'n1',
                                                       'object': 'n2'
                                                     }
        # link disease to target
        reasoner_std['query_graph']['edges']['e2'] = { 'predicate':BIOLINK_DISEASE_TO_PHENOTYPIC_FEATURE_PREDICATE,
                                                       'subject': 'n2',
                                                       'object': 'n3',
                                                       'properties': { 'qualifier':'>=',
                                                                       'days':970
                                                                     }
                                                     }

        # test input is TRAPI compliant
        validate_Message(reasoner_std) # doesn't return True/False for some reason... Will just present exception if not compliant

        handler = TrapiInterface(query=reasoner_std)
        queries = handler.build_chp_queries()
        queries = handler.run_chp_queries()
        reasoner_std_final = handler.construct_trapi_response()

        # test output is TRAPI compliant
        validate_Message(reasoner_std_final['message'])

        KG = reasoner_std_final['message']['knowledge_graph']
        for edge_key in KG['edges'].keys():
            edge = KG['edges'][edge_key]
            if edge['predicate'] == BIOLINK_DISEASE_TO_PHENOTYPIC_FEATURE_PREDICATE:
                p_survival = edge['attributes'][0]['value']
        print("probability of survival:",p_survival)

    def test_normal_two_genes_and_drug(self):
        """ Normal request with two genes and a drug
        """

        # empty response
        reasoner_std = { "query_graph": dict()
                       }

        # empty query graph
        reasoner_std["query_graph"] = { "edges": dict(),
                                        "nodes": dict()
                                      }

        # add in evidence gene
        gene1 = ('RAF1', 'ENSEMBL:ENSG00000132155')
        reasoner_std['query_graph']['nodes']['n0'] = { 'category':BIOLINK_GENE,
                                                       'id':'{}'.format(gene1[1])
                                                     }

        gene2 = ('BRCA1', 'ENSEMBL:ENSG00000012048')
        reasoner_std['query_graph']['nodes']['n1'] = { 'category':'biolink:Gene',
                                                       'id':'{}'.format(gene2[1])
                                                     }
        # add in evidence drug
        drug = ('CYCLOPHOSPHAMIDE', 'CHEMBL:CHEMBL88')
        reasoner_std['query_graph']['nodes']['n2'] = { 'category':BIOLINK_DRUG,
                                                       'id':'{}'.format(drug[1])
                                                     }
        # add in disease node
        disease = ('Breast_Cancer', 'MONDO:0007254')
        reasoner_std['query_graph']['nodes']['n3'] = { 'category':BIOLINK_DISEASE,
                                                       'id':'{}'.format(disease[1])
                                                     }
        # add target survival node
        phenotype = ('Survival_Time', 'EFO:0000714')
        reasoner_std['query_graph']['nodes']['n4'] = { 'category': BIOLINK_PHENOTYPIC_FEATURE,
                                                       'id': '{}'.format(phenotype[1]),
                                                     }
        # link genes/drugs to disease
        reasoner_std['query_graph']['edges']['e0'] = { 'predicate':BIOLINK_GENE_TO_DISEASE_PREDICATE,
                                                       'subject': 'n0',
                                                       'object': 'n3'
                                                     }
        reasoner_std['query_graph']['edges']['e1'] = { 'predicate':BIOLINK_GENE_TO_DISEASE_PREDICATE,
                                                       'subject': 'n1',
                                                       'object': 'n3'
                                                     }
        reasoner_std['query_graph']['edges']['e2'] = { 'predicate':BIOLINK_CHEMICAL_TO_DISEASE_OR_PHENOTYPIC_FEATURE_PREDICATE,
                                                       'subject': 'n2',
                                                       'object': 'n3'
                                                     }
        # link disease to target
        reasoner_std['query_graph']['edges']['e3'] = { 'predicate':BIOLINK_DISEASE_TO_PHENOTYPIC_FEATURE_PREDICATE,
                                                       'subject': 'n3',
                                                       'object': 'n4',
                                                       'properties': { 'qualifier':'>=',
                                                                       'days':970
                                                                     }
                                                     }

        # test input is TRAPI compliant
        validate_Message(reasoner_std) # doesn't return True/False for some reason... Will just present exception if not compliant

        handler = TrapiInterface(query=reasoner_std)
        queries = handler.build_chp_queries()
        queries = handler.run_chp_queries()
        reasoner_std_final = handler.construct_trapi_response()

        # test output is TRAPI compliant
        validate_Message(reasoner_std_final['message'])

        KG = reasoner_std_final['message']['knowledge_graph']
        for edge_key in KG['edges'].keys():
            edge = KG['edges'][edge_key]
            if edge['predicate'] == BIOLINK_DISEASE_TO_PHENOTYPIC_FEATURE_PREDICATE:
                p_survival = edge['attributes'][0]['value']
        print("probability of survival:",p_survival)

    def test_no_evidence(self):
        """ Test with no evidence
        """
	# empty response
        reasoner_std = { "query_graph": dict()
                       }

        # empty query graph
        reasoner_std["query_graph"] = { "edges": dict(),
                                        "nodes": dict()
                                      }

        # add in disease node
        disease = ('Breast_Cancer', 'MONDO:0007254')
        reasoner_std['query_graph']['nodes']['n0'] = { 'category':BIOLINK_DISEASE,
                                                       'id':'{}'.format(disease[1])
                                                     }
        # add target survival node
        phenotype = ('Survival_Time', 'EFO:0000714')
        reasoner_std['query_graph']['nodes']['n1'] = { 'category': BIOLINK_PHENOTYPIC_FEATURE,
                                                       'id': '{}'.format(phenotype[1]),
                                                     }
        # link disease to target
        reasoner_std['query_graph']['edges']['e0'] = { 'predicate':BIOLINK_DISEASE_TO_PHENOTYPIC_FEATURE_PREDICATE,
                                                       'subject': 'n0',
                                                       'object': 'n1',
                                                       'properties': { 'qualifier':'>=',
                                                                       'days':970
                                                                     }
                                                     }

        # test input is TRAPI compliant
        validate_Message(reasoner_std) # doesn't return True/False for some reason... Will just present exception if not compliant

        handler = TrapiInterface(query=reasoner_std)
        queries = handler.build_chp_queries()
        queries = handler.run_chp_queries()
        reasoner_std_final = handler.construct_trapi_response()

        # test output is TRAPI compliant
        validate_Message(reasoner_std_final['message'])

        KG = reasoner_std_final['message']['knowledge_graph']
        for edge_key in KG['edges'].keys():
            edge = KG['edges'][edge_key]
            if edge['predicate'] == BIOLINK_DISEASE_TO_PHENOTYPIC_FEATURE_PREDICATE:
                p_survival = edge['attributes'][0]['value']
        print("probability of survival:",p_survival)

    def test_no_evidence_omitting_KG_and_results(self):
        """ Test with no evidence, but omitting KG and Results (should be handled by handler)
        """
        # empty response
        reasoner_std = { "query_graph": dict(),
                       }

        # empty query graph
        reasoner_std["query_graph"] = { "edges": dict(),
                                        "nodes": dict()
                                      }
        # add in disease node
        disease = ('Breast_Cancer', 'MONDO:0007254')
        reasoner_std['query_graph']['nodes']['n0'] = { 'category':BIOLINK_DISEASE,
                                                       'id':'{}'.format(disease[1])
                                                     }
        # add target survival node
        phenotype = ('Survival_Time', 'EFO:0000714')
        reasoner_std['query_graph']['nodes']['n1'] = { 'category': BIOLINK_PHENOTYPIC_FEATURE,
                                                       'id': '{}'.format(phenotype[1]),
                                                     }
        # link disease to target
        reasoner_std['query_graph']['edges']['e0'] = { 'predicate':BIOLINK_DISEASE_TO_PHENOTYPIC_FEATURE_PREDICATE,
                                                       'subject': 'n0',
                                                       'object': 'n1',
                                                       'properties': { 'qualifier':'>=',
                                                                       'days':970
                                                                     }
                                                     }

        # test input is TRAPI compliant
        validate_Message(reasoner_std) # doesn't return True/False for some reason... Will just present exception if not compliant

        handler = TrapiInterface(query=reasoner_std)
        queries = handler.build_chp_queries()
        queries = handler.run_chp_queries()
        reasoner_std_final = handler.construct_trapi_response()

        # test output is TRAPI compliant
        validate_Message(reasoner_std_final['message'])

        KG = reasoner_std_final['message']['knowledge_graph']
        for edge_key in KG['edges'].keys():
            edge = KG['edges'][edge_key]
            if edge['predicate'] == BIOLINK_DISEASE_TO_PHENOTYPIC_FEATURE_PREDICATE:
                p_survival = edge['attributes'][0]['value']
        print("probability of survival:",p_survival)

    def test_one_gene(self):
        """ Test with one gene
        """
        # empty response
        reasoner_std = { "query_graph": dict()
                       }
        # empty query graph
        reasoner_std["query_graph"] = { "edges": dict(),
                                        "nodes": dict()
                                      }

        # add in evidence gene
        gene1 = ('RAF1', 'ENSEMBL:ENSG00000132155')
        reasoner_std['query_graph']['nodes']['n0'] = { 'category':BIOLINK_GENE,
                                                       'id':'{}'.format(gene1[1])
                                                     }
        # add in disease node
        disease = ('Breast_Cancer', 'MONDO:0007254')
        reasoner_std['query_graph']['nodes']['n1'] = { 'category':BIOLINK_DISEASE,
                                                       'id':'{}'.format(disease[1])
                                                     }
        # add target survival node
        phenotype = ('Survival_Time', 'EFO:0000714')
        reasoner_std['query_graph']['nodes']['n2'] = { 'category': BIOLINK_PHENOTYPIC_FEATURE,
                                                       'id': '{}'.format(phenotype[1]),
                                                     }
        # link genes/drugs to disease
        reasoner_std['query_graph']['edges']['e0'] = { 'predicate':BIOLINK_GENE_TO_DISEASE_PREDICATE,
                                                       'subject': 'n0',
                                                       'object': 'n1'
                                                     }
        # link disease to target
        reasoner_std['query_graph']['edges']['e1'] = { 'predicate':BIOLINK_DISEASE_TO_PHENOTYPIC_FEATURE_PREDICATE,
                                                       'subject': 'n1',
                                                       'object': 'n2',
                                                       'properties': { 'qualifier':'>=',
                                                                       'days':970
                                                                     }
                                                     }

        # test input is TRAPI compliant
        validate_Message(reasoner_std) # doesn't return True/False for some reason... Will just present exception if not compliant

        handler = TrapiInterface(query=reasoner_std)
        queries = handler.build_chp_queries()
        queries = handler.run_chp_queries()
        reasoner_std_final = handler.construct_trapi_response()

        # test output is TRAPI compliant
        validate_Message(reasoner_std_final['message'])

        KG = reasoner_std_final['message']['knowledge_graph']
        for edge_key in KG['edges'].keys():
            edge = KG['edges'][edge_key]
            if edge['predicate'] == BIOLINK_DISEASE_TO_PHENOTYPIC_FEATURE_PREDICATE:
                p_survival = edge['attributes'][0]['value']
        print("probability of survival:",p_survival)

    def test_one_drug(self):
        """ Test with one drug
        """
        # empty response
        reasoner_std = { "query_graph": dict()
                       }
        # empty query graph
        reasoner_std["query_graph"] = { "edges": dict(),
                                        "nodes": dict()
                                      }

        # add in evidence drug
        drug = ('CYCLOPHOSPHAMIDE', 'CHEMBL:CHEMBL88')
        reasoner_std['query_graph']['nodes']['n0'] = { 'category':BIOLINK_DRUG,
                                                       'id':'{}'.format(drug[1])
                                                     }
         # add in disease node
        disease = ('Breast_Cancer', 'MONDO:0007254')
        reasoner_std['query_graph']['nodes']['n1'] = { 'category':BIOLINK_DISEASE,
                                                       'id':'{}'.format(disease[1])
                                                     }
        # add target survival node
        phenotype = ('Survival_Time', 'EFO:0000714')
        reasoner_std['query_graph']['nodes']['n2'] = { 'category': BIOLINK_PHENOTYPIC_FEATURE,
                                                       'id': '{}'.format(phenotype[1]),
                                                     }
        # link genes/drugs to disease
        reasoner_std['query_graph']['edges']['e0'] = { 'predicate':BIOLINK_CHEMICAL_TO_DISEASE_OR_PHENOTYPIC_FEATURE_PREDICATE,
                                                       'subject': 'n0',
                                                       'object': 'n1'
                                                     }

        # link disease to target
        reasoner_std['query_graph']['edges']['e1'] = { 'predicate':BIOLINK_DISEASE_TO_PHENOTYPIC_FEATURE_PREDICATE,
                                                       'subject': 'n1',
                                                       'object': 'n2',
                                                       'properties': { 'qualifier':'>=',
                                                                       'days':970
                                                                     }
                                                     }

        # test input is TRAPI compliant
        validate_Message(reasoner_std) # doesn't return True/False for some reason... Will just present exception if not compliant

        handler = TrapiInterface(query=reasoner_std)
        queries = handler.build_chp_queries()
        queries = handler.run_chp_queries()
        reasoner_std_final = handler.construct_trapi_response()

        # test output is TRAPI compliant
        validate_Message(reasoner_std_final['message'])

        KG = reasoner_std_final['message']['knowledge_graph']
        for edge_key in KG['edges'].keys():
            edge = KG['edges'][edge_key]
            if edge['predicate'] == BIOLINK_DISEASE_TO_PHENOTYPIC_FEATURE_PREDICATE:
                p_survival = edge['attributes'][0]['value']
        print("probability of survival:",p_survival)

    def test_no_target(self):
        """ Test with no target (should crash)
        """
        # empty response
        reasoner_std = { "query_graph": dict()
                       }

        # empty query graph
        reasoner_std["query_graph"] = { "edges": dict(),
                                        "nodes": dict()
                                      }
        handler = TrapiInterface(query=reasoner_std)
        queries = handler.build_chp_queries()
        print(handler.error_msg)

    def test_no_disease(self):
        """ Test with no disease
        """
        # empty response
        reasoner_std = { "query_graph": dict()
                       }

        # empty query graph
        reasoner_std["query_graph"] = { "edges": dict(),
                                        "nodes": dict()
                                      }

        # add target survival node
        phenotype = ('Survival_Time', 'EFO:0000714')
        reasoner_std['query_graph']['nodes']['n2'] = { 'category': BIOLINK_PHENOTYPIC_FEATURE,
                                                       'id': '{}'.format(phenotype[1]),
                                                     }
        handler = TrapiInterface(query=reasoner_std)
        queries = handler.build_chp_queries()
        print(handler.error_msg)

    def test_default_survival(self):
        """ Test default survival
        """
        # empty response
        reasoner_std = { "query_graph": dict()
                       }

        # empty query graph
        reasoner_std["query_graph"] = { "edges": dict(),
                                        "nodes": dict()
                                      }

        # add in evidence gene
        gene1 = ('RAF1', 'ENSEMBL:ENSG00000132155')
        reasoner_std['query_graph']['nodes']['n0'] = { 'category':BIOLINK_GENE,
                                                       'id':'{}'.format(gene1[1])
                                                     }
        # add in disease node
        disease = ('Breast_Cancer', 'MONDO:0007254')
        reasoner_std['query_graph']['nodes']['n1'] = { 'category':BIOLINK_DISEASE,
                                                       'id':'{}'.format(disease[1])
                                                     }
        # add target survival node
        phenotype = ('Survival_Time', 'EFO:0000714')
        reasoner_std['query_graph']['nodes']['n2'] = { 'category': BIOLINK_PHENOTYPIC_FEATURE,
                                                       'id': '{}'.format(phenotype[1]),
                                                     }
        # link genes/drugs to disease
        reasoner_std['query_graph']['edges']['e0'] = { 'predicate':BIOLINK_GENE_TO_DISEASE_PREDICATE,
                                                       'subject': 'n0',
                                                       'object': 'n1'
                                                     }
        # link disease to target
        reasoner_std['query_graph']['edges']['e1'] = { 'predicate':BIOLINK_DISEASE_TO_PHENOTYPIC_FEATURE_PREDICATE,
                                                       'subject': 'n1',
                                                       'object': 'n2',
                                                     }

        # test input is TRAPI compliant
        validate_Message(reasoner_std) # doesn't return True/False for some reason... Will just present exception if not compliant

        handler = TrapiInterface(query=reasoner_std)
        queries = handler.build_chp_queries()
        queries = handler.run_chp_queries()
        reasoner_std_final = handler.construct_trapi_response()

        # test output is TRAPI compliant
        validate_Message(reasoner_std_final['message'])

        KG = reasoner_std_final['message']['knowledge_graph']
        for edge_key in KG['edges'].keys():
            edge = KG['edges'][edge_key]
            if edge['predicate'] == BIOLINK_DISEASE_TO_PHENOTYPIC_FEATURE_PREDICATE:
                p_survival = edge['attributes'][0]['value']
        print("probability of survival:",p_survival)

    def test_illegal_gene_to_disease(self):
        """ backwards gene to disease edge. Should return appropriate error message
        """

        # empty response
        reasoner_std = { "query_graph": dict()
                       }
        # empty query graph
        reasoner_std["query_graph"] = { "edges": dict(),
                                        "nodes": dict()
                                      }

        # add in evidence gene
        gene1 = ('RAF1', 'ENSEMBL:ENSG00000132155')
        reasoner_std['query_graph']['nodes']['n0'] = { 'category':BIOLINK_GENE,
                                                       'id':'{}'.format(gene1[1])
                                                     }
        # add in disease node
        disease = ('Breast_Cancer', 'MONDO:0007254')
        reasoner_std['query_graph']['nodes']['n1'] = { 'category':BIOLINK_DISEASE,
                                                       'id':'{}'.format(disease[1])
                                                     }
        # add target survival node
        phenotype = ('Survival_Time', 'EFO:0000714')
        reasoner_std['query_graph']['nodes']['n2'] = { 'category': BIOLINK_PHENOTYPIC_FEATURE,
                                                       'id': '{}'.format(phenotype[1]),
                                                     }
        # link genes/drugs to disease
        reasoner_std['query_graph']['edges']['e0'] = { 'predicate':BIOLINK_GENE_TO_DISEASE_PREDICATE,
                                                       'subject': 'n1',
                                                       'object': 'n0'
                                                     }
        # link disease to target
        reasoner_std['query_graph']['edges']['e1'] = { 'predicate':BIOLINK_DISEASE_TO_PHENOTYPIC_FEATURE_PREDICATE,
                                                       'subject': 'n1',
                                                       'object': 'n2',
                                                       'properties': { 'qualifier':'>=',
                                                                       'days':970
                                                                     }
                                                     }

        # test input is TRAPI compliant
        validate_Message(reasoner_std) # doesn't return True/False for some reason... Will just present exception if not compliant

        handler = TrapiInterface(query=reasoner_std)
        queries = handler.build_chp_queries()
        self.assertEqual('Default type detected. Gene to disease edge (edge id: e0) malformed.', handler.error_msg)

    def test_illegal_drug_to_disease(self):
        """ drug to disease edge is backwards. Should return appropriate error message.
        """
        # empty response
        reasoner_std = { "query_graph": dict()
                       }
        # empty query graph
        reasoner_std["query_graph"] = { "edges": dict(),
                                        "nodes": dict()
                                      }

        # add in evidence drug
        drug = ('CYCLOPHOSPHAMIDE', 'CHEMBL:CHEMBL88')
        reasoner_std['query_graph']['nodes']['n0'] = { 'category':BIOLINK_DRUG,
                                                       'id':'{}'.format(drug[1])
                                                     }
         # add in disease node
        disease = ('Breast_Cancer', 'MONDO:0007254')
        reasoner_std['query_graph']['nodes']['n1'] = { 'category':BIOLINK_DISEASE,
                                                       'id':'{}'.format(disease[1])
                                                     }
        # add target survival node
        phenotype = ('Survival_Time', 'EFO:0000714')
        reasoner_std['query_graph']['nodes']['n2'] = { 'category': BIOLINK_PHENOTYPIC_FEATURE,
                                                       'id': '{}'.format(phenotype[1]),
                                                     }
        # link genes/drugs to disease
        reasoner_std['query_graph']['edges']['e0'] = { 'predicate':BIOLINK_CHEMICAL_TO_DISEASE_OR_PHENOTYPIC_FEATURE_PREDICATE,
                                                       'subject': 'n1',
                                                       'object': 'n0'
                                                     }

        # link disease to target
        reasoner_std['query_graph']['edges']['e1'] = { 'predicate':BIOLINK_DISEASE_TO_PHENOTYPIC_FEATURE_PREDICATE,
                                                       'subject': 'n1',
                                                       'object': 'n2',
                                                       'properties': { 'qualifier':'>=',
                                                                       'days':970
                                                                     }
                                                     }
        # test input is TRAPI compliant
        validate_Message(reasoner_std) # doesn't return True/False for some reason... Will just present exception if not compliant

        handler = TrapiInterface(query=reasoner_std)
        queries = handler.build_chp_queries()
        self.assertEqual('Default type detected. Chemical to disease edge (edge id: e0) malformed.', handler.error_msg)

    def test_illegal_disease_to_phenotype(self):
        """ disease to phenotype edge is backwards. Should return appropriate error message.
        """
        # empty response
        reasoner_std = { "query_graph": dict()
                       }
        # empty query graph
        reasoner_std["query_graph"] = { "edges": dict(),
                                        "nodes": dict()
                                      }

        # add in evidence drug
        drug = ('CYCLOPHOSPHAMIDE', 'CHEMBL:CHEMBL88')
        reasoner_std['query_graph']['nodes']['n0'] = { 'category':BIOLINK_DRUG,
                                                       'id':'{}'.format(drug[1])
                                                     }
         # add in disease node
        disease = ('Breast_Cancer', 'MONDO:0007254')
        reasoner_std['query_graph']['nodes']['n1'] = { 'category':BIOLINK_DISEASE,
                                                       'id':'{}'.format(disease[1])
                                                     }
        # add target survival node
        phenotype = ('Survival_Time', 'EFO:0000714')
        reasoner_std['query_graph']['nodes']['n2'] = { 'category': BIOLINK_PHENOTYPIC_FEATURE,
                                                       'id': '{}'.format(phenotype[1]),
                                                     }
        # link genes/drugs to disease
        reasoner_std['query_graph']['edges']['e0'] = { 'predicate':BIOLINK_CHEMICAL_TO_DISEASE_OR_PHENOTYPIC_FEATURE_PREDICATE,
                                                       'subject': 'n0',
                                                       'object': 'n1'
                                                     }

        # link disease to target
        reasoner_std['query_graph']['edges']['e1'] = { 'predicate':BIOLINK_DISEASE_TO_PHENOTYPIC_FEATURE_PREDICATE,
                                                       'subject': 'n2',
                                                       'object': 'n1',
                                                       'properties': { 'qualifier':'>=',
                                                                       'days':970
                                                                     }
                                                     }
        # test input is TRAPI compliant
        validate_Message(reasoner_std) # doesn't return True/False for some reason... Will just present exception if not compliant

        handler = TrapiInterface(query=reasoner_std)
        queries = handler.build_chp_queries()
        self.assertEqual('Default type detected. Disease to phenotype edge (edge id: e1) malformed.', handler.error_msg)

    def test_illegal_edge(self):
        """ uses gene to drug edge that is incompatible with this query type. Should return appropriate error message.
        """
        # empty response
        reasoner_std = { "query_graph": dict()
                       }
        # empty query graph
        reasoner_std["query_graph"] = { "edges": dict(),
                                        "nodes": dict()
                                      }

        # add in evidence drug
        drug = ('CYCLOPHOSPHAMIDE', 'CHEMBL:CHEMBL88')
        reasoner_std['query_graph']['nodes']['n0'] = { 'category':BIOLINK_DRUG,
                                                       'id':'{}'.format(drug[1])
                                                     }
         # add in disease node
        disease = ('Breast_Cancer', 'MONDO:0007254')
        reasoner_std['query_graph']['nodes']['n1'] = { 'category':BIOLINK_DISEASE,
                                                       'id':'{}'.format(disease[1])
                                                     }
        # add target survival node
        phenotype = ('Survival_Time', 'EFO:0000714')
        reasoner_std['query_graph']['nodes']['n2'] = { 'category': BIOLINK_PHENOTYPIC_FEATURE,
                                                       'id': '{}'.format(phenotype[1]),
                                                     }
        # link genes/drugs to disease
        reasoner_std['query_graph']['edges']['e0'] = { 'predicate':BIOLINK_CHEMICAL_TO_DISEASE_OR_PHENOTYPIC_FEATURE_PREDICATE,
                                                       'subject': 'n0',
                                                       'object': 'n1'
                                                     }

        # link disease to target
        reasoner_std['query_graph']['edges']['e1'] = { 'predicate':BIOLINK_DISEASE_TO_PHENOTYPIC_FEATURE_PREDICATE,
                                                       'subject': 'n1',
                                                       'object': 'n2',
                                                       'properties': { 'qualifier':'>=',
                                                                       'days':970
                                                                     }
                                                     }

        reasoner_std['query_graph']['edges']['e2'] = { 'predicate':BIOLINK_CHEMICAL_TO_GENE_PREDICATE,
                                                       'subject': 'n0',
                                                       'object': 'n1'
                                                     }

        # test input is TRAPI compliant
        validate_Message(reasoner_std) # doesn't return True/False for some reason... Will just present exception if not compliant

        handler = TrapiInterface(query=reasoner_std)
        queries = handler.build_chp_queries()
        self.assertEqual('Default type detected. Received edge between gene and chemical (edge id: e2). This edge is incompatible with default type.', handler.error_msg)

    def test_unknown_edge(self):
        """ tests unknown edge. Should return appropriate error message.
        """
         # empty response
        reasoner_std = { "query_graph": dict()
                       }
        # empty query graph
        reasoner_std["query_graph"] = { "edges": dict(),
                                        "nodes": dict()
                                      }

        # add in evidence drug
        drug = ('CYCLOPHOSPHAMIDE', 'CHEMBL:CHEMBL88')
        reasoner_std['query_graph']['nodes']['n0'] = { 'category':BIOLINK_DRUG,
                                                       'id':'{}'.format(drug[1])
                                                     }
         # add in disease node
        disease = ('Breast_Cancer', 'MONDO:0007254')
        reasoner_std['query_graph']['nodes']['n1'] = { 'category':BIOLINK_DISEASE,
                                                       'id':'{}'.format(disease[1])
                                                     }
        # add target survival node
        phenotype = ('Survival_Time', 'EFO:0000714')
        reasoner_std['query_graph']['nodes']['n2'] = { 'category': BIOLINK_PHENOTYPIC_FEATURE,
                                                       'id': '{}'.format(phenotype[1]),
                                                     }
        # link genes/drugs to disease
        reasoner_std['query_graph']['edges']['e0'] = { 'predicate':BIOLINK_CHEMICAL_TO_DISEASE_OR_PHENOTYPIC_FEATURE_PREDICATE,
                                                       'subject': 'n0',
                                                       'object': 'n1'
                                                     }

        # link disease to target
        reasoner_std['query_graph']['edges']['e1'] = { 'predicate':BIOLINK_DISEASE_TO_PHENOTYPIC_FEATURE_PREDICATE,
                                                       'subject': 'n1',
                                                       'object': 'n2',
                                                       'properties': { 'qualifier':'>=',
                                                                       'days':970
                                                                     }
                                                     }

        reasoner_std['query_graph']['edges']['e2'] = { 'predicate':'biolink:unknown',
                                                       'subject': 'n0',
                                                       'object': 'n1'
                                                     }

        # test input is TRAPI compliant
        validate_Message(reasoner_std) # doesn't return True/False for some reason... Will just present exception if not compliant

        handler = TrapiInterface(query=reasoner_std)
        queries = handler.build_chp_queries()
        self.assertEqual('Default type detected. Unknown predicate type for edge (edge id: e2).', handler.error_msg)

if __name__ == '__main__':
    unittest.main()
