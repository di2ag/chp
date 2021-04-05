###########################################
# General Trapi Handler Exceptions
###########################################

class UnidentifiedQueryType(Exception):

    def __str__(self):
        return 'Unidentified query type. Please see https://github.com/di2ag/chp_client for details on our query types.'

class UnidentifiedGeneCurie(Exception):

    def __init__(self, *args):
        self.curie = args[0]

    def __str__(self):
        return 'Unidentified gene curie: {}.'.format(self.curie)

class UnidentifiedDrugCurie(Exception):

    def __init__(self, *args):
        self.curie = args[0]

    def __str__(self):
        return 'Unidentified chemical curie: {}.'.format(self.curie)

class UnidentifiedPhenotypeCurie(Exception):

    def __init__(self, *args):
        self.curie = args[0]

    def __str__(self):
        return 'Unidentified phenotypic curie {}.'.format(self.curie)

class TooManyContributionNodes(Exception):

    def __str__(self):
        return 'Can only have 1 node for contributions.'

class TooManyDiseaseNodes(Exception):

    def __str__(self):
        return 'Can only have 1 node for disease.'

class TooManyPhenotypeNodes(Exception):

    def __str__(self):
        return 'Can only have 1 or no node for disease.'


class MalformedSubjectObjectOnGeneToDisease(Exception):

    def __init__(self, *args):
        self.edge_id = args[0]

    def __str__(self):
        return 'Malformed subject and/or object nodes in gene to disease edge (edge id: {}).'.format(self.edge_id)

class MalformedSubjectObjectOnDrugToDisease(Exception):

    def __init__(self, *args):
        self.edge_id = args[0]

    def __str__(self):
        return 'Malformed subject and/or object nodes in drug to disease edge (edge id: {}).'.format(self.edge_id)


class MalformedSubjectObjectOnDiseaseToPhenotype(Exception):

    def __init__(self, *args):
        self.edge_id = args[0]

    def __str__(self):
        return 'Malformed subject and/or object nodes in disease to phenotype edge (edge id: {}).'.format(self.edge_id)

class MalformedSubjectObjectOnDrugGene(Exception):

    def __init__(self, *args):
        self.edge_id = args[0]

    def __str__(self):
        return 'Malformed subject and/or object nodes in Drug Gene edge (edge id: {}).'.format(self.edge_id)

class UnexpectedEdgeType(Exception):

    def __init__(self, *args):
        self.edge_id = args[0]

    def __str__(self):
        return 'Unexpected edge type (edge id: {}).'.format(self.edge_id)



#########################################
# Wildcard Query Exceptions
#########################################

class IncompatibleWildcardEdge(Exception):

    def __init__(self, *args):
        self.edge_id = args[0]

    def __str__(self):
        return 'Wildcard type detected. Received incompatible edge type for wildcard query (edge id: {}).'.format(self.edge_id)

#########################################
# OneHop Query Exceptions
#########################################

class IncompatibleDrugGeneOneHopEdge(Exception):

    def __init__(self, *args):
        self.edge_id = args[0]

    def __str__(self):
        return 'DrugGene Onehop detected. Received incompatible edge type for DrugGene onehop query (edge id: {}).'.format(self.edge_id)

#########################################
# Default Query Exceptions
#########################################

class IncompatibleDefaultEdge(Exception):

    def __init__(self, *args):
        self.edge_id = args[0]

    def __str__(self):
        return 'Default type detected. Received incompatible edge type for default query (edge id: {}).'.format(self.edge_id)
