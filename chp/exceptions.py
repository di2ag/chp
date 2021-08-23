###########################################
# General Trapi Handler Exceptions
###########################################

class NoSupportedQueriesFound(Exception):

    def __str__(self):
        return 'No CHP core supported queries where found in passed query.'

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

class TooManyContributionNodes(Exception):

    def __str__(self):
        return 'Can only have 1 node for contributions.'
