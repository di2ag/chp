#####################################
# General mixin errors
#####################################



#####################################
# One hop mixin errors
#####################################

class TooManyOneHopProxies(Exception):

    def __str__(self):
        return 'One hop queries that use a proxy for implicit two hop associations may only have 1 proxy type'
