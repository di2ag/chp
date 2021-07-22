#####################################
# General mixin errors
#####################################



#####################################
# One hop mixin errors
#####################################

class TooManyOneHopProxies(Exception):

    def __str__(self):
        return 'One hop queries that use a proxy for implicit two hop associations may only have 1 proxy type.'

class StandardOneHopNoContext(Exception):

    def __str__(self):
        return 'One hop standard probabilistic queries can not use proxies that contain no context.'
