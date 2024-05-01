from pgmpy.models import BayesianNetwork

class NoisyOrBayesianNetwork(BayesianNetwork):
    def __init__(self, ebunch=None, latents=set()):
        super().__init__(ebunch=ebunch, latents=latents)
        self.constants = dict()
