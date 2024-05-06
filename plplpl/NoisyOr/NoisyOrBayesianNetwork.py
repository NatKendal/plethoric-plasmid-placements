from pgmpy.models import BayesianNetwork
from pgmpy.factors.discrete import TabularCPD

class NoisyOrBayesianNetwork(BayesianNetwork):
    def __init__(self, ebunch=None, latents=set()):
        super().__init__(ebunch=ebunch, latents=latents)
        self.constants = dict() # An additional place to save relevant information to the network.

    #
    # Overrides:
    #

    # Adapted from BayesianNetwork.check_model, with the state name requirement removed.
    def check_model(self, progressBar=False):
        if progressBar:
            import tqdm
            iterator = tqdm.tqdm(self.nodes())
        else:
            iterator = self.nodes()
        if progressBar:
            print("Starting CPD completeness check.")
        for node in iterator:
            cpd = self.get_cpds(node=node)

            if cpd is None:
                raise ValueError(f"No CPD associated with {node}")
            elif isinstance(cpd, TabularCPD):
                evidence = cpd.get_evidence()
                parents = self.get_parents(node)

                if set(evidence) != set(parents):
                    raise ValueError(f"CPD associated with {node} doesn't have proper parents associated with it.")

                if not cpd.is_valid_cpd():
                    raise ValueError(f"Sum or integral of conditional probabilities for node {node} is not equal to 1.")

        if progressBar:
            print("Starting cardinality match check.")
        for node in iterator:
            cpd = self.get_cpds(node=node)
            for index, node in enumerate(cpd.variables[1:]):
                parent_cpd = self.get_cpds(node)
                if parent_cpd.cardinality[0] != cpd.cardinality[1+index]:
                    raise ValueError(f"The cardinality of {node} doesn't match in it's child nodes.")
