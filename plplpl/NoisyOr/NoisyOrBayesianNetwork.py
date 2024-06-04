from pgmpy.models import BayesianNetwork
from pgmpy.factors.discrete import TabularCPD
from pgmpy.factors.continuous import ContinuousFactor

class NoisyOrBayesianNetwork(BayesianNetwork):
    def __init__(self, ebunch=None, latents=set()):
        super().__init__(ebunch=ebunch, latents=latents)
        self.constants = dict() # An additional place to save relevant information to the network.
        self.cpds_dict = dict() # Overwrite CPD list for speed.
        self._validCPDTypes = (TabularCPD, ContinuousFactor)

    # Overwrite property from superclass to just convert internal dictionary to list.
    @property
    def cpds(self):
        return list(self.cpds_dict.values())

    #
    # Overrides:
    #

    # Adapted from pgmpy.BayesianNetwork.add_cpds. Add CPDs to dictionary instead of list.
    def add_cpds(self, *cpds):
        for cpd in cpds:
            if not isinstance(cpd, self._validCPDTypes):
                    raise ValueError("Only " + ", ".join(map(str, self._validCPDTypes)) + " are valid CPD types.")

            if set(cpd.scope()) - set(cpd.scope()).intersection(set(self.nodes())):
                raise ValueError("CPD defined on variable not in the model", cpd)

            if cpd.variable in self.cpds_dict:
                logger.warning(f"Replacing existing CPD for {cpd.variable}")
            self.cpds_dict[cpd.variable] = cpd

    # Adapted from pgmpy.BayesianNetwork.remove_cpds. Remove CPDs from dictionary instead of list.
    def remove_cpds(self, *cpds):
        for cpd in cpds:
            if isinstance(cpd, (str, int)):
                if cpd not in self.cpds_dict:
                    raise ValueError("CPD " + str(cpd) + " not in model.")
                self.cpds_dict.pop("cpd")
            elif isinstance(cpd, self._validCPDTypes):
                if cpd.variable in self.cpds_dict:
                    if cpd == self.cpds_dict[cpd.variable]:
                        self.cpds_dict.pop(cpd.variable)
                    else:
                        raise ValueError("Given CPD doesn't match saved CPD for the same variable.")
                else:
                    raise ValueError("No CPD was found for variable " + str(cpd.variable))
            else:
                raise ValueError("Requires str, int, or " + ", ".join(map(str, self._validCPDTypes)) + ", not " +str(type(cpd)))



    # Adapted from pgmpy.BayesianNetwork.get_cpds. Read CPDs from the dictionary instead of list.
    def get_cpds(self, node=None):
        if node == None:
            return list(self.cpds_dict.values())
        else:
            if node in self.cpds_dict:
                return self.cpds_dict[node]
            else:
                raise ValueError("Node not present in the Directed Graph")

    # Adapted from pgmpy.BayesianNetwork.check_model, with the state name requirement removed.
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
