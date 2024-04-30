from pgmpy.factors.discrete import DiscreteFactor

class NoisyOrFactor(DiscreteFactor):
    """
    Special factor generated when preforming operations on BinaryNoisyOrCPDs.
    """
    def __init__(self):
        super().__init__(variables=[], cardinality=[], values=[[]], state_names={})


