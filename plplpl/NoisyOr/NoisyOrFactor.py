from pgmpy.factors.discrete import DiscreteFactor

class NoisyOrFactor(DiscreteFactor):
    """
    Special factor generated when preforming operations on BinaryNoisyOrCPDs.

    considerations:
    is_valid_cpd
    product, marginalize
    """
    def __init__(self):
        super().__init__(variables=[], cardinality=[], values=[[]], state_names={})

    
