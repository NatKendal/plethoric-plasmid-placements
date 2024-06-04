from pgmpy.factors.discrete import DiscreteFactor

class NoisyOrFactor(DiscreteFactor):
    """
    Special factor generated when preforming operations on BinaryNoisyOrCPDs.

    Arguments:
    operation:
        one of "product" or "marginalize"
    *references:
        if product, then at least 2 references.
        if marginalize, then exactly one reference.

    Things to implement:
    is_valid_cpd
    product, marginalize
    """
    def __init__(self, operation, *references):
        super().__init__(variables=[], cardinality=[], values=[[]], state_names={})
        self.operation = operation
        if operation == "product":
            if len(references) < 2:
                raise ValueError("Expected at least two references to multiply.")
            self.variables = list(set().union(*[ref.variables for ref in references]))
            self.cardinality = None

        elif operation == "marginalize":
            if len(references) != 1:
                raise ValueError("Expected exactly one reference to marginalize.")
        else:
            raise ValueError("Expected operation to be either 'product' or 'marginalize'.")


    def get_value(self, **kwargs):
        pass

    def _get_variables(self):
        pass

    def _get_cardinality(self):
        pass
