import itertools

#from pgmpy.factors.discrete import DiscreteFactor
from plplpl.NoisyOr import BinaryNoisyOrCPD

class NoisyOrFactor(object):
    """
    Low overhead factor wrapper

    Arguments:
    operation:
        one of "product", "marginalize", or "reduce"
    references:
        if product, then at least one reference
        if marginalize, or reduce then exactly one reference
    argument:
        if product, ignored
        if marginalize, list of variables to marginalize
        if reduce, list of tuples (variable, state)

    """
    def __init__(self, operation, references, argument=[]):
        #super().__init__(variables=[], cardinality=[], values=[[]], state_names={})
        self.operation = operation
        self.references = []
        self.savedValues = dict()
        self.argument = argument
        self.references = references

        if len(self.references) == 0:
            raise ValueError("Factor is empty.")

        if operation == "product":
            self.variables = sorted(list(set().union(*[ref.variables for ref in references])))
        elif operation == "marginalize":
            if len(references) != 1:
                raise ValueError("Can only marginalize one reference.")
            self.variables = sorted(list(set(self.references[0].variables).difference(set(self.argument))))
        elif operation == "reduce":
            if len(references) != 1:
                raise ValueError("Can only reduce one reference.")
            self.variables = sorted(list(set(self.references[0].variables).difference(set([x[0] for x in self.argument]))))
        else:
            raise ValueError("Expected operation to be either 'product' or 'marginalize'.")

        self.cardinality = [2]*len(self.variables)

    def get_value(self, **kwargs):
        assignment = tuple([kwargs[var] for var in self.variables])
        if assignment in self.savedValues:
            return self.savedValues[assignment]

        newKwargs = kwargs.copy()

        if self.operation == "product":
            result = 1.0
            for reference in self.references:
                result = result * reference.get_value(**newKwargs)
            self.savedValues[assignment] = result
            return result
        elif self.operation == "marginalize":
            marginalizeAssignments = itertools.product([0,1], repeat=len(self.argument))
            total = 0
            for marginalizeAssignment in marginalizeAssignments:
                newKwargs.update({self.argument[i]:marginalizeAssignment[i] for i in range(len(self.argument))})
                total = total + self.references[0].get_value(**newKwargs.copy())
            self.savedValues[assignment] = total
            return total
        elif self.operation == "reduce":
            newKwargs.update({x[0]:x[1] for x in self.argument})
            result = self.references[0].get_value(**newKwargs)
            self.savedValues[assignment] = result
            return result
