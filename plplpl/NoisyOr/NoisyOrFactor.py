from collections import ChainMap
import itertools
from weakref import WeakValueDictionary

#from pgmpy.factors.discrete import DiscreteFactor
from plplpl.NoisyOr import BinaryNoisyOrCPD

class NoisyOrFactor(object):
    """
    Low overhead factor wrapper

    Arguments:
    operation:
        one of "constant", "sum", "product", "marginalize", "reduce", or "marginalizeGivenProbability"
    references:
        if constant, then no references
        if product or sum, then at least one reference
        if marginalize or reduce, then exactly one reference
    argument:
        if constant, list(value) (list with one element)
        if product, ignored
        if sum, ignored
        if marginalize, list of variables to marginalize
        if reduce, list of tuples (variable, state)
        if marginalizeGivenProbability, list of tuples (variable, probabilityOfBeingOne)

    """
    # Experimental pointer saving
    # Weak value should play better with multicore and also sequential runs.
    _registry = WeakValueDictionary()
    #_registry = dict()
    def __new__(cls, operation, references, argument=[]):
        factorID = (operation, tuple(sorted([id(ref) for ref in references])), tuple(sorted(argument)))
        if factorID in cls._registry:
            return cls._registry[factorID]
        instance = super().__new__(cls)
        cls._registry[factorID] = instance
        #if len(cls._registry) % 10000 == 0:
        #    print(len(cls._registry))
        return instance

    def __init__(self, operation, references, argument=[]):
        #super().__init__(variables=[], cardinality=[], values=[[]], state_names={})
        self.operation = operation
        self.references = []
        self.savedValues = dict()
        self.argument = argument.copy()
        self.factorPointers = dict()
        self.references = []

        for reference in references:
            self.references.append(reference)

        if (len(self.references) == 0) and (self.operation != "constant") and (self.operation != "simple"):
            raise ValueError("Factor is empty.")
        
        if operation == "constant":
            if len(self.references) > 0:
                raise ValueError("Constant factor had a reference, it shouldn't.")
            self.variables = []
        elif operation == "simple":
            if len(self.references) > 0:
                raise ValueError("Simple factor had a reference, it shouldn't.")
            if len(self.argument) != 3:
                raise ValueError("Expected argument to be of the form [variable, zeroProb, oneProb]")
            else:
                self.variables = [argument[0]]
        elif operation == "product":
            self.variables = sorted(list(set().union(*[ref.variables for ref in references])))
        elif operation == "sum":
            self.variables = sorted(list(set().union(*[ref.variables for ref in references])))
        elif operation == "marginalize":
            if len(references) != 1:
                raise ValueError("Can only marginalize one reference.")
            if len(self.argument) == 0:
                raise ValueError("Can't marginalize nothing.")
            self.variables = sorted(list(set(self.references[0].variables).difference(set(self.argument))))
        elif operation == "reduce":
            if len(references) != 1:
                raise ValueError("Can only reduce one reference.")
            self.variables = sorted(list(set(self.references[0].variables).difference(set([x[0] for x in self.argument]))))
        elif operation == "marginalizeGivenProbability":
            if len(references) != 1:
                raise ValueError("Can only marginalize* one reference.")
            if len(self.argument) == 0:
                raise ValueError("Can't marginalize nothing.")
            self.variables = sorted(list(set(self.references[0].variables).difference(set([x[0] for x in self.argument]))))
        else:
            raise ValueError("Expected operation from [constant, product, marginalize, reduce, marginalizeGivenProbability]")

    def __repr__(self):
        return "<NoisyOrFactor " + str(id(self)) + " over (" + ", ".join(self.variables) + ") doing " + self.operation + " with argument " + repr(self.argument) + " on " + str(len(self.references)) + " reference factors>"

    def marginalizeAll(self):
        if len(self.variables) >= 1:
            marginalizeAssignments = itertools.product([0,1], repeat=len(self.variables))
            total = 0
            v_assignment = {self.variables[i]:0 for i in range(len(self.variables))}
            for marginalizeAssignment in marginalizeAssignments:
                for i in range(len(self.variables)):
                    v_assignment[self.variables[i]] = marginalizeAssignment[i]
                    total = total + self.get_value(**v_assignment)
            return total
        else:
            return self.get_value()

    def get_value(self, **kwargs):
        if self.operation == "constant":
            return self.argument[0]
        elif self.operation == "simple":
            if self.variables[0] not in kwargs:
                raise ValueError("Missing variable " + self.variables[0])
            if kwargs[self.variables[0]] == 0:
                return self.argument[1]
            else:
                return self.argument[2]
        
        assignment = tuple([kwargs[var] for var in self.variables])
        if assignment in self.savedValues:
            return self.savedValues[assignment]
        
        if self.operation == "product":
            result = 1.0
            for reference in self.references:
                result = result * reference.get_value(**kwargs)
            self.savedValues[assignment] = result
            return result
        elif self.operation == "sum":
            result = 0.0
            for reference in self.references:
                result = result + reference.get_value(**kwargs)
            self.savedValues[assignment] = result
            return result
        elif self.operation == "marginalize":
            if len(self.argument) > 0:
                #if len(self.argument) > 6:
                #    print("About to marginalize " + str(len(self.argument)))
                marginalizeAssignments = itertools.product([0,1], repeat=len(self.argument))
                total = 0
                v_assignment = {self.argument[i]:0 for i in range(len(self.argument))}
                for marginalizeAssignment in marginalizeAssignments:
                    for i in range(len(self.argument)):
                        v_assignment[self.argument[i]] = marginalizeAssignment[i]
                    total = total + self.references[0].get_value(**ChainMap(v_assignment, kwargs))
                    #total = total + self.references[0].get_value(**v_assignment, **kwargs)
                self.savedValues[assignment] = total
                #if len(self.argument) > 6:
                #    print("Finished marginalizing " + str(len(self.argument)))
                return total
            else:
                return self.get_value(**kwargs)
        elif self.operation == "reduce":
            v_assignment = {x[0]:x[1] for x in self.argument}
            result = self.references[0].get_value(**ChainMap(v_assignment, kwargs))
            #result = self.references[0].get_value(**v_assignment, **kwargs)
            self.savedValues[assignment] = result
            return result
        elif self.operation == "marginalizeGivenProbability":
            if len(self.argument) > 0:
                marginalizeAssignments = itertools.product([0,1], repeat=len(self.argument))
                total = 0
                v_assignment = {self.argument[i][0]:0 for i in range(len(self.argument))}
                for marginalizeAssignment in marginalizeAssignments:
                    assignmentProb = 1.0
                    for i in range(len(marginalizeAssignment)):
                        v_assignment[self.argument[i][0]] = marginalizeAssignment[i]
                        assignmentProb = assignmentProb * (self.argument[i][1] if (marginalizeAssignment[i] == 1) else (1-self.argument[i][1]))
                    total = total + (assignmentProb * self.references[0].get_value(**ChainMap(v_assignment, kwargs)))
                    #total = total + (assignmentProb * self.references[0].get_value(**v_assignment, **kwargs))
                self.savedValues[assignment] = total
                return total
            else:
                return self.get_value(**kwargs)
