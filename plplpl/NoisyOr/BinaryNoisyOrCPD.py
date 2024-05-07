import itertools

import numpy as np

from pgmpy import config
from pgmpy.factors.discrete import TabularCPD
from pgmpy.global_vars import logger

class NoisyOrValues(object):
    def __init__(self, reference):
        self.reference = reference

    def __getitem__(self, key):
        if len(key) != len(self.reference.variables):
            raise KeyError("Expected " + str(len(self.reference.variables)) + " terms, got " + str(len(key)) + ".")
        return self.reference.get_value(**{self.reference.variables[i]:key[i] for i in range(len(key))})
    
    def __len__(self):
        return np.prod(self.reference.cardinality)

    def __array__(self):
        return np.reshape(self.reference.get_values(), self.reference.cardinality)


class BinaryNoisyOrCPD(TabularCPD):
    """
    Noisy OR conditional probability distribution for a binary variable.

    NOTE:
        - Implicitly assumes that variable and evidence variables are binary.
        - Calculates values of queries on the fly.

    TODO:
        - Memoize computed values? Perhaps up to maximumTableSize?

    variable: string that is a valid python variable name
        The variable whose CPD is defined.

    internalProbability: float
        Internal chance of being 1 that is ORed alongside the evidence.

    evidence: array-like of variable names
        List of variables with incoming edges.

    evidence_noise: array-like of floats
        Noise on the incoming edges from evidence.
        Must match length of evidence.
        Noise is given as the probability that a 1 is passed through edge.
        So interference chance is one minus given noise.

    maxTableSize: integer
        Maximum number of values to return in any case.
        If more than MaxTableSize values are requested, throw an error instead of running for a long time.

    state_names: see pgmpy.utils.StateNameMixin

    debug: integer
        Configure warning level.
            >= 2: Enables alert when using state number instead of name.
    """
    def __init__(self, variable, internalProbability, evidence=None, evidence_noise=None, maxTableSize=257, state_names={}, debug=0):
        # Setup initial internal probability.
        self.internalProbability = internalProbability

        # Catch weirdness with statenames having more variables than we pass to superclass.
        if variable in state_names:
            restricted_state_names = {variable: state_names[variable]}
        else:
            restricted_state_names = {}

        # Initialize superclass.
        super().__init__(variable, 2, [[1.0-internalProbability],[internalProbability]], evidence=None, evidence_card=None, state_names=restricted_state_names)
        # At this point we have:
        # self.variable = variable
        # self.variable_card = 2
        # self.variables = [variable]
        # self.cardinality = [2]
        # self.values = [[1.0-internalProbability],[internalProbability]]

        if len(evidence) != len(evidence_noise):
            raise ValueError("Length of evidence_noise doesn't match length of evidence")

        # Setup evidence.
        self.evidence = evidence
        self.variables.extend(evidence)
        self.cardinality = np.array([2]*(len(evidence) + 1))
        self.evidence_noise = evidence_noise

        # Setup maximum table size.
        self.maxTableSize = maxTableSize

        # Setup debug.
        self.debug = debug

        # Set values to be a special interface back into this object.
        self.values = NoisyOrValues(self)

    # Helper function to interpret a state name back to a state number.
    def get_state(self, variable, state):
        try:
            s = self.name_to_no[variable][state]
        except:
            if self.debug >= 2:
                logger.info(f"Using {variable} state as number instead of name.")
            s = state
        return s

    
    #
    # Calculation Functions:
    #

    # Calculate the probability of a complete set of state assignments.
    def get_value(self, **kwargs):
        # Check if all variables are assigned.
        if self.variable not in kwargs:
            raise ValueError("List of states is missing variable: " + str(self.variable))
        for variable in self.evidence:
            if variable not in kwargs:
                raise ValueError("List of states is missing variable: " + str(variable))

        # Initial failure probability.
        failureChance = 1.0 - self.internalProbability

        # Compute probability that each other check fails.
        for i in range(len(self.evidence)):
            if self.get_state(self.evidence[i], kwargs[self.evidence[i]]) == 1:
                failureChance = failureChance * (1.0 - self.evidence_noise[i])

        if self.get_state(self.variable, kwargs[self.variable]) == 1:
            return 1.0 - failureChance
        else:
            return failureChance

    # Calculate the probability of this variable given state assignments to evidence.
    def get_self_values(self, **kwargs):
        if self.variable in kwargs:
            raise ValueError("Can't get probability of each assignment to " + str(self.variable) + " if it's assigned.")
        kwargs[self.variable] = 0
        zeroProb = self.get_value(**kwargs)
        return [[zeroProb], [1-zeroProb]]

    # Calculates get_value for all possible assignments to variables.
    # Will fail if it would calculate more than self.maxTableSize values.
    # Should print in appropriate pgmpy format. #TODO check this
    def get_values(self):
        if np.prod(self.cardinality) > self.maxTableSize:
            raise Exception("Table too large to compute.")
        values = []
        for assignment in itertools.product(*[range(card) for card in self.cardinality]):
            values.append(self.get_value(**{self.variables[i]:assignment[i] for i in range(len(assignment))}))
        return np.reshape(np.array(values), (self.cardinality[0], np.prod(self.cardinality[1:])))

    
    #
    # Overrides:
    #

    # You can't set the value of a specific instance.
    def set_value(self, value, **kwargs):
        raise Exception("Not defined for NoisyOrCPD.")

    # A BinaryNoisyOrCPD is always a valid CPD.
    def is_valid_cpd(self):
        return True

    def to_factor(self):
        factor = DiscreteFactor.__new__(DiscreteFactor)
        factor.variables = self.variables.copy()
        factor.cardinality = self.cardinality.copy()
        factor.values = np.array(self.values).copy()
        factor.state_names = self.state_names.copy()
        factor.name_to_no = self.name_to_no.copy()
        factor.no_to_name = self.no_to_name.copy()
        return factor
    
    #
    # CPD Functions:
    #

    # String or array like to remove from evidence.
    def delete_evidence(self, evidenceVariable):
        if isinstance(evidenceVariable, (str, int)) and evidenceVariable in self.evidence:
            i = self.evidence.index(evidenceVariable)
            del self.evidence[i]
            del self.evidence_noise[i]
            del self.variables[i+1]
            self.cardinality = np.delete(self.cardinality, i+1)
        else:
            for e in evidenceVariable:
                i = self.evidence.index(e)
                del self.evidence[i]
                del self.evidence_noise[i]
                del self.variables[i+1]
                self.cardinality = np.delete(self.cardinality, i+1)
