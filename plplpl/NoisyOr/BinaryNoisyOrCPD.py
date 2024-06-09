import itertools

import numpy as np

from pgmpy import config
from pgmpy.factors.discrete import TabularCPD
from pgmpy.factors.discrete import DiscreteFactor
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
        return np.reshape(self.reference.get_values(), tuple(self.reference.cardinality))


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
    def __init__(self, variable, internalProbability, evidence=[], evidence_noise=[], maxTableSize=257, state_names={}, debug=0):
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
        #self.savedValues = dict()

        # Setup maximum table size.
        self.maxTableSize = maxTableSize

        # Setup debug.
        self.debug = debug

        # Set values to be a special interface back into this object.
        self.values = NoisyOrValues(self)

        # Properly setup state names.
        self.store_state_names(self.variables, self.cardinality, state_names)

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

        #assignment = tuple([kwargs[var] for var in self.variables])
        #if assignment in self.savedValues:
        #    return savedValues[assignment]

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
        return np.reshape(np.array(values), tuple([self.cardinality[0], int(np.prod(self.cardinality[1:]))]))

    # Calculate the probabilities of each state of cpd.variable given a probability of being 1 for each evidence variable. 
    def get_self_values_with_uncertain_evidence(self, **kwargs):
        # Check if all evidence variables are assigned.
        if self.variable in kwargs:
            raise ValueError("Can't get probability of each assignment to " + str(self.variable) + " if it's assigned.")
        for variable in self.evidence:
            if variable not in kwargs:
                raise ValueError("List of states is missing variable: " + str(variable))

        # Initial failure probability.
        failureChance = 1.0 - self.internalProbability

        # Compute weighted probability that each other check fails.
        for i in range(len(self.evidence)):
            failureChance = failureChance * (1.0 - self.evidence_noise[i] * kwargs[self.evidence[i]])

        return [[failureChance], [1-failureChance]]

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
        factor.cardinality = np.array(self.cardinality).copy()
        factor.values = np.array(self.values).copy()
        factor.state_names = self.state_names.copy()
        factor.name_to_no = self.name_to_no.copy()
        factor.no_to_name = self.no_to_name.copy()
        return factor

    def copy(self):
        return BinaryNoisyOrCPD(self.variable, self.internalProbability, evidence=self.evidence.copy(), evidence_noise=self.evidence_noise.copy(), maxTableSize=self.maxTableSize, state_names=self.state_names.copy(), debug=self.debug)

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

    # Reduce/restrict to a subset of the original variables.
    # Expects values as an arraylike of (variable_name, variable_state).
    def reduce(self, values, inplace=True, show_warnings=True):
        # Get the assignment from the values.
        assignment = dict()
        toDelete = set()
        for variable, state in values:
            if variable == self.variable:
                raise ValueError("Reduce not allowed on variable on which CPD is defined")
            assignment[variable] = state
            toDelete.add(variable)
        # Build the assignment.
        for variable in self.variables:
            if variable in assignment:
                continue
            elif variable == self.variable:
                assignment[variable] = 1
            else:
                assignment[variable] = 0
        # Calculate the assignment and replace the internal probability.
        phi = self if inplace else self.copy()
        phi.internalProbability = phi.get_value(**assignment)
        phi.variables = [var for var in phi.variables if var not in toDelete]
        newEvidence = []
        newEvidenceNoise = []
        for i in range(len(phi.evidence)):
            if phi.evidence[i] not in toDelete:
                newEvidence.append(phi.evidence[i])
                newEvidenceNoise.append(phi.evidence_noise[i])
        phi.evidence = newEvidence
        phi.evidence_noise = newEvidenceNoise
        phi.cardinality = [2]*len(phi.variables)
        return phi

    # Reduce/restrict to a subset of the original values.
    # Expects values as an arraylike of (variable_name, probability_of_being_one).
    def reduce_with_uncertain_evidence(self, values, inplace=True, show_warnings=True):
        # Get the assignment from the values.
        assignment = dict()
        toDelete = set()
        for variable, probability in values:
            if variable == self.variable:
                raise ValueError("Reduce not allowed on variable on which CPD is defined")
            assignment[variable] = probability
            toDelete.add(variable)
        # Build the assignment.
        for variable in self.variables:
            if variable in assignment:
                continue
            elif variable == self.variable:
                continue
            else:
                assignment[variable] = 0
        # Calculate the assignment and replace the internal probability.
        phi = self if inplace else self.copy()
        phi.internalProbability = phi.get_self_values_with_uncertain_evidence(**assignment)[1][0]
        phi.variables = [var for var in phi.variables if var not in toDelete]
        newEvidence = []
        newEvidenceNoise = []
        for i in range(len(phi.evidence)):
            if phi.evidence[i] not in toDelete:
                newEvidence.append(phi.evidence[i])
                newEvidenceNoise.append(phi.evidence_noise[i])
        phi.evidence = newEvidence
        phi.evidence_noise = newEvidenceNoise
        phi.cardinality = [2]*len(phi.variables)
        return phi
