import itertools
import numpy as np

from pgmpy import config
from pgmpy.factors.discrete.CPD import TabularCPD
from pgmpy.global_vars import logger

#TODO(1) is high priority.
#TODO(2) is also relevant.

# Memoize some probabilities?

class GeneralizedNoisyOrCPD(TabularCPD):
    def __init__(self, variable, varible_card, values, nullstate, conditions, evidence=None, evidence_card=None, evidence_noise=None, maxTableSize=257, state_names={}):
        """
        Initialize a GeneralizedNoisyOrCPD

        Parameters
        ----------
        variable: int, string (any hashable python object)
            The variable whose CPD is defined.

        variable_card: integer
            Cardinality/no. of states of `variable`

        values: 2D array, 2D list, or 2D tuple
            NOTE: This is the base likelihood for this cell to be in a given state.
                If there are no evidence variables, then this is the full CPD for this variable.
                Requires 2D array so it can be passed to superclasses, see TabularCPD for specific requirements.
                No evidence variables are passed to superclasses.

        evidence: array-like
            List of variables in evidences(if any) w.r.t. which CPD is defined.
            NOTE: No evidence variables are passed to superclasses.

        evidence_card: array-like of size len(evidence)
            cardinality/no. of states of variables in `evidence`(if any)
            NOTE: No evidence variables are passed to superclasses.

        evidence_noise: list of 2D array-like
            List of 2D arrays which give the noise on links given values for the evidence.
            i.e. each 2D array A has that
            A[i,j] = probability of reading state j given true state i.
            An identity matrix here has a noiseless link.
            All entries here should be in [0,1].

        nullstate: integer
            Value that this variable takes in the absence of any additional evidence.

        conditions: 2D list of tuples (evidence_variable, state)
            conditions[k] = [(x_1, i_1), (x_2, i_2), ...]
            denotes that this variable is in state k if any variable x_j is in state i_j. (After noise.)
            conditions[nullstate] is generally ignored, since it's the remaining probability after all checks.
            We still expect it to be there to offset the indexing.
            Repeats of the same evidence variable are allowed, even within a given conditions[k].
            No repeat of (evidence_variable, state) is permitted across all conditions.
            (variable, state) is permitted.
            i.e. the internal value of this state can be included as an additional modifier.

        maxTableSize: integer
            Maximum number of values to return in any case.
            If more than maxTableSize values are requested in any function, throw an error instead of running long.

        """
        # First initialize superclasses.
        super().__init__(variable, variable_card, values, evidence=None, evidence_card=None, state_names=state_names)
        # At this point we have:
        # self.variable=variable, self.variable_card=variable_card, self.variables=[variable]
        # self.cardinality=[variable_card]
        # self.values=values (Though only relevant if no evidence.)

        # Setup evidence and evidence noise.
        if evidence is not None:
            self.variables.extend(evidence)
            if not len(evidence_card) == len(evidence):
                raise ValueError("Length of evidence_card doesn't match length of evidence.")
            if not len(evidence_noise) == len(evidence):
                raise ValueError("iLength of evidence_noise doesn't match length of evidence.")
            self.cardinality.extend(evidence_card)
            self.evidence_noise = dict()
            for i in range(len(evidence)):
                # TODO: Check that evidence_noise is valid and/or normalize it here.
                self.evidence_noise[evidence[i]] = np.array(evidence_noise[i], dtype=config.get_dtype())

        # Setup nullstate.
        try:
            state = self.name_to_no[self.variable][nullstate]
        except KeyError:
            state = nullstate
        if state < 0 or state >= self.cardinality[0]:
            raise ValueError("nullstate is not a valid state for this variable.")
        self.nullstate = state

        # Verify that conditions is valid and setup conditions.
        self.conditions = dict()
        if len(conditions) != self.cardinality[0]:
            raise ValueError("Length of conditions doesn't match variable cardinality.")
        allConditions = set()
        for i in range(len(conditions)):
            stateConditions = dict()
            for condition in conditions[i]:
                if not isinstance(condition, tuple) or len(condition) != 2:
                    raise ValueError("Conditions must be tuples (evidence_variable, state).")
                try:
                    state = self.name_to_no[condition[0]][condition[1]]
                except KeyError:
                    #logger.info(f"Using {condition[0]} state as number instead of name.")
                    state = s
                if (condition[0], state) in allConditions:
                    raise ValueError("Conditions can't have (evidence_variable, state) repeats.")
                if condition[0] not in self.variables:
                    raise ValueError("Condition references variable that isn't this variable or in evidence.")
                if state < 0 or state >= self.cardinality[self.variables.index(condition[0])]:
                    raise ValueError("Condition references state that isn't valid for variable " + str(variable) + ".")
                allConditions.add((condition[0], state))
                if condition[0] in stateConditions:
                    stateConditions[condition[0]].append(state)
                else:
                    stateConditions[condition[0]] = [state]
            self.conditions[i] = stateConditions

        # Setup tablesize.
        self.maxTableSize = 257

    #
    # String Functions:
    #

    # TODO: Implement properly.
    def __str__(self):
        return object.__repr__(self)
    # TODO: Implement properly.
    def __repr__(self):
        return object.__repr__(self)
    # TODO: Check if this overwriting is necessary.
    def _str(self, phi_or_p="phi", tablefmt="grid", print_state_names=True):
        return self.__str__(self)
    # TODO: Check if this overwriting is necessary.
    def _make_table_str(self, tablefmt="fancy_grid", print_state_names=True, return_list=False):
        return self.__str__(self)
    # TODO: Check if this overwriting is necessary.
    def _truncate_strtable(self, cdf_str):
        return self.__str__(self)
    # TODO: Check if this overwriting is necessary:
    def to_csv(self, filename):
        pass

    #
    # CPD Functions:
    #

    # Calculates the value for a given set of state assignments.
    # Expects sufficiently many `variable=assignment` arguments.
    # Requires an assignment for the variable of this CPD.
    # Also requires assignment of all variables that are relevant to the conditions of the assignment of x.
    def get_value(self, **kwargs):
        if self.variable not in kwargs:
            raise ValueError("List of states is missing variable: " + str(self.variable))
        try:
            state = self.name_to_no[self.variable][kwargs[self.variable]]
        except KeyError:
            #logger.info(f"Using {self.variable} state as number instead of name.")
            state = kwargs[self.variable]
        if state == self.nullstate:
            failureChance = 1.0
            for i in range(self.cardinality[0]):
                if i == self.nullstate:
                    continue
                for var, condition in self.conditions[i]:
                    if var not in kwargs:
                        raise ValueError("List of states is missing variable: " + str(var))
                    try:
                        evidence_state = self.name_to_no[var][kwargs[var]]
                    except KeyError:
                        #logger.info(f"Using {var} state as number instead of name.")
                        evidence_state = kwargs[var]
                    failureChance = failureChance * (1.0 - sum(self.evidence_noise[var][evidence_state][read_state] for read_state in condition))
            return failureChance
        else:
            failureChance = 1.0
            for var, condition in self.conditions[state]:
                if var not in kwargs:
                    raise ValueError("List of states is missing variable: " + str(var))
                try:
                    evidence_state = self.name_to_no[var][kwargs[var]]
                except KeyError:
                    #logger.info(f"Using {var} state as number instead of name.")
                    evidence_state = kwargs[var]
                failureChance = failureChance * (1.0 - sum(self.evidence_noise[var][evidence_state][read_state] for read_state in condition))
            return 1.0 - failureChance

    # Calculates get_value for all possible assignments to variables.
    # Will fail if it would calculate more than self.maxTableSize values.
    # Should print in appropriate pgmpy format. #TODO check this
    def get_values(self):
        if np.prod(self.cardinalities) > self.maxTableSize:
            raise Exception("Table too large to compute.")
        values = []
        for assignment in itertools.product(*[range(card) for card in cardinalities]):
            values.append(self.get_value(**{self.variables[i]:assignment[i] for i in range(len(assignment))}))
        return values.reshape((self.cardinality[0], np.prod(self.cardinality[1:])))

    # Raises error, the Noisy Or model has different ways of changing CPD.
    def set_value(self, value, **kwargs):
        raise Exception("Not defined for GeneralizedNoisyOrCPD.")

    # TODO(2): Have to think about what this means in context.
    def identity_factor(self):
        pass

    # TODO(1): Implement.
    def marginalize(self, variables, inplace=True):
        pass

    # TODO(1): Implement.
    def maximize(self, variables, inplace=True):
        pass

    # TODO(1): Implement. Might not be relevant given how we do this?
    def normalize(self, inplace=True):
        pass

    # TODO(1): Implement.
    def reduce(self, values, inplace=True, show_warnings=True):
        pass

    # TODO(1): Implement.
    def sum(self, ph1, inplace=True):
        pass

    # TODO(1): Implement.
    def product(self, phi1, inplace=True):
        pass

    # TODO(1): Implement.
    def divide(self, phi1, inplace=True):
        pass

    # TODO(2): Implement.
    def sample(self, n):
        pass

    # TODO(2): Implement.
    def copy(self):
        pass

    # TODO(2): Implement.
    def is_valid_cpd(self):
        pass

    # TODO(2): May have to raise error. Calculating entire factor is costly.
    def to_factor(self):
        pass

    # TODO(1): Implement.
    def reorder_parents(self, new_order, inplace=True):
        pass

    # TODO(1): Implement.
    def get_evidence(self):
        pass

    # TODO: Implement. Unsure what this means in context. Note, this is static method.
    @staticmethod
    def get_random(variable, evidence=None, cardinality=None, state_name={}, seed=42):
        pass

    #
    # Operator overloading:
    #

    # TODO(2): Implement
    def __eq__(self, other, atol=1e-08):
        return NotImplemented

    # TODO: Implement?
    def __hash__(self):
        return NotImplemented
