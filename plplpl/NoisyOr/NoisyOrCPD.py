class NoisyOrCPD(object):
    """
    Required Arguments:
        variable: variable over which the CPD is defined.
    Optional Arguments:
        baseChance: chance of the variable being 1 regardless of evidence.
        evidence: list of variables on which the CPD depends.
        evidence_noise: list of noise values for the links to the evidence. Must be the same length as evidence.
    """
    def __init__(self, variable, baseChance=0, evidence=list(), evidence_noise=list()):
        # Initialize variable and base chance of being 1.
        self._variable = variable
        self._baseChance = baseChance

        # Initialize evidence.
        if not isinstance(evidence, list):
            raise TypeError("Expected evidence to be list but got " + str(type(evidence)) + ".")
        self._evidence = evidence
        if not isinstance(evidence_noise, list):
            raise TypeError("Expected evidence_noise to be list but got " + str(type(evidence_noise)) + ".")
        self._evidence_noise = evidence_noise
        if len(self._evidence) != len(self._evidence):
            raise ValueError("Expected evidence and evidence_noise to be the same length.")

        # Setup set of all variables.
        self.variables = [self._variable] + self._evidence
        self._savedValues = dict() # Saves probability of variable = 1 given tuple of evidence assignments.

    # Calculate a value given an assignment to the variable and evidence.
    # Evidence variables can have uncertain probability and it is handled appropriately.
    def get_value(self, **kwargs):
        # Check that we have everything we need in kwargs.
        for var in self.variables:
            if var not in kwargs:
                raise ValueError("**kwargs missing variable " + str(var) + ".")

        # Check that the value of variable in integral.
        if not isinstance(kwargs[self._variable], int):
            raise ValueError("Assignment to " + str(self._variable) + " must be 0 or 1.")

        # Check if we already computed this value.
        assignment = tuple([kwargs[var] for var in self._evidence])
        if assignment in self._savedValues:
            if kwargs[self._variable] == 1:
                return self._savedValues[assignment]
            elif kwargs[self._variable] == 0:
                return (1.0 - self._savedValues[assignment])
            else:
                raise ValueError("Assignment to " + str(self._variable) + " wasn't 0 or 1.")

        # Otherwise, we compute a new value.
        # The chance of being zero is the chance that every cause fails.
        zeroChance = (1.0 - self._baseChance)
        for i in range(len(self._evidence)):
            # Chance of failure to transmit is 1 - (evidenceProbability * evidence_noise).
            zeroChance *= 1.0 - (kwargs[self._evidence[i]] * self._evidence_noise[i])
        
        # Store the computed result.
        self._savedValues[assignment] = (1.0 - zeroChance) # Chance of one is 1-zeroChance.

        # Return the correct result based on what we wanted from kwargs.
        if kwargs[self._variable] == 1:
            return (1.0 - zeroChance)
        else:
            return zeroChance

    # Remove a specific node from evidence.
    # WARNING: If this is changed, other objects that have precomputed a result involving this may not be updated.
    def remove_evidence(self, vertex):
        # Check that the variable is in our lists.
        try:
            i = self._evidence.index(vertex)
        except ValueError as e:
            raise ValueError(str(vertex) + " is not an evidence node.") from e

        # Remove the given vertex from evidence and variables.
        self._evidence = self._evidence[:i] + self._evidence[(i+1):]
        self._evidence_noise = self._evidence_noise[:i] + self._evidence_noise[(i+1):]
        self.variables = self.variables[:(i+1)] + self.variables[(i+2):]

        # Save the savedValues that didn't involve the vertex and delete the rest.
        newSavedValues = dict()
        for assignment in self._savedValues.keys():
            if assignment[i] == 0: # If assignment i was zero, then removing that edge doesn't change anything.
                newSavedValues[assignment[:i] + assignment[(i+1):]] = self._savedValues[assignment]

        # Overwrite the old values.
        self._savedValues = newSavedValues

    # Add a new evidence variable to the CPD.
    # WARNING: If this is changed, other objects that have precomputed a result involving this may not be updated.
    def add_evidence(self, vertex, noise):
        # Check if the variable already exists.
        if vertex in self.variables:
            raise ValueError(str(vertex) + " already exists in this CPD.")

        # Add it to the end of the internal lists.
        self.variables.append(vertex)
        self._evidence.append(vertex)
        self._evidence_noise.append(noise)

        # Update the saved values.
        newSavedValues = dict()
        for assignment in self._savedValues.keys():
            newSavedValues[assignment + (0,)] = self._savedValues[assignment]
        self._savedValues = newSavedValues

    # Safely change evidence_noise.
    # WARNING: If this is changed, other objects that have precomputed a result involving this may not be updated.
    def update_evidence_noise(self, i, evidence_noise):
        # Update the evidence noise.
        self._evidence_noise[i] == evidence_noise
        # Clear out saved values that may be invalid now.
        newSavedValues = dict()
        for assignment in self._savedValues.keys():
            if assignment[i] == 0: # If it was zero, then the noise didn't change anything.
                newSavedValues[assignment] = self._savedValues[assignment]
        self._savedValues = newSavedValues

    # Return the evidence noise for a particular evidence variable.
    def get_evidence_noise(self, evidence):
        index = self._evidence.index(evidence)
        return self._evidence_noise[index]


    # Change the base chance for the variable being one.
    # WARNING: If this is changed, other objects that have precomputed a result involving this may not be updated.
    def updateBaseChance(self, baseChance):
        self._baseChance = baseChance
        self._savedValues = dict() # Wipe the saved values for simplicity, we could try to reuse but it's annoying.
        
        

