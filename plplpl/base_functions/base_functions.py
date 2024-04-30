import pickle

"""
Basic class to package a delay function. (Colour and Maturation.)

name: function name (no underscores)
nameIndex: index in model extension where we insert the function name _[index 0]_[index 1]_[index 2]
nodePrefix: prefix for nodes for which this is an applicable function
constants: constants to assert correctness
lower: lower bound for timestep deltas
upper: upper bound for timestep deltas
values: dictionary with weights for valid timestep deltas
"""
class BaseDelayFunction(object):
    def __init__(self, name, nameIndex, nodePrefix, constants, lower, upper, values):
        self.name = name.replace("_", "")
        self.nameIndex = nameIndex
        self.nodePrefix = nodePrefix
        self.constants = constants
        self.lower = lower
        self.upper = upper
        self.values = values.copy()

    def checkConstants(self, model, forceExists=False):
        for key in self.constants.keys():
            if key in model.constants:
                if self.constants[key] != model.constants[key]:
                    return False
            elif forceExists:
                return False
        return True

    def weight(self, delta, **kwargs):
        if delta < self.lower:
            return None
        elif delta > self.upper:
            return None
        else:
            return self.values[delta]

    def injectName(self, modelExtension):
        return "_".join(modelExtension.split("_")[:self.nameIndex+1] + [self.name] + modelExtension.split("_")[self.nameIndex+2:])

    def save(self, filePath):
        with open(filePath, "wb") as f:
            pickle.dump(self, f)

"""
Basic class to package a conjugation function.

name: function name (no underscores)
nameIndex: index in model extension where we insert the function name
nodePrefix: prefix for nodes for which this is an applicable function
constants: constants to assert correctness
parameters: cell properties required for weight function calculation
"""
class BaseConjugationFunction(object):
    def __init__(self, name, nameIndex, nodePrefix, parameters):
        self.name = name.replace("_", "")
        self.nameIndex = nameIndex
        self.nodePrefix = nodePrefix
        self.parameters = parameters
        self.data = dict()

    def checkConstants(self, model, forceExists=False):
        for key in self.constants.keys():
            if key in model.constants:
                if self.constants[key] != model.constants[key]:
                    return False
                elif forceExists:
                    return False
        return True
    
    def loadData(self, dataFolder, modelName, debug=0):
        for parameter in self.parameters:
            with open(dataFolder + modelName + "_" + parameter + ".pickle", "rb") as f:
                self.data[parameter] = pickle.load(f)

    def weight(self, cell1, cell2, debug=0, **kwargs):
        # cell1parameters[param] is defined for param in self.parameters.
        return 0

    def injectName(self, modelExtension):
        return "_".join(modelExtension.split("_")[:self.nameIndex-1] + [self.name] + modelExtension.split("_")[self.nameIndex:])
    
    def save(self, filePath):
        with open(filePath, "wb") as f:
            pickle.dump(self, f)


