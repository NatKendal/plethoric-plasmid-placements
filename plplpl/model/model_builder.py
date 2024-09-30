import pickle

from plplpl.base_functions import *
from plplpl.conjugation_functions import *
from plplpl.colour_delay_functions import *
from plplpl.maturation_delay_functions import *

from plplpl.NoisyOr import BayesianNetwork
from plplpl.NoisyOr import NoisyOrCPD


"""
Returns a NoisyOrBayesianNetwork

modelFolder: path to directory
dataFolder: path to directory
modelName: unique name for this model/computation
maturation_min: minimum number of timesteps for maturation 
maturation_max: maximum number of timesteps for maturation
save: if we should save the model to a file (pickle)
debug: 0 = nothing, 1 = status, 2 = verbose
progressBar: if we should show a progress bar on long for loops
doCheck: if we should do a graph check to ensure everything works and is acyclic
"""
def setupGraph(modelFolder, dataFolder, modelName, maturation_min, maturation_max, save=True, debug=0, progressBar=False, doCheck=True):
    # Assert correctness of constants.
    assert maturation_min >= 1
    assert maturation_max >= maturation_min

    # Load progress bar if requested.
    if progressBar:
        import tqdm

    # Get list of human friendly names to use as variable names.
    with open(dataFolder + modelName + "_humanFriendlyName.pickle", "rb") as f:
        name = pickle.load(f)

    # Save the list of unique cell ids for later.
    cells = list(name.keys())
    
    if debug >= 1:
        print("Computing the vertices.")
    
    # Calculate the set of vertices.
    vertices = set()
    for cell in cells:
        # Add two nodes per cell per timestep.
        # g[cellId]_[step] is 1 if cellId has the gene at timestep step, else 0.
        # m[cellId]_[step] is 1 if cellId is matured (can conjugate) at timestep step, else 0.
        vertices.add("g"+name[cell])
        vertices.add("m"+name[cell])
        # vertices.add("c"+name[cell])
        # We previously included colour, but we have simplified the model in v2 and now use colour implicitly.

    if debug >= 1:
        print("Computing the edges.")

    # Calculate the list of edges.
    edges = list()

    # We have to add edges one timestep forward for gene and maturation.
    with open(dataFolder + modelName + "_forwardLinks.pickle", "rb") as f:
        forwardLinks = pickle.load(f)

    # We also add edges from the maturation node of a cell to the gene node of its neighbours every timestep.
    with open(dataFolder + modelName + "_neighbours.pickle", "rb") as f:
        neighbours = pickle.load(f)

    if progressBar:
        iterable = tqdm.tqdm(cells)
    else:
        iterable = cells
    # Get the outgoing edges for each cell.
    for cell in iterable:
        if progressBar:
            iterable.set_description(desc="Adding edges for " + str(name[cell]))
        # First: add the simple downward edges. (g -> g, m -> m)
        for child in forwardLinks[cell]:
            edges.append(("g"+name[cell], "g"+name[child]))
            edges.append(("m"+name[cell], "m"+name[child]))
        # Second: add inter-cell edges to neighbours. (m -> g)
        for neighbour in neighbours[cell]:
            edges.append(("m"+name[cell], "g"+name[neighbour]))
        # Third: add the intra-cell edges to children. (g -> m)
        children = forwardLinks[cell].copy()
        depth = 1
        while depth <= maturation_max:
            newChildren = []
            for child in children: # Effectively do a BFS, and track how deep we've gone.
                newChildren.extend(forwardLinks[child])
                if depth >= maturation_min and depth <= maturation_max: # If we're in the range, add an edge.
                    edges.append(("g"+name[cell], "m"+name[child]))
            children = newChildren
            depth += 1

        if debug >= 2:
            print("Cell " + name[cell] + " complete.")

    if debug >= 1:
        print("Instantiating BayesianNetwork.")

    model = BayesianNetwork(vertices=vertices, edges=edges)

    # Save additional graph properties.
    model.constants["maturation_min"] = maturation_min
    model.constants["maturation_max"] = maturation_max

    if debug >= 1:
        print("Doing a structure check.")

    if doCheck:
        assert model.checkStructure(progressBar=progressBar)

    if save:
        if debug >= 1:
            print("Saving model.")
        with open(modelFolder + modelName + "_model_None_None_None.pickle", "wb") as f:
            pickle.dump(model, f)

    return model


"""
Adds CPDs to a model.

modelFolder: path to directory
dataFolder: path to directory
modelName: unique name for this model/computation
modelExtension: model extension in the form `_[conjugation function]_[colour function]_[maturation function]` (don't include .pickle)
delayFunctionPickleFile: pickle file with the BaseDelayFunction() instance to use
save: if we should save the model to a file (pickle)
debug: 0 = nothing, 1 = status, 2 = verbose
progressBar: if we should show a progress bar on long for loops
loadedModel: if we should use a model already in memory instead of loading one.
"""
def addDelayFunctionToModel(modelFolder, dataFolder, modelName, modelExtension, delayFunctionPickleFile, save=True, debug=0, progressBar=False, loadedModel=None):
    if progressBar:
        import tqdm

    if debug >= 1:
        print("Loading model.")
    if loadedModel:
        model = loadedModel
    else:
        with open(modelFolder + modelName + "_model" + modelExtension + ".pickle", "rb") as f:
            model = pickle.load(f)

    if debug >= 1:
        print("Loading weight function.")
    with open(delayFunctionPickleFile, "rb") as f:
        delayFunction = pickle.load(f)

    if not delayFunction.checkConstants(model, debug=True):
        raise ValueError("Constant mismatch between model and delayFunction." )

    if debug >= 1:
        print("Computing CPDs and adding to model.")
    if progressBar:
        iterator = tqdm.tqdm(model)
    else:
        iterator = model
    for node in iterator:
        if node[0] not in delayFunction.nodePrefix: # Check if the node is valid, otherwise skip it.
            continue
        if debug >= 2:
            print("Computing CPD for " + node)
        if progressBar:
            iterator.set_description(desc="Working on " + node)
        timestep = int(node.split("_")[1])
        evidence = []
        evidence_noise = []
        for parent in model.parents(node):
            evidence.append(parent)
            if parent[0] not in delayFunction.nodePrefix:
                evidence_noise.append(delayFunction.weight(timestep - int(parent.split("_")[1])))
            else:
                evidence_noise.append(1.0)

        model.add_cpd(node, NoisyOrCPD(node, baseChance=0, evidence=evidence, evidence_noise=evidence_noise))

    if debug >= 1:
        print("Finished adding CPDs.")
    
    model.constants["function" + str(delayFunction.nameIndex)] = delayFunction.name

    if save:
        if debug >= 1:
            print("Saving model.")
        with open(modelFolder + modelName + "_model" + delayFunction.injectName(modelExtension) + ".pickle", "wb") as f:
            pickle.dump(model, f)

    return model

"""
Add conjugation function CPDs to a model.

modelFolder: path to directory
dataFolder: path to directory
modelName: unique name for this model/computation
modelExtension: model extension in the form `_[conjugation function]_[colour function]_[maturation function]` (don't include .pickle)
conjugationFunctionPickleFile: pickle file with the BaseDelayFunction() instance to use
save: if we should save the model to a file (pickle)
debug: 0 = nothing, 1 = status, 2 = verbose
progressBar: if we should show a progress bar on long for loops
loadedModel: if we should use a model already in memory instead of loading one.
saveConjugationFunction: if we should save the conjugation function with loaded data for use later.
"""
def addConjugationFunctionToModel(modelFolder, dataFolder, modelName, modelExtension, conjugationFunctionPickleFile, save=True, debug=0, progressBar=False, loadedModel=None, saveConjugationFunction=True):
    if progressBar:
        import tqdm

    if debug >= 1:
        print("Loading model.")
    if loadedModel:
        model = loadedModel
    else:
        with open(modelFolder + modelName + "_model" + modelExtension + ".pickle", "rb") as f:
            model = pickle.load(f)

    if debug >= 1:
        print("Loading node unique ids.")
    with open(dataFolder + modelName + "_humanFriendlyNameLookup.pickle", "rb") as f:
        uid = pickle.load(f)

    if debug >= 1:
        print("Loading backward links.")
    with open(dataFolder + modelName + "_backwardLinks.pickle", "rb") as f:
        backwardLinks = pickle.load(f)

    if debug >= 1:
        print("Loading conjugation function.")
    with open(conjugationFunctionPickleFile, "rb") as f:
        conjugationFunction = pickle.load(f)

    if debug >= 1:
        print("Loading required data.")
    conjugationFunction.loadData(dataFolder, modelName, debug=debug)

    if debug >= 1:
        print("Computing CPDs and adding to model.")
    if progressBar:
        iterator = tqdm.tqdm(model)
    else:
        iterator = model
    for node in iterator:
        if node[0] not in conjugationFunction.nodePrefix:
            continue
        if debug >= 2:
            print("Computing CPD for " + node)
        if progressBar:
            iterator.set_description(desc="Working on " + node)

        evidence = []
        evidence_noise = []
        for parent in list(model.parents(node)):
            if parent[0] == "m":
                weight = conjugationFunction.weight(uid[parent[1:]], uid[node[1:]], debug=debug)
                if weight != 0:
                    evidence.append(parent)
                    evidence_noise.append(weight)
                else:
                    model.remove_edge((parent, node))
            else:
                evidence.append(parent)
                evidence_noise.append(1.0)

        model.add_cpd(node, NoisyOrCPD(node, baseChance=0, evidence=evidence, evidence_noise=evidence_noise))


    if debug >= 1:
        print("Finished adding CPDs.")

    if saveConjugationFunction:
        if debug >= 1:
            print("Saving conjugation function to " + dataFolder + modelName + "_" + conjugationFunction.name + ".pickle")
        with open(dataFolder + modelName + "_" + conjugationFunction.name + ".pickle", "wb") as f:
            pickle.dump(conjugationFunction, f)

    model.constants["function" + str(conjugationFunction.nameIndex)] = conjugationFunction.name

    if save:
        if debug >= 1:
            print("Saving model.")
        with open(modelFolder + modelName + "_model" + conjugationFunction.injectName(modelExtension) + ".pickle", "wb") as f:
            pickle.dump(model, f)

    return model

