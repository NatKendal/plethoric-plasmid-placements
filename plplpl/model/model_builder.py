import pickle

from plplpl.base_functions import *
from plplpl.conjugation_functions import *
from plplpl.colour_delay_functions import *
from plplpl.maturation_delay_functions import *

from plplpl.NoisyOr import NoisyOrBayesianNetwork
from plplpl.NoisyOr import BinaryNoisyOrCPD


"""
Returns a NoisyOrBayesianNetwork

modelFolder: path to directory
dataFolder: path to directory
modelName: unique name for this model/computation
colour_min: minimum number of timesteps for colour to appear
colour_max: maximum number of timesteps for colour to appear
maturation_min: minimum number of timesteps for maturation 
maturation_max: maximum number of timesteps for maturation
save: if we should save the model to a file (pickle)
debug: 0 = nothing, 1 = status, 2 = verbose
progressBar: if we should show a progress bar on long for loops
"""
#TESTING:
#c_min = 10
#c_max = 36
#m_min = 4
#m_max = 29
#mb.setupGraph("data/trap6test/", "data/trap6test/trap6test_data/", "trap6test", 10, 36, 4, 29, save=True, debug=1, progressBar=True)
def setupGraph(modelFolder, dataFolder, modelName, colour_min, colour_max, maturation_min, maturation_max, save=True, debug=0, progressBar=False):
    # Load progress bar if requested.
    if progressBar:
        import tqdm

    # Initialize a NoisyOrBayesianNetwork.
    model = NoisyOrBayesianNetwork()

    # Get list of human friendly names to use as variable names.
    with open(dataFolder + modelName + "_humanFriendlyName.pickle", "rb") as f:
        name = pickle.load(f)

    # Save the list of unique cell ids for later.
    cells = list(name.keys())
    
    if debug >= 1:
        print("Starting to add nodes to graph.")
    
    if progressBar:
        iterable = tqdm.tqdm(cells)
    else:
        iterable = cells
    # Add the nodes to the graph.
    for cell in iterable:
        # Add three nodes per cell per timestep.
        # g[cellId]_[step] is 1 if cellId has the gene at timestep step, else 0.
        # c[cellId]_[step] is 1 if cellId is coloured (not green) at timestep step, else 0.
        # m[cellId]_[step] is 1 if cellId is matured (can conjugate) at timestep step, else 0.
        model.add_nodes_from(["g"+name[cell], "c"+name[cell], "m"+name[cell]])
        if debug >= 2:
            print("Added cell " + name[cell])
        if progressBar:
            iterator.set_description(desc="Adding " + cell)

    # Add edges to the graph.
    # First: get the children of a cell, either true children or itself at the next timestep.
    # We have to add edges one timestep forward for the three properties, gene, colour, and maturation.
    with open(dataFolder + modelName + "_forwardLinks.pickle", "rb") as f:
        forwardLinks = pickle.load(f)
    # Second: get the neighbours of a cell.
    # We have to add edges one timestep forward for each possible mate.
    with open(dataFolder + modelName + "_neighbours.pickle", "rb") as f:
        neighbours = pickle.load(f)

    if debug >= 1:
        print("Finished adding nodes to graph, starting to add edges.")

    if progressBar:
        iterable = tqdm.tqdm(cells)
    else:
        iterable = cells
    # Add the edges, one cell at a time.
    for cell in iterable:
        if progressBar:
            iterator.set_description(desc="Adding edges to " + cell)
        # First: add the simple downward edges. (g -> g, c -> c, m -> m)
        for child in forwardLinks[cell]:
            model.add_edge("g"+name[cell], "g"+name[child])
            model.add_edge("c"+name[cell], "c"+name[child])
            model.add_edge("m"+name[cell], "m"+name[child])
        # Second: add inter-cell edges to neighbours. (m -> g)
        for neighbour in neighbours[cell]:
            model.add_edge("m"+name[cell], "g"+name[neighbour])
        # Third: add the intra-cell edges to children. (g -> c, g -> m)
        children = forwardLinks[cell].copy()
        depth = 1
        while depth <= max(colour_max, maturation_max):
            # If we're in the correct range for colour, add edges.
            if depth >= colour_min and depth <= colour_max:
                for child in children:
                    model.add_edge("g"+name[cell], "c"+name[child])
            # If we're in the correct range for maturation, add edges.
            if depth >= maturation_min and depth <= maturation_max:
                for child in children:
                    model.add_edge("g"+name[cell], "m"+name[child])
            # In any case, get all of the children's children and repeat until we hit the required depth.
            newChildren = []
            for child in children:
                newChildren.extend(forwardLinks[child])
            children = newChildren
            depth += 1

        if debug >= 2:
            print("Cell " + name[cell] + " complete.")

    if debug >= 1:
        print("Finished adding edges to graph.")

    # Save additional graph properties.
    model.constants["colour_min"] = colour_min
    model.constants["colour_max"] = colour_max
    model.constants["maturation_min"] = maturation_min
    model.constants["maturation_max"] = maturation_max

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
safeMode: if we should use the BayesianNetwork.add_cpds. It's slower but does checks.
loadedModel: if we should use a model already in memory instead of loading one.
"""
#TESTING:
#mb.addDelayFunctionToModel("data/trap6test/", "data/trap6test/trap6test_data/", "trap6test", "_None_None_None", "functions/colourDelayFunctionUniform.pickle", save=True, debug=1, progressBar=True, safeMode=False, loadedModel=None)
#mb.addDelayFunctionToModel("data/trap6test/", "data/trap6test/trap6test_data/", "trap6test", "_None_colourDelayFunctionUniform_None", "functions/maturationDelayFunctionNormal.pickle", save=True, debug=1, progressBar=True, safeMode=False, loadedModel=None)
def addDelayFunctionToModel(modelFolder, dataFolder, modelName, modelExtension, delayFunctionPickleFile, save=True, debug=0, progressBar=False, safeMode=False, loadedModel=None):
    if progressBar:
        import tqdm
        if debug >= 1:
            print("Loaded tqdm for progress bar.")

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

    if not delayFunction.checkConstants(model):
        raise ValueError("Constant mismatch between model and delayFunction.")

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
        for predecessor in model.predecessors(node):
            evidence.append(predecessor)
            if predecessor[0] == "g":
                evidence_noise.append(delayFunction.weight(timestep - int(predecessor.split("_")[1])))
            else:
                evidence_noise.append(1.0)
        # NOTE: We are forcefully adding CPDs to the model here. Trading safety for speed.
        # 8 hours -> 30 seconds
        if safeMode:
            model.add_cpds(BinaryNoisyOrCPD(node, [[1], [0]], evidence=evidence, evidence_noise=evidence_noise))
        else:
            model.cpds.append(BinaryNoisyOrCPD(node, [[1], [0]], evidence=evidence, evidence_noise=evidence_noise))

    if debug >= 1:
        print("Finished adding CPDs.")

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
safeMode: if we should use the BayesianNetwork.add_cpds. It's slower but does checks.
loadedModel: if we should use a model already in memory instead of loading one.
saveConjugationFunction: if we should save the conjugation function with loaded data for use later.
"""
#TESTING:
#mb.addConjugationFunctionToModel("data/trap6test/", "data/trap6test/trap6test_data/", "trap6test", "_None_colourDelayFunctionUniform_maturationDelayFunctionNormal", "functions/contactWeightsFixedNaive.pickle", save=True, debug=1, progressBar=True, safeMode=False, loadedModel=None, saveConjugationFunction=True)
def addConjugationFunctionToModel(modelFolder, dataFolder, modelName, modelExtension, conjugationFunctionPickleFile, save=True, debug=0, progressBar=False, safeMode=False, loadedModel=None, saveConjugationFunction=False):
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
        print("Loading parents.")
    with open(dataFolder + modelName + "_parent.pickle", "rb") as f:
        parent = pickle.load(f)

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
        for predecessor in model.predecessors(node):
            evidence.append(predecessor)
            if predecessor[0] == "m":
                evidence_noise.append(conjugationFunction.weight(parent[uid[node[1:]]], uid[predecessor[1:]], debug=debug))
            else:
                evidence_noise.append(1.0)

        # NOTE: We are forcefully adding CPDs to the model here. Trading safety for speed.
        if safeMode:
            model.add_cpds(BinaryNoisyOrCPD(node, [[1], [0]], evidence=evidence, evidence_noise=evidence_noise))
        else:
            model.cpds.append(BinaryNoisyOrCPD(node, [[1], [0]], evidence=evidence, evidence_noise=evidence_noise))

    if debug >= 1:
        print("Finished adding CPDs.")

    if saveConjugationFunction:
        if debug >= 1:
            print("Saving conjugation function to " + dataFolder + modelName + "_" + conjugationFunction.name + ".pickle")
        with open(dataFolder + modelName + "_" + conjugationFunction.name + ".pickle", "wb") as f:
            pickle.dump(conjugationFunction, f)

    if save:
        if debug >= 1:
            print("Saving model.")
        with open(modelFolder + modelName + "_model" + conjugationFunction.injectName(modelExtension) + ".pickle", "wb") as f:
            pickle.dump(model, f)

    return model

#import plplpl.model.model_builder as mb
#mb.addConjugationFunctionToModel("data/trap6test/", "data/trap6test/trap6test_data/", "trap6test", "_None_colourDelayFunctionUniform_maturationDelayFunctionNormal", "functions/contactWeightsFixedNaive.pickle", save=True, debug=1, progressBar=True, safeMode=False, loadedModel=mb.addDelayFunctionToModel("data/trap6test/", "data/trap6test/trap6test_data/", "trap6test", "_None_None_None", "functions/colourDelayFunctionUniform.pickle", save=True, debug=1, progressBar=True, safeMode=False, loadedModel=mb.addDelayFunctionToModel("data/trap6test/", "data/trap6test/trap6test_data/", "trap6test", "_None_colourDelayFunctionUniform_None", "functions/maturationDelayFunctionNormal.pickle", save=True, debug=1, progressBar=True, safeMode=False, loadedModel=mb.setupGraph("data/trap6test/", "data/trap6test/trap6test_data/", "trap6test", 10, 36, 4, 29, save=True, debug=1, progressBar=True))), saveConjugationFunction=False)


# c23621_103
