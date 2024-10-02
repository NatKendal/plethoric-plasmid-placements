import pickle

from plplpl.NoisyOr import BayesianNetwork

"""
modelFolder: path to directory
dataFolder: path to directory
modelName: unique name for this model/computation
modelExtension: model extension in the form `_[conjugation function]_[colour function]_[maturation function]` (don't include .pickle)
normalizeTo: factor to normalize the sum of edge weights to. Defaults to 1.0
save: if we should save the model to a file (pickle)
debug: 0 = nothing, 1 = status, 2 = verbose
progressBar: if we should show a progress bar on long for loops
loadedModel: if we should use a model already in memory instead of loading one.
edgeList: if we should normalize only particular edges instead of the entire graph.
    - (NOTE: This would have to be the same between different runs in order to have an unbiased comparison.)
    - (It could be accomplished by taking the union of all non-trivial edges across different runs.)
skipTrivial: if we should skip edges into the gene node of red and yellow cells.
    - We can't trim off anything else, since different maturation/colour ranges could change the edge list.
    - Edges into yellow or red cells are always trivial, because we assume it takes at least one step to change colour.
"""
# To fairly compare different conjugation functions based on their relative merit,
# we normalize the total conjugation edge weight to one.
# We could consider discounting trivial edges, but that limits comparisons between colour and maturation functions.
def normalizeConjugation(modelFolder, dataFolder, modelName, modelExtension, modelExtensionExtra, normalizeTo=1.0, save=True, debug=0, progressBar=False, loadedModel=None, edgeList=None, skipTrivial=True):
    if progressBar:
        import tqdm

    if debug >= 1:
        print("Loading model.")
    if loadedModel:
        model = loadedModel
    else:
        with open(modelFolder + modelName + "_model" + modelExtension + modelExtensionExtra + ".pickle", "rb") as f:
            model = pickle.load(f)

    if debug >= 1:
        print("Loading cell names and unique ids.")
    with open(dataFolder + modelName + "_humanFriendlyNameLookup.pickle", "rb") as f:
        uid = pickle.load(f)
    with open(dataFolder + modelName + "_humanFriendlyName.pickle", "rb") as f:
        name = pickle.load(f)

    if debug >= 1:
        print("Loading cell names and unique ids.")
    with open(dataFolder + modelName + "_humanFriendlyNameLookup.pickle", "rb") as f:
        uid = pickle.load(f)
    with open(dataFolder + modelName + "_humanFriendlyName.pickle", "rb") as f:
        name = pickle.load(f)

    if debug >= 1:
        print("Loading colours.")
    with open(dataFolder + modelName + "_colours.pickle", "rb") as f:
        colours = pickle.load(f)

    if debug >= 1:
        print("Getting edge list.")
    if edgeList == None:
        edgeList = set()
        if progressBar:
            iterator = tqdm.tqdm(model)
        else:
            iterator = model
        for vertex in iterator:
            if progressBar:
                iterator.set_description(desc="Working on " + str(vertex))
            if vertex[0] == "g":
                if colours[uid[vertex[1:]]] != 1:
                    continue
                for parent in model.parents(vertex):
                    if parent[0] == "m":
                        edgeList.add((parent, vertex))

    if debug >= 1:
        print("Starting to calculate total edge weight.")

    totalWeight = 0.0

    if progressBar:
        iterator = tqdm.tqdm(model.vertices())
    else:
        iterator = model.vertices()
    for vertex in iterator:
        if debug >= 2:
            print("Working on " + str(vertex))
        if progressBar:
            iterator.set_description(desc="Working on " + str(vertex))
        if vertex[0] == "g":
            for i in range(len(model.cpd(vertex)._evidence)):
                if (model.cpd(vertex)._evidence[i], vertex) in edgeList:
                    totalWeight += model.cpd(vertex)._evidence_noise[i]

    normalizationFactor = normalizeTo/totalWeight

    if debug >= 1:
        print("Total Weight: " + str(totalWeight))
        print("Normalization factor: " + str(normalizationFactor))
        print("Normalizing all transconjugant edges.")

    if progressBar:
        iterator = tqdm.tqdm(model.vertices())
    else:
        iterator = model.vertices()
    for vertex in iterator:
        if debug >= 2:
            print("Working on " + str(vertex))
        if progressBar:
            iterator.set_description(desc="Working on " + str(vertex))
        if vertex[0] == "g":
            for i in range(len(model.cpd(vertex)._evidence)):
                #if (model.cpd(vertex)._evidence[i], vertex) in edgeList:
                model.cpd(vertex).update_evidence_noise(i, normalizationFactor * model.cpd(vertex)._evidence_noise[i])

    if debug >= 1:
        print("Finished normalizing.")

    if save:
        if debug >= 1:
            print("Saving model.")
        with open(modelFolder + modelName + "_model" + modelExtension + "_conjugationNormalized.pickle", "wb") as f:
            pickle.dump(model, f)

    return model


"""
modelFolder: path to directory
dataFolder: path to directory
modelName: unique name for this model/computation
modelExtension: model extension in the form `_[conjugation function]_[colour function]_[maturation function]` (don't include .pickle)
normalizeTo: factor to normalize the sum of edge weights to. Defaults to 1.0
save: if we should save the model to a file (pickle)
debug: 0 = nothing, 1 = status, 2 = verbose
progressBar: if we should show a progress bar on long for loops
loadedModel: if we should use a model already in memory instead of loading one.
edgeList: if we should normalize only particular edges instead of the entire graph.
    - (NOTE: This would have to be the same between different runs in order to have an unbiased comparison.)
    - (It could be accomplished by taking the union of all non-trivial edges across different runs.)
skipTrivial: if we should skip edges into the gene node of red and yellow cells.
    - We can't trim off anything else, since different maturation/colour ranges could change the edge list.
    - Edges into yellow or red cells are always trivial, because we assume it takes at least one step to change colour.
"""
# A second slightly different normalization step designed to help moderate the impact of maturation functions that mature faster.
# More matured cells is strictly better, so we normalize edge weights a second time based on the naive probability calculation.
def normalizeMaturation(modelFolder, dataFolder, modelName, modelExtension, modelExtensionExtra, naiveProbabilitiesExtension, normalizeTo=1.0, save=True, debug=0, progressBar=False, loadedModel=None, edgeList=None, skipTrivial=True):
    if progressBar:
        import tqdm

    if debug >= 1:
        print("Loading model.")
    if loadedModel:
        model = loadedModel
    else:
        with open(modelFolder + modelName + "_model" + modelExtension + modelExtensionExtra + ".pickle", "rb") as f:
            model = pickle.load(f)

    if debug >= 1:
        print("Loading cell names and unique ids.")
    with open(dataFolder + modelName + "_humanFriendlyNameLookup.pickle", "rb") as f:
        uid = pickle.load(f)
    with open(dataFolder + modelName + "_humanFriendlyName.pickle", "rb") as f:
        name = pickle.load(f)

    if debug >= 1:
        print("Loading cell names and unique ids.")
    with open(dataFolder + modelName + "_humanFriendlyNameLookup.pickle", "rb") as f:
        uid = pickle.load(f)
    with open(dataFolder + modelName + "_humanFriendlyName.pickle", "rb") as f:
        name = pickle.load(f)

    if debug >= 1:
        print("Loading colours.")
    with open(dataFolder + modelName + "_colours.pickle", "rb") as f:
        colours = pickle.load(f)

    if debug >= 1:
        print("Loading naive probabilities.")

    with open(modelFolder + modelName + "_modeldata" + naiveProbabilitiesExtension + modelExtensionExtra + "_naiveProbabilities.pickle", "rb") as f:
        naiveProbabilities = pickle.load(f)

    if debug >= 1:
        print("Getting edge list.")
    if edgeList == None:
        edgeList = set()
        if progressBar:
            iterator = tqdm.tqdm(model)
        else:
            iterator = model
        for vertex in iterator:
            if progressBar:
                iterator.set_description(desc="Working on " + str(vertex))
            if vertex[0] == "g":
                if colours[uid[vertex[1:]]] != 1:
                    continue
                for parent in model.parents(vertex):
                    if parent[0] == "m":
                        edgeList.add((parent, vertex))

    if debug >= 1:
        print("Starting to calculate adjusted total edge weight.")

    totalWeight = 0.0

    if progressBar:
        iterator = tqdm.tqdm(model)
    else:
        iterator = model
    for vertex in iterator:
        if debug >= 2:
            print("Working on " + str(vertex))
        if progressBar:
            iterator.set_description(desc="Working on " + str(vertex))
        if vertex[0] == "g":
            for i in range(len(model.cpd(vertex)._evidence)):
                if (model.cpd(vertex)._evidence[i], vertex) in edgeList:
                    totalWeight += naiveProbabilities[model.cpd(vertex)._evidence[i]] * model.cpd(vertex)._evidence_noise[i]

    normalizationFactor = normalizeTo/totalWeight

    if debug >= 1:
        print("Total Weight: " + str(totalWeight))
        print("Normalization factor: " + str(normalizationFactor))
        print("Normalizing all transconjugant edges.")

    if progressBar:
        iterator = tqdm.tqdm(model)
    else:
        iterator = model
    for vertex in iterator:
        if debug >= 2:
            print("Working on " + str(vertex))
        if progressBar:
            iterator.set_description(desc="Working on " + str(vertex))
        if vertex[0] == "g":
            for i in range(len(model.cpd(vertex)._evidence)):
                #if (model.cpd(vertex)._evidence[i], vertex) in edgeList:
                model.cpd(vertex).update_evidence_noise(i, normalizationFactor * model.cpd(vertex)._evidence_noise[i])

    if debug >= 1:
        print("Finished normalizing.")

    if save:
        if debug >= 1:
            print("Saving model.")
        with open(modelFolder + modelName + "_model" + naiveProbabilitiesExtension + "_maturationNormalized.pickle", "wb") as f:
            pickle.dump(model, f)

    return model
