import pickle

from plplpl.NoisyOr import NoisyOrBayesianNetwork

"""
modelFolder: path to directory
dataFolder: path to directory
modelName: unique name for this model/computation
modelExtension: model extension in the form `_[conjugation function]_[colour function]_[maturation function]` (don't include _contradictionsPruned.pickle)
conjugationFunctionPickleFile: full path to a conjugation function pickle file to calculate edge weights to get normalization.
normalizeTo: factor to normalize the sum of edge weights to. Defaults to 1.0
save: if we should save the model to a file (pickle)
debug: 0 = nothing, 1 = status, 2 = verbose
progressBar: if we should show a progress bar on long for loops
loadedModel: if we should use a model already in memory instead of loading one.
loadedEdgeList: if we should use an edge list already in memory instead of loading one.
loadedBackwardLinks: if we should use backward links already in memory instead of loading them.
"""

# For simplicity, we normalize all m -> g edges in the model. This could be fixed to only modify relevant edges if desired.
def normalizeModel(modelFolder, dataFolder, modelName, modelExtension, conjugationFunctionPickleFile, normalizeTo=1.0, save=True, debug=0, progressBar=False, loadedModel=None, loadedEdgeList=None, loadedBackwardLinks=None):
    if progressBar:
        import tqdm

    if debug >= 1:
        print("Loading model.")
    if loadedModel:
        model = loadedModel
    else:
        with open(modelFolder + modelName + "_model" + modelExtension + "_contradictionsPruned.pickle", "rb") as f:
            model = pickle.load(f)

    if debug >= 1:
        print("Loading cell names and unique ids.")
    with open(dataFolder + modelName + "_humanFriendlyNameLookup.pickle", "rb") as f:
        uid = pickle.load(f)
    with open(dataFolder + modelName + "_humanFriendlyName.pickle", "rb") as f:
        name = pickle.load(f)

    if debug >= 1:
        print("Loading edge list.")
    if loadedEdgeList:
        edgeList = loadedEdgeList
    else:
        with open(modelFolder + modelName + "_modeldata" + modelExtension + "_queryEdgeList.pickle", "rb") as f:
            edgeList = pickle.load(f)

    if debug >= 1:
        print("Loading backward links. (pruned)")
    if loadedBackwardLinks:
        backwardLinks = loadedBackwardLinks
    else:
        with open(modelFolder + modelName + "_modeldata" + modelExtension + "_backwardLinksPostEvidence.pickle", "rb") as f:
            backwardLinks = pickle.load(f)

    if debug >= 1:
        print("Loading conjugation function.")
    with open(conjugationFunctionPickleFile, "rb") as f:
        conjugationFunction = pickle.load(f)

    if debug >= 1:
        print("Starting to calculate total edge weights.")

    totalWeight = 0.0

    if progressBar:
        iterator = tqdm.tqdm(edgeList)
    else:
        iterator = edgeList
    for edge in iterator:
        if debug >= 2:
            print("Working on edge " + edge[0] + " to " + edge[1])
        if progressBar:
            iterator.set_description(desc="Working on edge " + edge[0] + " -> " + edge[1])
        totalWeight += conjugationFunction.weight(uid[edge[0][1:]], backwardLinks[uid[edge[1][1:]]], debug=debug)

    normalizationFactor = normalizeTo/totalWeight

    if debug >= 1:
        print("Total Weight: " + str(totalWeight))
        print("Normalization factor: " + str(normalizationFactor))
        print("Normalizing all transconjugant edges.")
 
    if progressBar:
        iterator = tqdm.tqdm(model.get_cpds())
    else:
        iterator = model.get_cpds()
    for cpd in iterator:
        if cpd.variable[0] != "g":
            continue
        if debug >= 2:
            print("Working on " + cpd.variable)
        if progressBar:
            iterator.set_description(desc="Working on " + cpd.variable)
        for i in range(len(cpd.evidence)):
            if cpd.evidence[i][0] == "m":
                cpd.evidence_noise[i] = cpd.evidence_noise[i] * normalizationFactor

    if debug >= 1:
        print("Finished normalizing.")

    if save:
        if debug >= 1:
            print("Saving model.")
        with open(modelFolder + modelName + "_model" + modelExtension + "_contradictionsPruned_normalized.pickle", "wb") as f:
            pickle.dump(model, f)

    return model
