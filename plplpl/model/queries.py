import pickle

from plplpl.NoisyOr import BayesianNetwork

"""
modelFolder: path to directory
dataFolder: path to directory
modelName: unique name for this model/computation
modelExtension: model extension in the form `_[conjugation function]_[colour function]_[maturation function]` (don't include _contradictionsPruned.pickle)
save: if we should save the model to a file (pickle)
debug: 0 = nothing, 1 = status, 2 = verbose
progressBar: if we should show a progress bar on long for loops
loadedModel: if we should use a model already in memory instead of loading one.
"""
def find_queries(modelFolder, dataFolder, modelName, modelExtension, modelExtensionExtra, save=True, debug=0, progressBar=False, loadedModel=None):
    if progressBar:
        import tqdm

    if debug >= 1:
        print("Loading model.")
    if loadedModel:
        model = loadedModel
    else:
        with open(modelFolder + modelName + "_model" + modelExtension + modelExtensionExtra + ".pickle", "rb") as f:
            model = pickle.load(f)

    # Get constants from model.
    colour_min = model.constants["colour_min"]
    colour_max = model.constants["colour_max"]
    maturation_min = model.constants["maturation_min"]
    maturation_max = model.constants["maturation_max"]

    if debug >= 1:
        print("Loading evidence.")
    with open(modelFolder + modelName + "_modeldata" + modelExtension + "_evidence.pickle", "rb") as f:
        evidence = pickle.load(f)

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

    directParent = model.getDirectParents()
    directChildren = model.getDirectChildren()
    
    if debug >= 1:
        print("Looking for lightup events.")

    # Find lightup events
    lightups = list()
    if progressBar:
        iterator = tqdm.tqdm(model)
    else:
        iterator = model
    for vertex in iterator:
        if vertex[0] != "g":
            continue
        # If the colour changed from green to yellow, this was a lightup event.
        if directParent[vertex] and (colours[uid[vertex[1:]]] == 2) and (colours[uid[directParent[vertex][1:]]] == 1):
            lightups.append(vertex)
            if progressBar:
                iterator.set_description(desc="Found new lightup at " + str(vertex))
            if debug >= 2:
                print("Found new lightup at " + str(vertex))

    #if debug >= 1:
    #    print("Found " + str(len(lightups)) + " lightups. Finding all queries.")

    if debug >= 1:
        print("Found " + str(len(lightups)) + " lightups. Finding all critical segments.")

    # Find all critical segments. (Chains of undefined nodes with no splits.)
    criticalSegmentSet = set()
    criticalSegmentMap = dict()
    seen = set()
    if progressBar:
        iterator = tqdm.tqdm(model)
    else:
        iterator = model
    for vertex in iterator:
        if vertex[0] != "g":
            continue
        if vertex in seen:
            continue
        if vertex in evidence:
            continue
        if progressBar:
            iterator.set_description(desc="Working on " + str(vertex))
        if debug >= 2:
            print("Working on " + str(vertex))

        segment = []
        segment.append(vertex)
        seen.add(vertex)
        parent = vertex
        while (directParent[parent]) and (directParent[parent] not in evidence) and (len(directChildren[directParent[parent]]) == 1):
            parent = directParent[parent]
            segment.append(parent)
            seen.add(parent)
        child = vertex
        while (len(directChildren[child]) == 1) and (directChildren[child][0] not in evidence):
            child = directChildren[child][0]
            segment.append(child)
            seen.add(child)
        sortedSegment = tuple(sorted(segment, key=lambda x: int(x.split("_")[1])))
        criticalSegmentSet.add(sortedSegment)
        for part in sortedSegment:
            criticalSegmentMap[part] = sortedSegment

    if debug >= 1:
        print("Found " + str(len(criticalSegmentSet)) + " critical segments. Building critical segment links.")

    criticalSegmentParents = dict()
    criticalSegmentChildren = dict()
    criticalSegmentLive = dict()
    criticalSegmentRoot = dict()

    if progressBar:
        iterator = tqdm.tqdm(criticalSegmentSet)
    else:
        iterator = criticalSegmentSet
    for segment in iterator:
        if progressBar:
            iterator.set_description(desc="Working on segment starting at " + str(segment[0]))
        if debug >= 2:
            print("Working on segment starting at " + str(segment[0]))
        for child in model.children(segment[-1]):
            if child[0] != "g":
                continue
            elif child in evidence:
                if evidence[child] == 0:
                    print("Something went wrong at " + str(child) + " " + str(segment))
                elif evidence[child] == 1:
                    if segment not in criticalSegmentChildren:
                        criticalSegmentChildren[segment] = []
                    criticalSegmentChildren[segment].append(child)
                    criticalSegmentLive[segment] = True
            else:
                for otherSegment in criticalSegmentSet:
                    if otherSegment == segment:
                        continue
                    if child in otherSegment:
                        if segment not in criticalSegmentChildren:
                            criticalSegmentChildren[segment] = []
                        criticalSegmentChildren[segment].append(otherSegment)
                        criticalSegmentParents[otherSegment] = segment

    if debug >= 1:
        print("Built critical segment links. Computing root segments.")

    if progressBar:
        iterator = tqdm.tqdm(criticalSegmentSet)
    else:
        iterator = criticlSegmentSet
    for segment in iterator:
        if progressBar:
            iterator.set_description(desc="Working on segment starting at " + str(segment[0]))
        if debug >= 2:
            print("Working on segment starting at " + str(segment[0]))
        parent = segment
        while parent in criticalSegmentParents:
            parent = criticalSegmentParents[parent]
        criticalSegmentRoot[segment] = parent

    if debug >= 1:
        print("Computed root segments. Computing live segments.")

    if progressBar:
        iterator = tqdm.tqdm(criticalSegmentSet)
    else:
        iterator = criticlSegmentSet
    for segment in iterator:
        if segment not in criticalSegmentLive:
            criticalSegmentLive[segment] = False
            continue
        if progressBar:
            iterator.set_description(desc="Working on segment starting at " + str(segment[0]))
        if debug >= 2:
            print("Working on segment starting at " + str(segment[0]))
        if criticalSegmentLive[segment] == True:
            parent = segment
            while parent in criticalSegmentParents:
                parent = criticalSegmentParents[parent]
                criticalSegmentLive[parent] = True

    if debug >= 1:
        print("Computed live segments. Finding query around each lightup event.")

    queries = dict()
    if progressBar:
        iterator = tqdm.tqdm(lightups)
    else:
        iterator = lightups
    for lightup in iterator:
        if progressBar:
            iterator.set_description(desc="Working on lightup " + str(lightup))
        if debug >= 2:
            print("Working on lightup " + str(lightup))

        # Setup query dictionary.
        queries[lightup] = dict()

        # Find critical region.
        critical = list()
        depth = 0
        parent = lightup
        while directParent[parent] and (depth < colour_max):
            parent = directParent[parent]
            depth += 1
            if (depth >= colour_min) and (depth <= colour_max):
                # Critical region.
                if parent not in evidence:
                    critical.append(parent)
                # Bottom of the critical region.
                elif (parent in evidence) and (evidence[parent] == 1) and directParent[parent]:
                    if (directParent[parent] not in evidence) or (evidence[directParent[parent]] == 0):
                        critical.append(parent)

        queries[lightup]["critical"] = sorted(critical, key=lambda x: int(x.split("_")[1]))

    """

    if debug >= 1:
        print("Computing incoming queries to each lightup event.")

    if progressBar:
        iterator = tqdm.tqdm(queries.keys())
    else:
        iterator = queries.keys()
    for query in iterator:
        if progressBar:
            iterator.set_description(desc="Working on lightup " + str(lightup))
        if debug >= 2:
            print("Working on lightup " + str(lightup))

        incomingFixed = set() # Maturation nodes in evidence.
        incomingVariable = set() # Maturation nodes not in evidence but have parents in a critical region.
        incomingCritical = set() # Incoming critical regions.
        incomingSpecial = set() # Maturation nodes not in evidence but also not from a critical region.

        for critical in queries[query]["critical"]:
            for maturation in [parent for parent in model.parents(critical) if parent[0] == "m"]:
                if maturation in evidence:
                    incomingFixed.add(maturation)
                else:
                    for gene in [parent for parent in model.parents(maturation) if parent[0] == "g"]:
                        if gene in evidence:
                            continue
                        for query in queries.keys():
                            if gene in queries[query]["critical"]:
                                incomingVariable.add(maturation)
                                incomingCritical.add(tuple(queries[query]["critical"]))
                                break
                        else:
                            incomingSpecial.add(maturation)

    """
    if debug >= 1:
        print("Finished building queries.")

    criticalSegments = {"set":criticalSegmentSet, "map":criticalSegmentMap, "root":criticalSegmentRoot, "live":criticalSegmentLive, "parents":criticalSegmentParents, "children":criticalSegmentChildren}

    if save:
        if debug >= 1:
            print("Saving queries and critical segments.")
        with open(modelFolder + modelName + "_modeldata" + modelExtension + "_queries.pickle", "wb") as f:
            pickle.dump(queries, f)
        with open(modelFolder + modelName + "_modeldata" + modelExtension + "_criticalSegments.pickle", "wb") as f:
            pickle.dump(criticalSegments, f)

    return queries, criticalSegments
