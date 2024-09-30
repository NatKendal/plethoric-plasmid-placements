from collections import deque

import pickle

def naiveCalcNode(model, directParents, naiveProbabilities, treatAsZero, node, previousProbability):
    kwargs = dict()
    missing = set()
    for parent in model.parents(node):
        if parent == directParents[node]:
            kwargs[parent] = previousProbability
        elif parent in treatAsZero:
            kwargs[parent] = 0
        else:
            if parent in naiveProbabilities:
                kwargs[parent] = naiveProbabilities[parent]
            else:
                kwargs[parent] = 0
                missing.add(parent)
    kwargs[node] = 1
    return model.cpd(node).get_value(**kwargs), missing

# Calculate chance of getting plasmid at each point, including probability of correct lightups.
def calcSegmentReceptionChances(model, directParents, directChildren, queries, colourFunction, naiveProbabilities, treatAsZero, relevantLightups, segment):
    segmentProbability = dict()
    missing = set()
    for node in segment:
        # Get the probability assuming that node's parent didn't have it.
        workingProbability, newMissing = naiveCalcNode(model, directParents, naiveProbabilities, treatAsZero, node, 0)
        missing = missing.union(newMissing)
        # Adjust by the probability of lightups happening at the given delta.
        for lightup in relevantLightups[segment]:
            workingProbability *= colourFunction.value(abs(int(lightup.split("_")[1]) - int(node.split("_")[1])))
        # Save the probability.
        segmentProbability[node] = workingProbability
    return segmentProbability, missing

# Recursive subprocess for calcSegmentFromRoot
def recursiveCalcSegmentHelper(model, directParents, queries, criticalSegments, colourFunction, naiveProbabilities, relevantLightups, live, treatAsZero, allSegmentProbabilities, precondition, segment, includeSelf):
    satisfactionProbability = 1.0
    missing = set()
    if includeSelf:
        for node in segment:
            precondition = precondition + ((1-precondition) * allSegmentProbabilities[node])
    if segment in criticalSegments["children"]:
        for child in criticalSegments["children"][segment]:
            if child in queries:
                prob, newMissing = naiveCalcNode(model, directParents, naiveProbabilities, treatAsZero, child, precondition)
                satisfactionProbability *= prob
                missing = missing.union(newMissing)
                for lightup in relevantLightups[segment]:
                    satisfactionProbability *= colourFunction.value(abs(int(lightup.split("_")[1]) - int(child.split("_")[1])))
            elif child in live:
                prob, newMissing = recursiveCalcSegmentHelper(model, directParents, queries, criticalSegments, colourFunction, naiveProbabilities, relevantLightups, live, treatAsZero, allSegmentProbabilities, precondition, child, True)
                satisfactionProbability *= prob
                missing = missing.union(newMissing)
    return satisfactionProbability, missing

# Check if some critical regions should be rescaled since we know they are satisfied.
def calcSegmentFromRoot(model, directParents, directChildren, evidence, queries, criticalSegments, colourFunction, naiveProbabilities, treatAsZero, rootSegment):
    colour_max = model.constants["colour_max"]
    # Given a root critical segment, first determine if it has children that have to be forced to 1.
    live = set()
    seen = set()
    endpoints = set()
    relevantLightups = dict()
    missing = set()
    queue = deque()
    queue.append((rootSegment, [rootSegment]))
    while queue:
        segment, path = queue.popleft()
        seen.add(segment)
        if segment in criticalSegments["children"]:
            for child in criticalSegments["children"][segment]:
                if (child in evidence) and (evidence[child] == 1):
                    endpoints.add(child)
                    # Start by identifying the lightup event(s) downstream from this certain gene node.
                    depth = 0
                    grandChildren = [child]
                    newLightups = []
                    while depth <= colour_max:
                        newGrandChildren = []
                        for grandChild in grandChildren:
                            for greatGrandChild in directChildren[grandChild]:
                                if greatGrandChild in queries.keys():
                                    newLightups.append(greatGrandChild)
                                else:
                                    newGrandChildren.append(greatGrandChild)
                        grandChildren = newGrandChildren
                        depth += 1
                    # Note each lightup event as relevant for each segment on the path towards here.
                    # Also note those paths as live.
                    for pathSegment in path:
                        live.add(pathSegment)
                        if pathSegment not in relevantLightups:
                            relevantLightups[pathSegment] = []
                        for lightup in newLightups:
                            relevantLightups[pathSegment].append(lightup)
                elif child in criticalSegments["set"]:
                    queue.append((child, path + [child]))
                else:
                    print("Something went wrong at " + str(rootSegment) + "/" + str(segment))

    # Handle the live segments, the ones that lead to a 1.
    # Compute the chances of reception at each live point.
    allSegmentProbabilities = dict()
    for segment in live:
        segmentProbability, newMissing = calcSegmentReceptionChances(model, directParents, directChildren, queries, colourFunction, naiveProbabilities, treatAsZero, relevantLightups, segment)
        missing = missing.union(newMissing)
        for node in segment:
            allSegmentProbabilities[node] = segmentProbability[node]

    # Calculate the new naive probabilities.
    newNaiveProbabilities = dict()
    queue = deque()
    queue.append((rootSegment, 0.0))
    while queue:
        segment, probability = queue.popleft()
        if segment in live: # Segments with a child that has evidence = 1
            for node in segment:
                probability = probability + ((1-probability) * allSegmentProbabilities[node])
                newNaiveProbabilities[node] = probability
            childrenSatisfactionChance, newMissing = recursiveCalcSegmentHelper(model, directParents, queries, criticalSegments, colourFunction, naiveProbabilities, relevantLightups, live, treatAsZero, allSegmentProbabilities, 0.0, segment, False)
            missing = missing.union(newMissing)
            constantModifier = newNaiveProbabilities[segment[-1]] + ((1-newNaiveProbabilities[segment[-1]]) * childrenSatisfactionChance)
            for node in segment:
                newNaiveProbabilities[node] = newNaiveProbabilities[node]/constantModifier
            if segment in criticalSegments["children"]:
                for child in criticalSegments["children"][segment]:
                    if child in evidence:
                        newNaiveProbabilities[child] = 1
                    else:
                        queue.append((child, newNaiveProbabilities[segment[-1]]))
        else: # Segments that don't lead to any evidence.
            for node in segment:
                prob, newMissing = naiveCalcNode(model, directParents, naiveProbabilities, treatAsZero, node, 0.0)
                missing = missing.union(newMissing)
                probability = probability + ((1-probability) * prob)
                newNaiveProbabilities[node] = probability
            if segment in criticalSegments["children"]:
                for child in criticalSegments["children"][segment]:
                    queue.append((child, newNaiveProbabilities[segment[-1]]))

    return newNaiveProbabilities, missing

def recursivelySolve(model, directParents, directChildren, evidence, queries, criticalSegments, colourFunction, segmentGroupMaturationNodes, naiveProbabilities, treatAsZero, rootSegment, depth, maxDepth):
    tempNaiveProbabilities = naiveProbabilities.copy()
    newTreatAsZero = treatAsZero.union(segmentGroupMaturationNodes[rootSegment])

    newNaiveProbabilities, missing = calcSegmentFromRoot(model, directParents, directChildren, evidence, queries, criticalSegments, colourFunction, tempNaiveProbabilities, newTreatAsZero, rootSegment)

    if depth < maxDepth:
        otherRoots = set()
        for node in missing:
            for otherRoot in segmentGroupMaturationNodes.keys():
                if node in segmentGroupMaturationNodes[otherRoot]:
                    otherRoots.add(otherRoot)
                    break
            else:
                raise ValueError("Failed to find missing value " + str(missing) + " when handling " + str(rootSegment) + ".")
        for otherRoot in otherRoots:
            recursiveNaiveProbabilities = recursivelySolve(model, directParents, directChildren, evidence, queries, criticalSegments, colourFunction, segmentGroupMaturationNodes, naiveProbabilities, newTreatAsZero, otherRoot, depth+1, maxDepth)
            tempNaiveProbabilities.update(recursiveNaiveProbabilities)

        newNaiveProbabilities, missing = calcSegmentFromRoot(model, directParents, directChildren, evidence, queries, criticalSegments, colourFunction, tempNaiveProbabilities, newTreatAsZero, rootSegment)
    
    for maturation in sorted(segmentGroupMaturationNodes[rootSegment], key=lambda x: int(x.split("_")[1])):
        parentVal = 0.0
        if directParents[maturation]:
            if directParents[maturation] in newNaiveProbabilities:
                parentVal = newNaiveProbabilities[directParents[maturation]]
            elif directParents[maturation] in naiveProbabilities:
                parentVal = naiveProbabilities[directParents[maturation]]
            else:
                print("SOMETHING WENT WRONG!")
                breakpoint()
        newNaiveProbabilities[maturation], _ = naiveCalcNode(model, directParents, newNaiveProbabilities, treatAsZero, maturation, parentVal)
    return newNaiveProbabilities

def computeNaiveProbabilities(modelFolder, dataFolder, modelName, modelExtension, modelExtensionExtra, colourFunctionPickleFile, depth=1, save=True, debug=0, progressBar=False, loadedModel=None):
    if progressBar:
        import tqdm

    if debug >= 1:
        print("Loading model.")
    if loadedModel:
        model = loadedModel
    else:
        with open(modelFolder + modelName + "_model" + modelExtension + modelExtensionExtra + ".pickle", "rb") as f:
            model = pickle.load(f)

    directParents = model.getDirectParents()
    directChildren = model.getDirectChildren()

    if debug >= 1:
        print("Loading evidence.")
    with open(modelFolder + modelName + "_modeldata" + modelExtension + "_evidence.pickle", "rb") as f:
        evidence = pickle.load(f)

    if debug >= 1:
        print("Loading queries.")
    with open(modelFolder + modelName + "_modeldata" + modelExtension + "_queries.pickle", "rb") as f:
        queries = pickle.load(f)

    if debug >= 1:
        print("Loading critical segments.")
    with open(modelFolder + modelName + "_modeldata" + modelExtension + "_criticalSegments.pickle", "rb") as f:
        criticalSegments = pickle.load(f)

    if debug >= 1:
        print("Loading colour function.")
    with open(colourFunctionPickleFile, "rb") as f:
        colourFunction = pickle.load(f)

    if not colourFunction.checkConstants(model):
        print("CONSTANTS CHECK FAILED!")
        exit()

    if debug >= 1:
        print("Initializing naive probabilities and critical segment groups.")
    naiveProbabilities = dict()

    for vertex in evidence.keys():
        # Add evidence to naive probabilities, except if the evidence node is the end of a critical region.
        if (vertex[0] == "g") and (evidence[vertex] == 1) and directParents[vertex] and (directParents[vertex] not in evidence):
            continue
        else:
            naiveProbabilities[vertex] = evidence[vertex]
    
    for vertex in sorted(model.vertices(), key=lambda x: int(x.split("_")[1])):
        if vertex[0] != "m":
            continue
        if vertex in evidence:
            continue
        defined = True
        for parent in model.parents(vertex):
            if parent in naiveProbabilities:
                continue
            else:
                defined = False
        if defined:
            naiveProbabilities[vertex], _ = naiveCalcNode(model, directParents, naiveProbabilities, set(), vertex, naiveProbabilities[directParents[vertex]] if directParents[vertex] else 0.0)

    # Get the root critical segments, then sort them from earliest for processing.
    criticalRoots = list()
    for segment in criticalSegments["set"]:
        if segment not in criticalSegments["parents"]:
            criticalRoots.append(segment)
    criticalRoots.sort(key=lambda x: int(x[0].split("_")[1]))

    # Calculate all the variable maturation nodes downstream from each critical segment group.
    segmentGroupMaturationNodes = dict()
    for rootSegment in criticalRoots:
        segmentGroupMaturationNodes[rootSegment] = set()
        queue = deque()
        queue.append(rootSegment)
        while queue:
            segment = queue.popleft()
            if segment in criticalSegments["children"]:
                for child in criticalSegments["children"][segment]:
                    if child in criticalSegments["set"]:
                        queue.append(child)
                    else:
                        for grandChild in model.children(child):
                            if (child[0] == "m") and (child not in evidence):
                                segmentGroupMaturationNodes[rootSegment].add(child)
            for node in segment:
                for child in model.children(node):
                    if (child[0] == "m") and (child not in evidence):
                        segmentGroupMaturationNodes[rootSegment].add(child)

    if debug >= 1:
        print("Finished initializing naive probabilities and critical segment groups. Now calculating naive probabilities.")

    if progressBar:
        iterator = tqdm.tqdm(criticalRoots)
    else:
        iterator = criticalRoots
    for rootSegment in iterator:
        if progressBar:
            iterator.set_description(desc="Working on critical segment starting at " + str(rootSegment[0]))
        if debug >= 2:
            print("Working on critical segment starting at " + str(rootSegment[0]))
        naiveProbabilities.update(recursivelySolve(model, directParents, directChildren, evidence, queries, criticalSegments, colourFunction, segmentGroupMaturationNodes, naiveProbabilities, set(), rootSegment, 0, depth))

    if debug >= 1:
        print("Finished calculating naive probabilities.")

    for vertex in model:
        if vertex not in naiveProbabilities:
            print("ERROR OCCURRED! MISSING " + vertex)

    if save:
        if debug >= 1:
            print("Saving naive probabilities.")
        modelExtension = colourFunction.injectName(modelExtension)
        with open(modelFolder + modelName + "_modeldata" + modelExtension + modelExtensionExtra + "_naiveProbabilities.pickle", "wb") as f:
            pickle.dump(naiveProbabilities, f)

    return naiveProbabilities
