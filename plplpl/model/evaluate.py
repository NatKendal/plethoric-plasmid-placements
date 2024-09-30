import itertools
import pickle
import signal

from collections import deque

from plplpl.NoisyOr import NoisyOrCPD
from plplpl.NoisyOr import NoisyOrFactor

class TimeLimitException(Exception):
    pass

def assignmentHelper(model, directParents, directChildren, naiveProbabilities, groupChildren, root):
    zeros = []
    current = root
    assignments = dict()
    zeroSoFar = 1.0
    while (current in groupChildren) and (len(groupChildren[current]) == 1):
        prob = model.cpd(current).get_value(**{current:1}, **{parent:(naiveProbabilities[parent] if parent != directParents[current] else 0) for parent in model.parents(current)})
        assignments[tuple(zeros)] = zeroSoFar * prob
        zeros.append(current)

        child = current
        while child != groupChildren[current][0]:
            zeroSoFar *= model.cpd(child).get_value(**{child:0}, **{parent:(naiveProbabilities[parent] if parent != directParents[child] else 0) for parent in model.parents(child)})
            if len(directChildren[child]) > 1:
                queue = deque()
                for grandChild in directChildren[child]:
                    queue.append((grandChild, 1.0))
                while queue:
                    node, prob = queue.popleft()
                    if node == groupChildren[current][0]:
                        child = node
                        zeroSoFar *= prob
                        break
                    else:
                        for grandChild in directChildren[node]:
                            queue.append((grandChild, prob * model.cpd(node).get_value(**{node:0}, **{parent:(naiveProbabilities[parent] if parent != directParents[node] else 0) for parent in model.parents(node)})))
            else:
                child = directChildren[child][0]
        current = child
    if current not in groupChildren:
        prob = model.cpd(current).get_value(**{current:1}, **{parent:(naiveProbabilities[parent] if parent != directParents[current] else 0) for parent in model.parents(current)})
        assignments[tuple(zeros)] = zeroSoFar * prob
        zeros.append(current)
        assignments[tuple(zeros)] = zeroSoFar * (1-prob)
        return assignments
    elif len(groupChildren[current]) > 1:
        prob = model.cpd(current).get_value(**{current:1}, **{parent:(naiveProbabilities[parent] if parent != directParents[current] else 0) for parent in model.parents(current)})
        assignments[tuple(zeros)] = zeroSoFar * prob
        zeros.append(current)

        found = []
        queue = deque()
        queue.append((current, []))
        while (len(found) < len(groupChildren[current])):
            node, path = queue.popleft()
            if node in groupChildren[current]:
                found.append((node, path))
            else:
                for child in directChildren[node]:
                    queue.append((child, path + [node]))
        shared = set(found[0][1])
        for i in range(1, len(found)):
            shared = shared.intersection(found[i][1])
        for node in shared:
            zeroSoFar *= model.cpd(node).get_value(**{node:0}, **{parent:(naiveProbabilities[parent] if parent != directParents[node] else 0) for parent in model.parents(node)})

        subassignments = []
        for node, path in found:
            subassignment = assignmentHelper(model, directParents, directChildren, naiveProbabilities, groupChildren, node)
            specificZerosBefore = 1.0
            for previous in path:
                if previous in shared:
                    continue
                specificZerosBefore *= model.cpd(previous).get_value(**{previous:0}, **{parent:(naiveProbabilities[parent] if parent != directParents[previous] else 0) for parent in model.parents(previous)})
            subassignmentList = [(assignment, specificZerosBefore * subassignment[assignment]) for assignment in subassignment.keys()]
            subassignments.append(subassignmentList)

        allAssignments = itertools.product(*subassignments)
        for assignment in allAssignments:
            fullAssignment = zeros.copy()
            fullProb = zeroSoFar
            for nodes, prob in assignment:
                fullAssignment.extend(nodes)
                fullProb *= prob
            assignments[tuple(fullAssignment)] = fullProb
        return assignments
    else:
        raise RuntimeError("Something went very wrong in assignmentHelper.")

def evaluateQueryPoint(model, directParents, directChildren, evidence, queries, criticalSegments, naiveProbabilities, queryPoint, queryPointValues):
    # Fetch the query point CPD and reduce out where possible.
    allIncoming = set() # All incoming maturation nodes to lineage.
    argument = [(queryPoint, 1)]
    if directParents[queryPoint]:
        argument.append((directParents[queryPoint], 0))
    for incoming in model.parents(queryPoint):
        if incoming in evidence:
            if incoming[0] == "m":
                argument.append((incoming, evidence[incoming]))
        else:
            allIncoming.add(incoming)
    queryPointFactor = NoisyOrFactor("reduce", [model.cpd(queryPoint)], argument=argument)

    # Get all maturation nodes we will have to marginalize over.
    criticalMaturation = []
    criticalMaturationChildren = dict()
    querySpecialEvidence = dict()
    queryTime = int(queryPoint.split("_")[1])
    for child in model.children(queryPoint):
        if child[0] == "g":
            continue
        criticalMaturation.append(child)
        for parent in model.parents(child):
            if parent != "m":
                if int(parent.split("_")[1]) >= queryTime:
                    querySpecialEvidence[parent] = 1
                else:
                    querySpecialEvidence[parent] = 0
        relevantChildren = [grandChild for grandChild in directChildren[child] if grandChild in model.children(queryPoint)]
        if len(relevantChildren) > 0:
            criticalMaturationChildren[child] = relevantChildren
    criticalMaturation.sort(key=lambda x: int(x.split("_")[1]))

    # Calculate the fixed maturation nodes for this lineage.
    outgoingZeros = set()
    outgoingOnes = set()
    minTimestep = int(criticalMaturation[0].split("_")[1])
    maxTimestep = int(criticalMaturation[-1].split("_")[1])
    for maturation in criticalMaturation:
        if int(maturation.split("_")[1]) == minTimestep:
            parent = maturation
            while directParents[parent]:
                parent = directParents[parent]
                outgoingZeros.add(parent)
        if int(maturation.split("_")[1]) == maxTimestep:
            queue = deque()
            for child in directChildren[maturation]:
                queue.append(child)
            while queue:
                child = queue.popleft()
                outgoingOnes.add(child)
                for grandChild in directChildren[child]:
                    queue.append(grandChild)


    # Build the factor giving the chance of all previous gene nodes being 0.
    zerosFactors = [] # List of non-trivial CPDs
    zerosConstant = 1.0 # Product of CPDs that are constants.

    parent = queryPoint
    while directParents[parent]:
        parent = directParents[parent]

        # Build factor.
        argument = [(parent, 0)]
        if directParents[parent]:
            argument.append((directParents[parent], 0))
        for incoming in model.parents(parent):
            if incoming in evidence:
                if incoming[0] == "m":
                    argument.append((incoming, evidence[incoming]))
        factor = NoisyOrFactor("reduce", [model.cpd(parent)], argument=argument)
        if len(factor.variables) == 0:
            zerosConstant *= factor.get_value()
        else:
            zerosFactors.append(factor)
            for variable in factor.variables:
                allIncoming.add(variable)

    # Calculate the penalty from forced one maturation nodes pointing into zero nodes.
    # (Chance that a failure to transmit occurred where it must have.)
    forcedOneTransmitFailureProbability = 1.0 # Chance all possible transmissions failed.
    for node in outgoingOnes:
        for child in model.children(node):
            if (child[0] == "g") and (child in evidence) and (evidence[child] == 0):
                forcedOneTransmitFailureProbability *= (1-model.cpd(child).get_evidence_noise(node))

    # Build a factor to represent the penalty for criticalMaturation nodes pointing into zero nodes.
    criticalMaturationNoise = []
    for node in criticalMaturation:
        failureChance = 1.0
        for child in model.children(node):
            if (child[0] == "g") and (child in evidence) and (evidence[child] == 0):
                failureChance *= (1 - model.cpd(child).get_evidence_noise(node))
        # failureChance is probability that all transmissions failed.
        # So 1-failureChance is probability that at least one transmission succeeded.
        criticalMaturationNoise.append(1-failureChance)
    criticalMaturationTransmitFailureFactor = NoisyOrFactor("reduce", [NoisyOrCPD("temp", baseChance=0, evidence=criticalMaturation.copy(), evidence_noise=criticalMaturationNoise.copy())], argument=[("temp", 0)])


    # Locate critical segments relevant to our queryPoint.
    relevantSegments = set()
    singlePointGeneLightupPoints = set() # Edge case where there is no critical segment.
    for node in outgoingZeros.union(criticalMaturation):
        for child in model.children(node):
            if child[0] == "m":
                continue
            if (child not in evidence):
                segment = criticalSegments["map"][child]
                if criticalSegments["live"][segment]:
                    # If it doesn't lead to a definite gene node, then it can't impact the probability in any way.
                    relevantSegments.add(segment)
            elif (child in evidence) and (evidence[child] == 1) and directParents[child]:
                if (directParents[child] not in evidence):
                    segment = criticalSegments["map"][directParents[child]]
                elif evidence[directParents[child]] == 0:
                    singlePointGeneLightupPoints.add(child)

    # Locate transmission points relevant to our queryPoint.
    transmissionPoints = set()
    seen = set()
    for segment in relevantSegments:
        queue = deque()
        queue.append(segment)
        while queue:
            currentSegment = queue.popleft()
            if currentSegment in seen:
                continue
            seen.add(currentSegment)
            for child in criticalSegments["children"][currentSegment]:
                if child in evidence:
                    transmissionPoints.add(child)
                elif criticalSegments["live"][child]:
                    queue.append(child)
    # Locate the root of each transmission point.
    roots = set()
    for transmissionPoint in transmissionPoints:
        roots.add(criticalSegments["root"][criticalSegments["map"][directParents[transmissionPoint]]])

    # For each root, build collapsed CPDs for each live segment.
    depth1CriticalRegions = []
    for root in roots:
        if not criticalSegments["live"][root]:
            raise RuntimeError("Root node wasn't live, this shouldn't happen.")
        queue = deque()
        queue.append((root, None))
        tempFactors = []
        tempToMarginalize = []
        while queue:
            segment, previous = queue.popleft()
            if not criticalSegments["live"][segment]:
                raise RuntimeError("Dead segments should be filtered out. Something went wrong.")
            nodes = list(segment)
            while (segment in criticalSegments["children"]) and (len([x for x in criticalSegments["children"][segment] if x in evidence]) == 0) and (len([x for x in criticalSegments["children"][segment] if criticalSegments["live"][x]]) == 1):
                # There's only one relevant child, so expand through it.
                segment = [x for x in criticalSegments["children"][segment] if criticalSegments["live"][x]][0]
                nodes.extend(list(segment))
            if (segment in criticalSegments["children"]) and (len([x for x in criticalSegments["children"][segment] if (x in evidence) or (criticalSegments["live"][x])]) >= 2):
                # Build an intermediate CPD.
                nonTrivialIncoming = []
                nonTrivialIncomingNoise = []
                if previous:
                    nonTrivialIncoming.append(previous)
                    nonTrivialIncomingNoise.append(1.0)
                failureChance = 1.0
                for node in nodes:
                    cpd = model.cpd(node)
                    for parent in model.parents(node):
                        if parent[0] == "g":
                            continue
                        if (parent in criticalMaturation) or (parent in outgoingZeros) or (parent in allIncoming):
                            nonTrivialIncoming.append(parent)
                            nonTrivialIncomingNoise.append(cpd.get_evidence_noise(parent))
                        else:
                            failureChance *= (1 - (cpd.get_evidence_noise(parent) * naiveProbabilities[parent]))
                intermediateCPD = NoisyOrCPD(nodes[-1], baseChance=(1-failureChance), evidence=nonTrivialIncoming, evidence_noise=nonTrivialIncomingNoise)
                tempFactors.append(intermediateCPD)
                tempToMarginalize.append(nodes[-1])
                for childSegment in criticalSegments["children"][segment]:
                    if childSegment in evidence:
                        # Make a small terminal factor.
                        nonTrivialIncoming = [nodes[-1]]
                        nonTrivialIncomingNoise = [1.0]
                        failureChance = 1.0
                        cpd = model.cpd(childSegment)
                        for parent in model.parents(childSegment):
                            if parent[0] == "g":
                                continue
                            if (parent in criticalMaturation) or (parent in outgoingZeros) or (parent in allIncoming):
                                nonTrivialIncoming.append(parent)
                                nonTrivialIncomingNoise.append(cpd.get_evidence_noise(parent))
                            else:
                                failureChance *= (1 - (cpd.get_evidence_noise(parent) * naiveProbabilities[parent]))
                        terminalCPD = NoisyOrCPD(childSegment, baseChance=(1-failureChance), evidence=nonTrivialIncoming, evidence_noise=nonTrivialIncomingNoise)
                        tempFactors.append(NoisyOrFactor("reduce", [terminalCPD], argument=[(childSegment, 1)]))
                    else:
                        # Move down the tree.
                        if criticalSegments["live"][childSegment]:
                            queue.append((childSegment, nodes[-1]))

            else:
                # Build a terminal factor.
                terminal = [x for x in criticalSegments["children"][segment] if (x in evidence)][0]
                nonTrivialIncoming = []
                nonTrivialIncomingNoise = []
                if previous:
                    nonTrivialIncoming.append(previous)
                    nonTrivialIncomingNoise.append(1.0)
                failureChance = 1.0
                nodes.append(terminal)
                for node in nodes:
                    cpd = model.cpd(node)
                    for parent in model.parents(node):
                        if parent[0] == "g":
                            continue
                        if (parent in criticalMaturation) or (parent in outgoingZeros) or (parent in allIncoming):
                            nonTrivialIncoming.append(parent)
                            nonTrivialIncomingNoise.append(cpd.get_evidence_noise(parent))
                        else:
                            failureChance *= (1 - (cpd.get_evidence_noise(parent) * naiveProbabilities[parent]))
                terminalCPD = NoisyOrCPD(nodes[-1], baseChance=(1-failureChance), evidence=nonTrivialIncoming, evidence_noise=nonTrivialIncomingNoise)
                tempFactors.append(NoisyOrFactor("reduce", [terminalCPD], argument=[(nodes[-1], 1)]))
        depth1CriticalRegions.append((tempFactors, tempToMarginalize))
        # Future work:
        # Dynamically adjust this so that we model at greater depth if time permits.

    # Handle single point gene lightup edge cases.
    for node in singlePointGeneLightupPoints:
        nonTrivialIncoming = []
        nonTrivialIncomingNoise = []
        failureChance = 1.0
        cpd = model.cpd(node)
        for parent in model.parents(node):
            if parent[0] == "g":
                continue
            if (parent in criticalMaturation) or (parent in outgoingZeros) or (parent in allIncoming):
                nonTrivialIncoming.append(parent)
                nonTrivialIncomingNoise.append(cpd.get_evidence_noise(parent))
            else:
                failureChance *= (1 - (cpd.get_evidence_noise(parent) * naiveProbabilities[parent]))
        depth1CriticalRegions.append(([NoisyOrFactor("reduce", [NoisyOrCPD(node, baseChance=(1-failureChance), evidence=nonTrivialIncoming, evidence_noise=nonTrivialIncomingNoise)], argument=[(node, 1)])], []))

    # Construct factors that resolve to the penalty from pointing into a critical region.
    downstreamLightupFactors = []
    individualToMarginalize = set()
    for factorList, marginalizeList in depth1CriticalRegions:
        individualToMarginalize = individualToMarginalize.union(marginalizeList)
        for factor in factorList:
            onesArgument = [(variable, 1) for variable in factor.variables if (variable in criticalMaturation) or (variable in outgoingZeros)]
            zerosArgument = [(variable, 0) for variable in factor.variables if (variable in outgoingZeros)]
            downstreamLightupFactors.append(NoisyOrFactor("quotient", [NoisyOrFactor("reduce", [factor], argument=zerosArgument), NoisyOrFactor("reduce", [factor], argument=onesArgument)]))

    # Find maturation node groups.
    allIncomingParent = dict()
    allIncomingChildren = dict()
    maturationNodeGroups = []
    for node in allIncoming:
        parent = node
        while directParents[parent]:
            parent = directParents[parent]
            if parent in allIncoming:
                allIncomingParent[node] = parent
                if parent not in allIncomingChildren:
                    allIncomingChildren[parent] = []
                allIncomingChildren[parent].append(node)
                break
            elif parent in evidence:
                break
    for root in allIncoming:
        if root in allIncomingParent:
            continue
        # Otherwise, this is a root node.
        incomingGroup = []
        queue = deque()
        queue.append(root)
        while queue:
            node = queue.popleft()
            incomingGroup.append(node)
            if node in allIncomingChildren:
                for child in allIncomingChildren[node]:
                    queue.append(child)
        maturationNodeGroups.append(incomingGroup)    

    # Combine all factors.

    allFactors = [queryPointFactor, criticalMaturationTransmitFailureFactor] + zerosFactors + downstreamLightupFactors

    # Marginalize out over each maturation node group.
    for group in sorted(maturationNodeGroups, key=lambda x: len(x)):
        nodeSet = set(group)
        workingFactors = []
        newFactors = []
        for factor in allFactors:
            if len(nodeSet.intersection(factor.variables)) > 0:
                workingFactors.append(factor)
            else:
                newFactors.append(factor)
        if len(workingFactors) > 1:
            productFactor = NoisyOrFactor("product", workingFactors)
        elif len(workingFactors) == 1:
            productFactor = workingFactors[0]
        else:
            continue

        assignments = assignmentHelper(model, directParents, directChildren, naiveProbabilities, allIncomingChildren, group[0])
        reducedFactors = []
        for assignment in assignments.keys():
            reducedFactors.append(NoisyOrFactor("product", [NoisyOrFactor("constant", [], argument=[assignments[assignment]]), NoisyOrFactor("reduce", [productFactor], argument=[(node, 0 if node in assignment else 1) for node in group])]))
        sumFactor = NoisyOrFactor("sum", reducedFactors)
        newFactors.append(sumFactor)
        allFactors = newFactors

    # Marginalize out the remaining things to marginalize.
    for node in individualToMarginalize:
        workingFactors = []
        newFactors = []
        for factor in allFactors:
            if node in factor.variables:
                workingFactors.append(factor)
            else:
                newFactors.append(factor)
        if len(workingFactors) > 1:
            productFactor = NoisyOrFactor("product", workingFactors)
        elif len(workingFactors) == 1:
            productFactor = workingFactors[0]
        else:
            continue
        newFactors.append(NoisyOrFactor("marginalize", [productFactor], argument=[node]))
        allFactors = newFactors

    # Marginalize out the maturation nodes of the query lineage.
    productFactor = NoisyOrFactor("product", allFactors)
    assignments = assignmentHelper(model, directParents, directChildren, querySpecialEvidence, criticalMaturationChildren, criticalMaturation[0])
    reducedFactors = []
    for assignment in assignments.keys():
        if assignments[assignment] > 0:
            reducedFactors.append(NoisyOrFactor("product", [NoisyOrFactor("constant", [], argument=[assignments[assignment]]), NoisyOrFactor("reduce", [productFactor], argument=[(node, 0 if node in assignment else 1) for node in criticalMaturation])]))
    finalFactor = NoisyOrFactor("sum", reducedFactors)

    finalFactorValue = finalFactor.get_value()

    finalResult = finalFactorValue * zerosConstant * forcedOneTransmitFailureProbability

    if queryPointValues:
        queryPointValues[queryPoint] = finalResult
    return finalResult

def handler(signum, frame):
    raise TimeLimitException("queryPoint evaluation timeout")

"""
Evaluates every query point one at a time.
"""
def evaluateModel(modelFolder, dataFolder, modelName, modelExtension, modelExtensionExtra, timeout=10, save=True, debug=0, progressBar=False, loadedModel=None):
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
        print("Loading naive probabilities.")
    with open(modelFolder + modelName + "_modeldata" + modelExtension + modelExtensionExtra + "_naiveProbabilities.pickle", "rb") as f:
        naiveProbabilities = pickle.load(f)

    # Get set of query points to evaluate.
    queryPoints = set()
    for query in queries.keys():
        for queryPoint in queries[query]["critical"]:
            queryPoints.add(queryPoint)

    if timeout > 0:
        # Try to use signals to catch querypoints that take too long.
        signal.signal(signal.SIGALRM, handler)

    # Calculate each query point.
    queryPointValues = dict()
    if debug >= 1:
        print("Calculating each query point.")
    if progressBar:
        iterator = tqdm.tqdm(sorted(queryPoints))
    else:
        iterator = sorted(queryPoints)
    for queryPoint in iterator:
        if progressBar:
            iterator.set_description(desc="Starting query point " + str(queryPoint))
        if debug >= 2:
            print("Starting query point " + str(queryPoint))
        
        # Multiprocessing takes too long to pass the data to the subprocess.
        """
        # Start evaluation on second thread and kill if it takes longer than timeout.
        p = multiprocessing.Process(target=evaluateQueryPoint, args=(model, directParents, directChildren, evidence, queries, criticalSegments, naiveProbabilities, queryPoint, queryPointValues))

        p.start()

        if timeout > 0:
            p.join(timeout)
            if p.is_alive():
                p.terminate()
        else:
            p.join()
        """

        #queryPointValues[queryPoint] = evaluateQueryPoint(model, directParents, directChildren, evidence, queries, criticalSegments, naiveProbabilities, queryPoint, queryPointValues)

        if timeout > 0:
            signal.alarm(timeout)

        try:
            queryPointValues[queryPoint] = evaluateQueryPoint(model, directParents, directChildren, evidence, queries, criticalSegments, naiveProbabilities, queryPoint, queryPointValues)
            if timeout > 0:
                signal.alarm(0) # Cancel the previous alarm.
        except TimeLimitException:
            queryPointValues[queryPoint] = -1

    if debug >= 1:
        print("Finished calculating all query points!")

    if save:
        if debug >= 1:
            print("Saving query point values.")
        with open(modelFolder + modelName + "_modeldata" + modelExtension + modelExtensionExtra + "_querypoints.pickle", "wb") as f:
            pickle.dump(queryPointValues, f)

    return queryPointValues
