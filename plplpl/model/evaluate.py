import collections
import itertools
import math
import pickle

from pgmpy.factors import factor_product

from plplpl.NoisyOr import BinaryNoisyOrCPD
from plplpl.NoisyOr import NoisyOrBayesianNetwork
from plplpl.NoisyOr import NoisyOrFactor

"""
modelFolder: path to directory
dataFolder: path to directory
modelName: unique name for this model/computation
modelExtension: model extension in the form `_[conjugation function]_[colour function]_[maturation function]` (don't include _contradictionsPruned_normalized.pickle)
save: if we should save the model to a file (pickle)
debug: 0 = nothing, 1 = status, 2 = verbose
progressBar: if we should show a progress bar on long for loops
cacheObjects: if we should keep a reference to used objects while computing queries.
loadedModel: if we should use a model already in memory instead of loading one.
loadedEvidence: if we should use evidence already in memory instead of loading it.
loadedConjugateQueries: if we should use conjugate queries already in memory instead of loading them.
loadedNonConjugateQueries: if we should use nonconjugate queries already in memory instead of loading them.
loadedFullQueries: if we should use full queries already in memory instead of loading them.
loadedNaiveProbabilities: if we should use naive probabilities already in memory instead of loading them.
loadedIncomingProbabilities: if we should use incoming probabilities already in memory instead of loading them.

TODO: Detailed description of the evaluation calculated.
"""

# Helper function for recursively calculating assignments.
# Skips dead ends.
def assignmentHelper(model, node, endpoints, prefix):
    if node in endpoints:
        return [[]]

    assignments = []
    children = [child for child in model.successors(node) if child[0] == prefix]
    if len(children) == 0:
        return None
    subAssignments = []
    for child in children:
        subAssignment = assignmentHelper(model, child, endpoints, prefix)
        if subAssignment != None:
            subAssignments.append(subAssignment)
    if len(subAssignments) == 0:
        return None
    for subAssignment in itertools.product(*subAssignments):
        newAssignment = [node]
        for componentAssignment in subAssignment:
            for newNode in componentAssignment:
                newAssignment.append(newNode)
        assignments.append(newAssignment)
    if (len(assignments) > 0):
        assignments.append([])
    return assignments

# Helper function for recursively building assignments.
# Uses a dictionary of children of nodes instead of the model itself.
def assignmentHelperFromDictionary(children, node):
    assignments = [[]]
    if len(children[node]) == 0:
        return [[node], []]
    subAssignments = []
    for child in children[node]:
        subAssignments.append(assignmentHelperFromDictionary(children, child))
    for subAssignment in itertools.product(*subAssignments):
        newAssignment = [node]
        for componentAssignment in subAssignment:
            for newNode in componentAssignment:
                newAssignment.append(newNode)
        assignments.append(newAssignment)
    return assignments

# Helper function for chunk based approach to assignments.
def assignmentHelperByChunk(model, node, endpoints, prefix, includeDeadEnds=False):
    chunk = []
    while node not in endpoints:
        chunk.append(node)
        children = [child for child in model.successors(node) if child[0] == prefix]
        if len(children) == 1:
            node = children[0]
        elif len(children) == 0:
            if includeDeadEnds:
                return (tuple(chunk), [])
            else:
                return (tuple(), [])
        elif len(children) > 1:
            nextChunks = []
            for child in children:
                nextChunk = assignmentHelperByChunk(model, child, endpoints, prefix, includeDeadEnds=includeDeadEnds)
                if len(nextChunk[0]) == 0:
                    continue
                nextChunks.append(nextChunk)
            return (tuple(chunk), nextChunks)
    return (tuple(chunk), [])

# We marginalize out all nodes in the given chunk.
def chunkCalculator(model, geneAssignment, endpoints, variableDownwardZeroPrecomputed, allOneChunkProbabilities, downstreamFromChunk, upstreamFromChunk, childChunks, chunk, factors, preconditionProbability):
    # chunk = (tuple(nodes), [childChunks])
    # Assumes nodes are ordered correctly.
    
    # Get coefficient of all ones beneath this chunk.
    childChunkAllOneConstant = 1.0
    for child in chunk[1]:
        childChunkAllOneConstant = childChunkAllOneConstant * allOneChunkProbabilities[child[0]]

    # Get P(node = 1 | parents(node)) for node in chunk.
    probZeroAt = dict()
    remainingAssignment = {}
    for node in chunk[0]:
        remainingAssignment[node] = 0
        for parent in model.predecessors(node):
            if parent not in geneAssignment:
                remainingAssignment[parent] = 0
        probZeroAt[node] = model.get_cpds(node).get_value(**geneAssignment, **remainingAssignment)

    # Compute the probability of each assignment to chunk that includes maturation in this chunk.
    maturationAt = dict()
    runningDownwardZeroProbability = 1.0
    for node in chunk[0]:
        runningDownwardZeroProbability = runningDownwardZeroProbability * variableDownwardZeroPrecomputed[node]
    runningProbability = preconditionProbability
    marginalizeProb = dict()
    for node in chunk[0]:
        # P(maturationAt(node)) = prod(
        #   P(not mature yet) = runningProbability
        #   P(matured here) = 1 - probZeroAt[node]
        #   P(downwardZeros correct) = runningDownwardZeroProbability
        # )
        maturationAt[node] = runningProbability * (1 - probZeroAt[node]) * runningDownwardZeroProbability * childChunkAllOneConstant
        runningProbability = runningProbability * probZeroAt[node]
        marginalizeProb[node] = runningProbability
        runningDownwardZeroProbability = runningDownwardZeroProbability / variableDownwardZeroPrecomputed[node]

    # Combine all factors (only for maturation in this chunk)
    if len(factors) > 0:
        productFactor = NoisyOrFactor("product", factors)
    else:
        productFactor = NoisyOrFactor("constant", [], [1.0])
    toReduce = []
    if chunk[0] in downstreamFromChunk:
        for variable in productFactor.variables:
            if variable in downstreamFromChunk[chunk[0]]:
                toReduce.append((variable, 1))
    if len(toReduce) > 0:
        reducedProductFactor = NoisyOrFactor("reduce", [productFactor], argument=toReduce)
    else:
        reducedProductFactor = productFactor

    # Compute the factor for each maturation assignment.
    maturationAtFactors = []
    zeroSoFar = set()
    for node in chunk[0]:
        toReduce = []
        for variable in reducedProductFactor.variables:
            if variable in chunk[0]:
                if variable in zeroSoFar:
                    toReduce.append((variable, 0))
                else:
                    toReduce.append((variable, 1))
        if len(toReduce) > 0:
            tempFactor = NoisyOrFactor("reduce", [reducedProductFactor], argument=toReduce)
        else:
            tempFactor = reducedProductFactor
        constantFactor = NoisyOrFactor("constant", [], argument=[maturationAt[node]])
        maturationAtFactors.append(NoisyOrFactor("product", [constantFactor, tempFactor]))

    # Reduce all the factors to this chunk being zero.
    newFactors = []
    for factor in factors:
        toReduce = []
        for variable in factor.variables:
            if variable in chunk[0]:
                toReduce.append((variable, 0))
        if len(toReduce) > 0:
            newFactors.append(NoisyOrFactor("reduce", [factor], argument=toReduce))
        else:
            newFactors.append(factor)

    # Calculate recursive chunks.
    for nextChunk in chunk[1]:
        factorsToSend = []
        unused = []
        for factor in newFactors:
            for variable in factor.variables:
                if (variable in nextChunk[0]) or (variable in downstreamFromChunk[nextChunk[0]]):
                    factorsToSend.append(factor)
                    break
            else:
                unused.append(factor)
        recursiveFactor = chunkCalculator(model, geneAssignment, endpoints, variableDownwardZeroPrecomputed, allOneChunkProbabilities, downstreamFromChunk, upstreamFromChunk, childChunks, nextChunk, factorsToSend, runningProbability)
        if recursiveFactor:
            unused.append(recursiveFactor)
        newFactors = unused

    if len(newFactors) > 0:
        # Once all recursive calls are complete, multiply any remaining factors.
        nonMatureFactor = NoisyOrFactor("product", newFactors)

        # Final factor to return is sum of this and maturing in this chunk factors.
        maturationAtFactors.append(nonMatureFactor)

    if len(maturationAtFactors) == 0:
        return None
    return NoisyOrFactor("sum", maturationAtFactors)

# Helper function to construct blocks for connectedDown queries.
# Returns [blocks], [parents_of_blocks]
def blockBuilder(model, node, endpoints, source=None):
    block = []
    current = node
    while True:
        block.append(current)
        if current in endpoints:
            return [block], [None]
        children = [x for x in model.successors(current) if x[0] == "g"]
        if len(children) == 0:
            return [], []
        elif len(children) >= 2:
            blocks = []
            parents = []
            for child in children:
                newBlocks, newParents = blockBuilder(model, child, endpoints)
                blocks = blocks + newBlocks
                parents = parents + newParents
            if len(blocks) > 0:
                return [block] + blocks, [None] + [(current if x == None else x) for x in parents]
            else:
                return [], []
        else:
            current = children[0]

# Helper function to merge queries together recursively.
def queryMergeHelper(conQueries, evaluations, precons, query):
    if len(conQueries[query]["childQueries"]) > 0:
        children = []
        for child in conQueries[query]["childQueries"]:
            children.append(queryMergeHelper(conQueries, evaluations, precons, child))
        childChance = precons[query] * math.prod(children)
        #return (evaluations[query] + childChance - (evaluations[query] * childChance))
        # This not the above because it isn't an OR, it's the sum of probability of all query assignments that we accept.
        return evaluations[query] + childChance
    else:
        return evaluations[query]

def evaluateModel(modelFolder, dataFolder, modelName, modelExtension, save=True, debug=0, progressBar=False, cacheObjects=True, loadedModel=None, loadedEvidence=None, loadedConjugateQueries=None, loadedNonConjugateQueries=None, loadedFullQueries=None, loadedNaiveProbabilities=None, loadedIncomingProbabilities=None):

    # Set this to a specific query to skip directly to it and then enter debugger just before calculating it.
    manualDebug = None

    if progressBar:
        import tqdm

    if debug >= 1:
        print("Loading model.")
    if loadedModel:
        model = loadedModel
    else:
        with open(modelFolder + modelName + "_model" + modelExtension + "_contradictionsPruned_normalized.pickle", "rb") as f:
            model = pickle.load(f)

    previous = dict()
    for node in model.nodes:
        parent = None
        for predecessor in model.predecessors(node):
            if predecessor[0] == node[0]:
                parent = predecessor
        previous[node] = parent

    if debug >= 1:
        print("Loading cell names and unique ids.")
    with open(dataFolder + modelName + "_humanFriendlyNameLookup.pickle", "rb") as f:
        uid = pickle.load(f)
    with open(dataFolder + modelName + "_humanFriendlyName.pickle", "rb") as f:
        name = pickle.load(f)

    if debug >= 1:
        print("Loading evidence.")
    if loadedEvidence:
        evidence = loadedEvidence
    else:
        with open(modelFolder + modelName + "_modeldata" + modelExtension + "_evidence.pickle", "rb") as f:
            evidence = pickle.load(f)

    if debug >= 1:
        print("Loading queries.")
    if loadedConjugateQueries:
        conQueries = loadedConjugateQueries
    else:
        with open(modelFolder + modelName + "_modeldata" + modelExtension + "_completeConjugateQueries.pickle", "rb") as f:
            conQueries = pickle.load(f)
    if loadedNonConjugateQueries:
        ncQueries = loadedNonConjugateQueries
    else:
        with open(modelFolder + modelName + "_modeldata" + modelExtension + "_nonConjugateQueries.pickle", "rb") as f:
            ncQueries = pickle.load(f)
    if loadedFullQueries:
        fullQueries = loadedFullQueries
    else:
        with open(modelFolder + modelName + "_modeldata" + modelExtension + "_fullQueries.pickle", "rb") as f:
            fullQueries = pickle.load(f)

    if debug >= 1:
        print("Loading naive probabilities.")
    if loadedNaiveProbabilities:
        naiveProbabilities = loadedNaiveProbabilities
    else:
        with open(modelFolder + modelName + "_modeldata" + modelExtension + "_naiveProbabilities.pickle", "rb") as f:
            naiveProbabilities = pickle.load(f)

    if debug >= 1:
        print("Loading incoming probabilities.")
    if loadedIncomingProbabilities:
        incomingProbabilities = loadedIncomingProbabilities
    else:
        with open(modelFolder + modelName + "_modeldata" + modelExtension + "_incomingProbabilities.pickle", "rb") as f:
            incomingProbabilities = pickle.load(f)

    #
    # Conjugate Queries:
    # 

    conQueryEvaluation = dict()
    conQueryPreconditionConstants = dict()
    cachedObjects = list()

    if debug >= 1:
        print("Starting evaluation of conjugate queries.")

    if progressBar:
        #iterator = tqdm.tqdm(sorted(conQueries.keys(), key=lambda x: int(x.split("_")[1]), reverse=True))
        #iterator = tqdm.tqdm(sorted(conQueries.keys(), key=lambda x: int(x.split("_")[1])))
        iterator = tqdm.tqdm(conQueries.keys())
    else:
        iterator = conQueries.keys()

    for query in iterator:
        if manualDebug:
            if query != manualDebug:
                continue
        #NoisyOrFactor._registry = dict()
        if debug >= 2:
            print(query + " - Initializing")
        if progressBar:
            iterator.set_description(desc=query + " - Initializing")

        if debug >= 2:
            print(query + " - calculating greater query and implications")
        if progressBar:
            iterator.set_description(query + " - calculating greater query and implications")

        # Discern which full query we are a part of.
        fullQuery = None
        for q in fullQueries.keys():
            if query in fullQueries[q]:
                fullQuery = q

        # Get other queries upward.
        parentQueries = []
        if len(conQueries[query]["parentQueries"]) > 0:
            parent = conQueries[query]["parentQueries"][0]
            parentQueries.append(parent)
            while len(conQueries[parent]["parentQueries"]) > 0:
                parent = conQueries[parent]["parentQueries"][0]
                parentQueries.append(parent)

        # Get other queries downward.
        childQueries = []
        if len(conQueries[query]["childQueries"]) > 0:
            children = conQueries[query]["childQueries"]
            while children:
                newChildren = []
                for child in children:
                    childQueries.append(child)
                    newChildren = newChildren + conQueries[child]["childQueries"]
                children = newChildren

        # The rest are siblings.
        siblingQueries = []
        for q in fullQueries[fullQuery]:
            if (q == query) or (q in childQueries) or (q in parentQueries):
                continue
            siblingQueries.append(q)

        # Get evidence implied by the current query calculation.
        pseudoEvidence = dict()
        # Parents are implicitly zero.
        for parentQuery in parentQueries:
            for node in (conQueries[parentQuery]["query"]+conQueries[parentQuery]["critical"]):
                pseudoEvidence[node] = 0

        # Siblings are fixed, but have to be calculated.
        siblingNodes = []
        for siblingQuery in siblingQueries:
            for node in (conQueries[siblingQuery]["query"]+conQueries[siblingQuery]["critical"]):
                siblingNodes.append(node)
        for node in sorted(siblingNodes, key=lambda x: int(x.split("_")[1])):
            parent = previous[node] 
            pseudoEvidence[node] = incomingProbabilities[node] + pseudoEvidence[parent] - (incomingProbabilities[node] * pseudoEvidence[parent])

        # Clean up children coming off that die out.
        queue = collections.deque()

        for node in pseudoEvidence.keys():
            if node[0] == "g":
                queue.append(node)

        while queue:
            node = queue.popleft()
            for child in model.successors(node):
                if (child[0] == "g") and (child not in pseudoEvidence) and (child not in conQueries[query]["critical"]) and (child not in conQueries[query]["query"]) and (child not in evidence):
                    pseudoEvidence[child] = incomingProbabilities[child] + pseudoEvidence[node] - (incomingProbabilities[child] * pseudoEvidence[node])
                    queue.append(child)

        # Calculate pseudoevidence for maturation nodes coming off of siblings etc.
        for node in sorted(conQueries[fullQuery]["hidden"], key=lambda x: int(x.split("_")[1])):
            if node in conQueries[query]["hidden"]:
                continue
            if (node in evidence) and (previous[node] not in fullQueries[fullQuery]):
                continue
            arguments = dict()
            for parent in model.predecessors(node):
                if parent in pseudoEvidence:
                    arguments[parent] = pseudoEvidence[parent]
                else:
                    # There was an exception here when this was evidence[parent].
                    # I think it's a dead end child acting weird, might be worth examining in the future..?
                    arguments[parent] = naiveProbabilities[parent]
            pseudoEvidence[node] = model.get_cpds(node).get_self_values_with_uncertain_evidence(**arguments)[1][0]

        # Only using precondition constant in query merging later.
        preconditionConstant = 1.0
        for node in (conQueries[query]["query"] + conQueries[query]["critical"]):
            cpd = model.get_cpds(node)
            for i in range(len(cpd.evidence)):
                if cpd.evidence[i][0] == "g":
                    continue
                elif cpd.evidence[i][0] == "m":
                    if cpd.evidence[i] in pseudoEvidence:
                        preconditionConstant = preconditionConstant * (1 - (pseudoEvidence[cpd.evidence[i]] * cpd.evidence_noise[i]))
                    else:
                        preconditionConstant = preconditionConstant * (1 - (naiveProbabilities[cpd.evidence[i]] * cpd.evidence_noise[i]))
                else:
                    print("Something weird happened. Shouldn't get here.")
        
        # Find all the possible assignments to query and critical nodes.
        firstNode = sorted(conQueries[query]["query"] + conQueries[query]["critical"], key=lambda x: int(x.split("_")[1]))[0]
        QCzeroAssignments = assignmentHelper(model, firstNode, conQueries[query]["query"], "g")

        if progressBar:
            iterator.set_description(desc=query + " - Constructing simplified downward queries")

        # Keep track of everything that is incoming to this query.
        allIncoming = set(conQueries[query]["incoming"])
        allIncomingQueries = set()
        for incoming in allIncoming:
            for q in fullQueries.keys():
                if incoming in conQueries[q]["hidden"]:
                    allIncomingQueries.add(q)

        # Build "connectedDown" query blocks.
        # We only care about blocks where this query directly interacts with them.
        downwardQueryCPDs = list()
        downwardQueryEvidence = list()
        downwardQueryConstant = 1.0
        seenConnectedDown = set()
        #seenIncomingQueries = set()
        for downwardQuery in sorted(conQueries[query]["connectedDown"].keys(), key=lambda x: int(x.split("_")[1])):
            if downwardQuery in seenConnectedDown:
                continue
            seenConnectedDown.add(downwardQuery)
            if downwardQuery in fullQueries[fullQuery]:
                continue # This is pointing downward into another arm of itself, we don't have to handle it here.
            fixed = True
            for hidden in conQueries[query]["connectedDown"][downwardQuery]:
                # Check if this downward query is not connected to us directly.
                for parent in model.predecessors(hidden):
                    if parent[0] != "g":
                        continue
                    if parent not in conQueries[query]["lineage"]:
                        fixed = False
                        break
                if fixed == False:
                    break

            # If so, it resolves to trivial.
            # This is because we will divide by probability given the best assignment to this query, which is always true.
            if fixed:
                continue

            # Setup the block builder.
            endpoints = []
            # Going backwards and starting at the parent is unnecessary, since we include those in the precondition or already see them.
            """
            parent = downwardQuery
            while len(conQueries[parent]["parentQueries"]) > 0:
                parent = conQueries[parent]["parentQueries"][0]
            queue = collections.deque()
            queue.append(parent)
            firstBlockNode = sorted(conQueries[parent]["query"] + conQueries[parent]["critical"], key=lambda x: int(x.split("_")[1]))[0]
            """
            queue.append(downwardQuery)
            firstBlockNode = sorted(conQueries[downwardQuery]["query"] + conQueries[downwardQuery]["critical"], key=lambda x: int(x.split("_")[1]))[0]
            while queue:
                q = queue.popleft()
                seenConnectedDown.add(q)
                if len(conQueries[q]["childQueries"]) > 0:
                    for child in conQueries[q]["childQueries"]:
                        queue.append(child)
                else:
                    for end in conQueries[q]["query"]:
                        endpoints.append(end)
            blocks, parents = blockBuilder(model, firstBlockNode, endpoints)
            blockCPDs = []
            for i in range(len(blocks)):
                # Build a new BinaryNoisyOrCPD of all possible causes of each block.
                variable = blocks[i][-1]
                internalProbability = 0
                cpd_evidence = []
                cpd_evidence_noise = []
                if parents[i]:
                    cpd_evidence.append(parents[i])
                    cpd_evidence_noise.append(1.0)
                for node in blocks[i]:
                    cpd = model.get_cpds(node)
                    for parent in model.predecessors(node):
                        if parent[0] == "m":
                            cpd_evidence.append(parent)
                            cpd_evidence_noise.append(cpd.evidence_noise[cpd.evidence.index(parent)])
                blockCPD = BinaryNoisyOrCPD(variable, internalProbability, evidence=cpd_evidence, evidence_noise=cpd_evidence_noise)
                values = dict()
                for variable in blockCPD.evidence:
                    if variable[0] == "g":
                        continue
                    # If it's part of a sibling, child, or parent query, we immediately add it to be reduced out.
                    # This is so that we only consider the difference this query makes to the result.
                    if variable in pseudoEvidence:
                        values[variable] = pseudoEvidence[variable]
                        continue
                    # If it's relevant to us, continue, we'll marginalize it out later.
                    if variable in conQueries[query]["hidden"]:
                        continue
                    if variable in allIncoming:
                        continue
                    skip = False
                    # Don't add other queries if they aren't directly connected to the current ones.
                    # It's significantly slower for a correlation occurring two queries away.
                    # But also only include it if it isn't in *their* fixed lineage.
                    if variable in evidence:
                        values[variable] = evidence[variable]
                        continue
                    #for q in fullQueries:
                    for q in allIncomingQueries:
                        if variable in conQueries[q]["hidden"]:
                            found = False
                            for incoming in allIncoming:
                                if incoming not in conQueries[q]["hidden"]:
                                    continue
                                predecessor = incoming
                                while (previous[predecessor] not in conQueries[query]["incoming"]) and (previous[predecessor] not in evidence):
                                    if previous[predecessor] == variable:
                                        # We only care if it's upstream from a direct incoming node.
                                        found = True
                                        break
                                    else:
                                        predecessor = previous[predecessor]
                                if found == True:
                                    break
                            if found:
                                allIncoming.add(variable)
                                #seenIncomingQueries.add(q)
                                skip = True
                                break
                            else:
                                # Otherwise, we can just use naive probabilities as an estimate and reduce it out.
                                pass
                            #allIncomingQueries.add(q)
                            #seenIncomingQueries.add(q)
                    if skip:
                        continue
                    
                    # Otherwise, reduce it.
                    values[variable] = naiveProbabilities[variable]

                blockCPD.reduce_with_uncertain_evidence([(key, value) for key, value in values.items()], inplace=True)
                blockCPDs.append(blockCPD)
            downwardQueryCPDs.append(blockCPDs)
            downwardQueryEvidence.append(endpoints.copy())

        if progressBar:
            iterator.set_description(desc=query + " - Computing downwardQueryConditionConstant")

        downwardQueryConditionConstant = 1.0

        downwardQueryFactors = []
        downwardQueryToMarginalize = []
        for i in range(len(downwardQueryCPDs)):
            # Figure out what the maximum value of these blocks are alongside the cleanup.
            productFactorParts = []
            toMarg = []
            nextDownwardQueryFactors = []
            nextDownwardQueryToMarginalize = []
            for cpd in downwardQueryCPDs[i]:
                if cpd.variable in downwardQueryEvidence[i]:
                    reducedFactor = NoisyOrFactor("reduce", [cpd], argument=[(cpd.variable, 1)])
                    nextDownwardQueryFactors.append(reducedFactor)
                    productFactorParts.append(reducedFactor)
                else:
                    nextDownwardQueryToMarginalize.append(cpd.variable)
                    nextDownwardQueryFactors.append(cpd)
                    toMarg.append(cpd.variable)
                    productFactorParts.append(cpd)
            productFactor = NoisyOrFactor("product", productFactorParts)
            if len(toMarg) > 0:
                # Future work: efficiently do this marginalization.
                productFactor = NoisyOrFactor("marginalize", [productFactor], argument=toMarg)
            lastAppearance = True
            for variable in productFactor.variables:
                for childQuery in childQueries:
                    if variable in conQueries[childQuery]["hidden"]:
                        lastAppearance = False
                        break
                if lastAppearance == False:
                    break
            toReduce = []
            for v in productFactor.variables:
                if v in conQueries[query]["hidden"]:
                    continue
                if v in evidence:
                    toReduce.append((v, evidence[v]))
                else:
                    # This is an estimate, but it shouldn't make a significant difference as long as it's consistent.
                    # It's just significantly faster.
                    # Future work: find a refined estimate/quick calculation.
                    toReduce.append((v, 1))
            if len(toReduce) > 1:
                productFactor = NoisyOrFactor("reduce", [productFactor], argument=toReduce)
            """
            if len([v for v in productFactor.variables if (v not in conQueries[query]["hidden"]) and (v in evidence)]) > 0:
                productFactor = NoisyOrFactor("reduce", [productFactor], argument=[(v, evidence[v]) for v in productFactor.variables if (v not in conQueries[query]["hidden"]) and (v in evidence)])
            if len([v for v in productFactor.variables if v not in conQueries[query]["hidden"]]) > 0:
                productFactor = NoisyOrFactor("marginalizeGivenProbability", [productFactor], argument=[(v, naiveProbabilities[v]) for v in productFactor.variables if v not in conQueries[query]["hidden"]])
            """
            value = NoisyOrFactor("reduce", [productFactor], argument=[(v, 1) for v in productFactor.variables]).get_value()
            if value == 0:
                print("A downwardQuery " + " can't be satisfied, it's being skipped.")
            else:
                if lastAppearance:
                    zeroValue = NoisyOrFactor("reduce", [productFactor], argument=[(v, 0) for v in productFactor.variables]).get_value()
                    preconditionConstant = preconditionConstant * (zeroValue / value)

                downwardQueryConditionConstant = downwardQueryConditionConstant * value
                for factor in nextDownwardQueryFactors:
                    downwardQueryFactors.append(factor)
                for v in nextDownwardQueryToMarginalize:
                    downwardQueryToMarginalize.append(v)

        if progressBar:
            iterator.set_description(desc=query + " - Computing shared downward zeroes")

        # Calculate the downwardZero values that are shared for all assignments.
        # Since the downward zero stuff only matters for non-zero, we only have to consider this query and child queries.
        variableDownwardZero = set()
        for downwardZero in conQueries[query]["downwardZero"]:
            for parent in model.predecessors(downwardZero):
                if parent in conQueries[query]["hidden"]:
                    for grandParent in model.predecessors(parent):
                        if grandParent[0] == "g" and grandParent not in conQueries[query]["lineage"]:
                            variableDownwardZero.add(downwardZero)
                            break
                if downwardZero in variableDownwardZero:
                    break
        fixedDownwardZeroMultiplier = 1.0
        for downwardZero in conQueries[query]["downwardZero"]:
            if downwardZero in variableDownwardZero:
                continue
            cpd = model.get_cpds(downwardZero)
            for parent in model.predecessors(downwardZero):
                if parent in conQueries[query]["hidden"]:
                    fixedDownwardZeroMultiplier = fixedDownwardZeroMultiplier * (1.0 - cpd.evidence_noise[cpd.evidence.index(parent)])

        # Precompute all variable downward zero values to be used as needed later.
        variableDownwardZeroPrecomputed = dict()
        for hidden in conQueries[query]["hidden"]:
            variableDownwardZeroPrecomputed[hidden] = 1.0
        for downwardZero in variableDownwardZero:
            cpd = model.get_cpds(downwardZero)
            for parent in model.predecessors(downwardZero):
                if parent in conQueries[query]["hidden"]:
                    variableDownwardZeroPrecomputed[hidden] = variableDownwardZeroPrecomputed[hidden] * (1.0 - cpd.evidence_noise[cpd.evidence.index(parent)])

        if progressBar:
            iterator.set_description(desc=query + " - Examining incoming queries") 

        # Figure out which queries are incoming
        incomingQueries = dict() # Dictionary with a list of relevant variables, as opposed to set allIncomingQueries
        querylessIncoming = set()
        for incoming in allIncoming:
            for incomingQuery in fullQueries.keys():
                if incoming in conQueries[incomingQuery]["hidden"]:
                    if incomingQuery not in incomingQueries:
                        incomingQueries[incomingQuery] = set()
                    incomingQueries[incomingQuery].add(incoming)
                    break
            else:
                querylessIncoming.add(incoming)
        # Find the structure among each set of incoming queries.
        incomingParents = dict()
        incomingChildren = dict()
        # Get structure for everything, then use this to build "queries" for the queryless incoming nodes.
        for node in allIncoming:
            incomingChildren[node] = set()
            children = collections.deque(model.successors(node))
            while children:
                child = children.popleft()
                if child[0] != "m":
                    continue
                if child in evidence:
                    continue
                if child in allIncoming:
                    incomingParents[child] = node
                    incomingChildren[node].add(child)
                    continue
                for grandChild in model.successors(child):
                    if grandChild[0] != "m":
                        continue
                    if grandChild in evidence:
                        continue
                    children.append(grandChild)
        
        # Add a pseudoquery for all the hidden nodes connected to each other but aren't in a query.
        # Just so they are appropriately computed together.
        done = set()
        trivialIncomingMaturation = set()
        for incoming in querylessIncoming:
            if incoming in done:
                continue
            done.add(incoming)
            pseudoquery = set()
            pseudoquery.add(incoming)
            queue = collections.deque()
            for child in incomingChildren[incoming]:
                if child not in done:
                    queue.append(child)
            if incoming in incomingParents:
                parent = incomingParents[incoming]
                if parent not in done:
                    queue.append(parent)
            while queue:
                node = queue.popleft()
                if node in done:
                    continue
                done.add(node)
                pseudoquery.add(node)
                for child in incomingChildren[node]:
                    if child not in done:
                        queue.append(child)
                if node in incomingParents:
                    parent = incomingParents[node]
                    if parent not in done:
                        queue.append(parent)
            # If this is a standalone node, check if it appears in any incoming query.
            # Otherwise, we can reduce by it in the critical and switchpoint cpds.
            if len(pseudoquery) == 1:
                for factor in downwardQueryFactors:
                    if incoming in factor.variables:
                        incomingQueries[incoming] = pseudoquery
                        break
                else:
                    # It can be reduced into the critical/query cpds.
                    trivialIncomingMaturation.add(incoming)
                    allIncoming.remove(incoming)
                    if incoming in incomingParents:
                        incomingParents.pop(incoming)
                    if incoming in incomingChildren:
                        incomingChildren.pop(incoming)
            else:
                incomingQueries[incoming] = pseudoquery
        # All incoming nodes should be in incomingQueries now.

        if progressBar:
            iterator.set_description(desc=query + " - Building simplified CPDs for incoming nodes")

        # Build simplified CPDs for incoming hidden nodes.
        incomingCPDs = []
        if progressBar:
            nestedIterator = tqdm.tqdm(sorted(allIncoming, key=lambda x: int(x.split("_")[1])), leave=False)
        else:
            nestedIterator = sorted(allIncoming, key=lambda x: int(x.split("_")[1]))
        for incoming in nestedIterator:
            if incoming in pseudoEvidence:
                if progressBar:
                    nestedIterator.set_description(desc="Working on " + incoming + " (pseudo)")
                if incoming in incomingParents:
                    arguments = dict()
                    for predecessor in model.predecessors(incoming):
                        if predecessor[0] == "m":
                            arguments[predecessor] = 0
                        elif predecessor in pseudoEvidence:
                            arguments[predecessor] = pseudoEvidence[predecessor]
                        else:
                            arguments[predecessor] = evidence[predecessor]
                    incomingProbability = model.get_cpds(incoming).get_self_values_with_uncertain_evidence(**arguments)[1][0]
                    parent = previous[incoming]
                    while (parent != incomingParents[incoming]):
                        arguments = dict()
                        for predecessor in model.predecessors(parent):
                            if predecessor[0] == "m":
                                arguments[predecessor] = 0
                            elif predecessor in pseudoEvidence:
                                arguments[predecessor] = pseudoEvidence[predecessor]
                            else:
                                arguments[predecessor] = evidence[predecessor]
                        parentIncomingProbability = model.get_cpds(parent).get_self_values_with_uncertain_evidence(**arguments)[1][0]
                        incomingProbability = incomingProbability + parentIncomingProbability - (incomingProbability * parentIncomingProbability)
                        parent = previous[parent]
                    incomingCPDs.append(BinaryNoisyOrCPD(incoming, incomingProbability, evidence=[incomingParents[incoming]], evidence_noise=[1.0]))
                else:
                    incomingCPDs.append(BinaryNoisyOrCPD(incoming, pseudoEvidence[incoming]))
            else:
                if progressBar:
                    nestedIterator.set_description(desc="Working on " + incoming)
                if incoming in incomingParents:
                    incomingProbability = incomingProbabilities[incoming]
                    parent = None
                    for predecessor in model.predecessors(incoming):
                        if predecessor[0] == "m":
                            parent = predecessor
                    while (parent != incomingParents[incoming]):
                        incomingProbability = incomingProbability + incomingProbabilities[parent] - (incomingProbability * incomingProbabilities[parent])
                        for predecessor in model.predecessors(parent):
                            if predecessor[0] == "m":
                                parent = predecessor
                    incomingCPDs.append(BinaryNoisyOrCPD(incoming, incomingProbability, evidence=[incomingParents[incoming]], evidence_noise=[1.0]))
                else:
                    incomingCPDs.append(BinaryNoisyOrCPD(incoming, naiveProbabilities[incoming]))

        if progressBar:
            iterator.set_description(desc=query + " - Reducing out lineage nodes")

        # Calculate the parents of every hidden node.
        # Also decide which hidden nodes are always fixed due to only having parents in lineage.
        hiddenParents = dict()
        lineageHidden = set()
        for hidden in conQueries[query]["hidden"]:
            allLineage = True
            for parent in model.predecessors(hidden):
                if parent[0] == "m":
                    if parent in conQueries[query]["hidden"]:
                        hiddenParents[hidden] = parent
                    else:
                        hiddenParents[hidden] = None
                elif parent[0] == "g":
                    if (parent in conQueries[query]["query"]) or (parent in conQueries[query]["critical"]):
                        allLineage = False
            if allLineage:
                lineageHidden.add(hidden)

        # Reduce by the lineage hidden nodes, since they are the same across all assignments.
        newFactors = []
        for factor in downwardQueryFactors:
            toReduce = []
            for variable in factor.variables:
                if variable in lineageHidden:
                    toReduce.append((variable, 1))
            if len(toReduce) > 0:
                newFactors.append(NoisyOrFactor("reduce", [factor], toReduce))
            else:
                newFactors.append(factor)
        downwardQueryFactors = newFactors

        # For each assignment, build the query to evaluate.
        if progressBar:
            iterator.set_description(desc=query + " - Evaluating each assignment to gene nodes")
            nestedIterator = tqdm.tqdm([set(x) for x in QCzeroAssignments], leave=False)
        else:
            nestedIterator = [set(x) for x in QCzeroAssignments]

        finishedFactors = []

        totalAssignmentProbability = 0.0
        for assignment in nestedIterator:
            if progressBar:
                nestedIterator.set_description(desc="Calculating probability of colour evidence    ")
            # Calculate the probability that the colours stayed zero until they changed.
            colourEvidenceProbability = 1.0
            for node in conQueries[query]["colourEvidence"]:
                cpd = model.get_cpds(node)
                colourAssignment = dict()
                colourAssignment[node] = 0
                for parent in model.predecessors(node):
                    if parent[0] == "c":
                        colourAssignment[parent] = 0
                    elif parent in assignment:
                        colourAssignment[parent] = 0
                    elif parent in conQueries[query]["query"]:
                        colourAssignment[parent] = 1
                    elif parent in conQueries[query]["critical"]:
                        colourAssignment[parent] = 1
                    else:
                        colourAssignment[parent] = 0
                colourEvidenceProbability = colourEvidenceProbability * cpd.get_value(**colourAssignment)

            if progressBar:
                nestedIterator.set_description(desc="Calculating reduced cpds for critical region  ")

            criticalNodeReducedCPDs = list()
            switchPointCPDs = list()

            maturationStartpoints = set()
            maturationEndpoints = set()
            # Get reduced cpd for each query/critical node set to zero.
            for node in assignment:
                hardEvidence = dict()
                for parent in model.predecessors(node):
                    if parent[0] == "g":
                        hardEvidence[parent] = 0
                    elif parent in conQueries[query]["hardEvidence"]:
                        hardEvidence[parent] = evidence[parent]
                    elif parent in trivialIncomingMaturation:
                        hardEvidence[parent] = naiveProbabilities[parent]
                criticalNodeReducedCPDs.append(model.get_cpds(node).reduce_with_uncertain_evidence([(key, value) for key, value in hardEvidence.items()], inplace=False))
                # Check if this is the last in the assignment, if so, also do the next node.
                # This is the point where it's one.
                for child in model.successors(node):
                    if (child[0] == "g") and ((child in conQueries[query]["critical"]) or (child in conQueries[query]["query"])) and (child not in assignment):
                        hardEvidence = dict()
                        for parent in model.predecessors(child):
                            if parent[0] == "g":
                                hardEvidence[parent] = 0
                            elif parent in conQueries[query]["hardEvidence"]:
                                hardEvidence[parent] = evidence[parent]
                            elif parent in trivialIncomingMaturation:
                                hardEvidence[parent] = naiveProbabilities[parent]
                        switchPointCPDs.append(model.get_cpds(child).reduce_with_uncertain_evidence([(key, value) for key, value in hardEvidence.items()], inplace=False))
                        # Also locate maturation assignments from here.
                        for grandChild in model.successors(child):
                            if grandChild[0] != "m":
                                continue
                            if int(grandChild.split("_")[1])-int(child.split("_")[1]) == model.constants["maturation_min"]:
                                maturationStartpoints.add(grandChild)
                            if int(grandChild.split("_")[1])-int(child.split("_")[1]) == model.constants["maturation_max"]:
                                maturationEndpoints.add(grandChild)

            # Special case for all ones.
            if len(assignment) == 0:
                for node in conQueries[query]["critical"] + conQueries[query]["query"]:
                    hasParent = False
                    for parent in model.predecessors(node):
                        if (parent in conQueries[query]["critical"]) or (parent in conQueries[query]["query"]):
                            hasParent = True
                            break
                    if hasParent == False:
                        hardEvidence = dict()
                        for parent in model.predecessors(node):
                            if parent[0] == "g":
                                hardEvidence[parent] = 0
                            elif parent in conQueries[query]["hardEvidence"]:
                                hardEvidence[parent] = evidence[parent]
                            elif parent in trivialIncomingMaturation:
                                hardEvidence[parent] = naiveProbabilities[parent]
                        switchPointCPDs.append(model.get_cpds(node).reduce_with_uncertain_evidence([(key, value) for key, value in hardEvidence.items()], inplace=False))
                        for grandChild in model.successors(node):
                            if grandChild[0] != "m":
                                continue
                            if int(grandChild.split("_")[1])-int(node.split("_")[1]) == model.constants["maturation_min"]:
                                maturationStartpoints.add(grandChild)
                            if int(grandChild.split("_")[1])-int(node.split("_")[1]) == model.constants["maturation_max"]:
                                maturationEndpoints.add(grandChild)

            if progressBar:
                nestedIterator.set_description(desc="Calculating downward zero nodes               ")

            # Get all fixed zero hidden nodes.
            fixedMaturationZeroes = set()
            queue = collections.deque(maturationStartpoints)
            while queue:
                node = queue.popleft()
                for parent in model.predecessors(node):
                    if (parent in conQueries[query]["hidden"]) and (parent not in fixedMaturationZeroes):
                        fixedMaturationZeroes.add(parent)
                        queue.append(parent)

            if progressBar:
                nestedIterator.set_description(desc="Calculating gene node factors                 ")
            # Setup the factors dependent on the assignment to gene nodes.
            fixedFactors = list()
            for cpd in switchPointCPDs:
                #for parent in cpd.evidence:
                #    fixedToMarginalize.add(parent)
                #factor = cpd.to_factor()
                #factor.reduce([(cpd.variable, 1)], inplace=True)
                factor = NoisyOrFactor("reduce", [cpd], argument=[(cpd.variable, 1)])
                fixedFactors.append(factor)

            # Add in the incoming factors
            for cpd in incomingCPDs:
                fixedFactors.append(cpd)

            # Add in the downward query factors.
            for factor in downwardQueryFactors:
                fixedFactors.append(factor)
           
            if progressBar:
                nestedIterator.set_description(desc="Calculating maturation chunks                 ")

            maturationChunks = [assignmentHelperByChunk(model, x, maturationEndpoints, "m", includeDeadEnds=True) for x in maturationStartpoints]

            if progressBar:
                nestedIterator.set_description(desc="Calculating downward zero nodes for each chunk")

            # Precompute the probability of all ones downstream from some chunk.
            allOneChunkProbability = dict()
            allChunks = list()
            # It's actually a stack now, meh.
            queue = collections.deque(maturationChunks)
            while queue:
                chunk = queue.pop()
                nextChunks = list()
                for nextChunk in chunk[1]:
                    if nextChunk[0] not in allOneChunkProbability:
                        nextChunks.append(nextChunk)
                if len(nextChunks) > 0:
                    queue.append(chunk)
                    for nextChunk in nextChunks:
                        queue.append(nextChunk)
                else:
                    prob = 1.0
                    for node in chunk[0]:
                        prob = prob * variableDownwardZeroPrecomputed[node]
                    for nextChunk in chunk[1]:
                        prob = prob * allOneChunkProbability[nextChunk[0]]
                    allOneChunkProbability[chunk[0]] = prob
                    allChunks.append(chunk)

            if progressBar:
                nestedIterator.set_description(desc="Calculating relationships between chunks      ")

            # Calculate all nodes upstream and downstream from each chunk
            upstreamFromChunk = dict()
            downstreamFromChunk = dict()
            childChunks = dict()
            allChunkNodes = set()
            for chunk in allChunks:
                for node in chunk[0]:
                    allChunkNodes.add(node)
                downstreamFromChunk[chunk[0]] = set()
                childChunks[chunk[0]] = set()
                children = collections.deque(chunk[1])
                while children:
                    child = children.popleft()
                    for node in child[0]:
                        downstreamFromChunk[chunk[0]].add(node)
                    childChunks[chunk[0]].add(child[0])
                    if child[0] not in upstreamFromChunk:
                        upstreamFromChunk[child[0]] = set()
                    for node in chunk[0]:
                        upstreamFromChunk[child[0]].add(node)
                    for grandChild in child[1]:
                        children.append(grandChild)
 
            if progressBar:
                nestedIterator.set_description(desc="Reducing out all forced values                ")

            # Reduce factors by forced values for the assignment
            newFactors = []
            forcedHiddenNodes = dict()
            for factor in fixedFactors:
                toReduce = list()
                for variable in factor.variables:
                    if variable in allChunkNodes:
                        continue
                    if variable in lineageHidden:
                        toReduce.append((variable, 1))
                    elif variable in conQueries[query]["hidden"]:
                        if variable not in forcedHiddenNodes:
                            seen = set()
                            node = variable
                            while node:
                                seen.add(node)
                                if node in maturationEndpoints:
                                    for s in seen:
                                        forcedHiddenNodes[s] = 1
                                    break
                                node = hiddenParents[node]
                            if node == None:
                                for s in seen:
                                    forcedHiddenNodes[s] = 0
                        toReduce.append((variable, forcedHiddenNodes[variable]))
                    elif variable in evidence:
                        toReduce.append((variable, evidence[variable]))
                    elif (variable[0] == "g") and (variable in pseudoEvidence): # Catch weirdness around links.
                        toReduce.append((variable, pseudoEvidence[variable]))
                if len(toReduce) > 0:
                    newFactors.append(NoisyOrFactor("reduce", [factor], argument=toReduce))
                else:
                    newFactors.append(factor)
            fixedFactors = newFactors
            
            if progressBar:
                nestedIterator.set_description(desc="Marginalizing out incoming nodes              ")

            # Marginalize out all the incoming nodes from a query.
            # New approach: do so by reducing to each assignment and summing together.
            done = set()
            for incomingQuery in incomingQueries.keys():
                for root in [node for node in incomingQueries[incomingQuery] if node not in incomingParents]:
                    workingFactors = []
                    newFactors = []
                    incomingTree = set()
                    queue = collections.deque()
                    queue.append(root)
                    while queue:
                        v = queue.popleft()
                        incomingTree.add(v)
                        if v in incomingChildren:
                            for child in incomingChildren[v]:
                                queue.append(child)

                    for factor in fixedFactors:
                        for v in factor.variables:
                            if v in incomingTree:
                                workingFactors.append(factor)
                                break
                        else:
                            newFactors.append(factor)
                    if len(workingFactors) == 0:
                        continue
                    elif len(workingFactors) == 1:
                        productFactor = workingFactors[0]
                    else:
                        productFactor = NoisyOrFactor("product", workingFactors)
                    reducedFactors = []

                    for m_assignment in [[x] for x in assignmentHelperFromDictionary(incomingChildren, root)]:
                        zeroVariables = set()
                        for subAssignment in m_assignment:
                            for variable in subAssignment:
                                zeroVariables.add(variable)
                        constant = 1.0
                        # Add in constants from critical nodes set to zero, here.
                        for critical in assignment:
                            cpd = model.get_cpds(critical)
                            for i in range(len(cpd.evidence)):
                                if (cpd.evidence[i] in incomingTree) and (cpd.evidence[i] not in zeroVariables):
                                    constant = constant * (1 - cpd.evidence_noise[i])
                        if constant != 1:
                            reducedFactors.append(NoisyOrFactor("reduce", [NoisyOrFactor("product", [NoisyOrFactor("constant", [], argument=[constant]), productFactor])], argument=[(var, (0 if var in zeroVariables else 1)) for var in incomingTree]))
                        else:
                            reducedFactors.append(NoisyOrFactor("reduce", [productFactor], argument=[(var, (0 if var in zeroVariables else 1)) for var in incomingTree]))
                    newFactors.append(NoisyOrFactor("sum", reducedFactors))
                    fixedFactors = newFactors
 
            if progressBar:
                nestedIterator.set_description(desc="Marginalizing out intermediate nodes          ")

            # Marginalize out all the intermediate nodes in the downward queries
            for toMarginalize in downwardQueryToMarginalize:
                newFactors = []
                workingFactors = []
                for factor in fixedFactors:
                    if toMarginalize in factor.variables:
                        workingFactors.append(factor)
                    else:
                        newFactors.append(factor)
                if len(workingFactors) > 0:
                    newFactors.append((NoisyOrFactor("marginalize", [NoisyOrFactor("product", workingFactors)], argument=[toMarginalize])))
                fixedFactors = newFactors

            if progressBar:
                nestedIterator.set_description(desc="Calculating factor for all chunks             ")

            # Get assignment to all possible parents of maturation nodes.
            # Anything unspecified is zero.
            geneAssignment = {x:(0 if (x in assignment) or (x in conQueries[query]["lineage"]) else 1) for x in (conQueries[query]["query"]+conQueries[query]["critical"]+conQueries[query]["lineage"])}

            # Get the factor representing all assignments to maturation nodes.
            completeAssignmentFactor = chunkCalculator(model, geneAssignment, maturationEndpoints, variableDownwardZeroPrecomputed, allOneChunkProbability, downstreamFromChunk, upstreamFromChunk, childChunks, (tuple(), maturationChunks), fixedFactors, 1.0)
 
            finalFactor = completeAssignmentFactor

            finishedFactors.append(finalFactor)

            if cacheObjects:
                cachedObjects.append(finalFactor)

            if progressBar:
                nestedIterator.set_description(desc="Doing final calculation                       ")

            if manualDebug:
                if manualDebug == query:
                    breakpoint()
                    pass
           
            finalFactorValue = finalFactor.get_value()
            
            totalAssignmentProbability = totalAssignmentProbability + (colourEvidenceProbability * finalFactorValue)
        
        totalAssignmentProbability = fixedDownwardZeroMultiplier * totalAssignmentProbability / downwardQueryConditionConstant
 
        conQueryEvaluation[query] = totalAssignmentProbability
        conQueryPreconditionConstants[query] = preconditionConstant

    #
    # Nonconjugate Queries:
    #

    ncQueryEvaluation = dict()

    if debug >= 1:
        print("Finished evaluation for conjugate queries. Starting evaluation of nonconjugate queries.")

    if progressBar:
        iterator = tqdm.tqdm(ncQueries.keys())
    else:
        iterator = ncQueries.keys()

    for query in iterator:
        if debug >= 2:
            print("Working on query " + query)
        if progressBar:
            iterator.set_description(desc="Working on query " + query)

        ncQueryEvaluation[query] = model.get_cpds(query).get_self_values_with_uncertain_evidence(**{parent:naiveProbabilities[parent] for parent in ncQueries[query]})[1][0]

    #
    # Full Queries
    #

    if debug >= 1:
        print("Finished nonconjugate queries. Merging together full queries.")

    if progressBar:
        iterator = tqdm.tqdm(fullQueries.keys())
    else:
        iterator = fullQueries.keys()

    fullQueryEvaluation = dict()
    for query in iterator:
        fullQueryEvaluation[query] = queryMergeHelper(conQueries, conQueryEvaluation, conQueryPreconditionConstants, query)

    if debug >= 1:
        print("All queries evaluated!")

    if save:
        if debug >= 1:
            print("Saving all query evaluations.")
        with open(modelFolder + modelName + "_modelevaluation" + modelExtension + "_completeConjugateQueries.pickle", "wb") as f:
            pickle.dump(conQueryEvaluation, f)
        with open(modelFolder + modelName + "_modelevaluation" + modelExtension + "_nonConjugateQueries.pickle", "wb") as f:
            pickle.dump(ncQueryEvaluation, f)
        with open(modelFolder + modelName + "_modelevaluation" + modelExtension + "_conjugateQueryPreconditionConstants.pickle", "wb") as f:
            pickle.dump(conQueryPreconditionConstants, f)
        with open(modelFolder + modelName + "_modelevaluation" + modelExtension + "_fullQueries.pickle", "wb") as f:
            pickle.dump(fullQueryEvaluation, f)

    return conQueryEvaluation, ncQueryEvaluation, conQueryPreconditionConstants, fullQueryEvaluation



    

   


