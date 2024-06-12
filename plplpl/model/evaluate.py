import collections
import itertools
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
loadedModel: if we should use a model already in memory instead of loading one.
loadedEvidence: 
loadedConjugateQueries:
loadedNonConjugateQueries:
loadedNaiveProbabilities:

This function calculates two things:
    For each identified non-conjugate query node gC_T:
        P(gC_T = 1 | naiveCalculation for all nodes of timestep < T), lower is better.

    For each identified conjugate query given as
        - "query": Query point gene nodes
        - "critical": Critical gene nodes
        - "unknown": Forced unknown nodes
        - "colourEvidence": Evidence colour nodes
        - "lineage": Fixed downstream consequence gene nodes
        - "hidden": Downstream maturation nodes to eliminate
        - "hardEvidence": Upstream nodes given as evidence that are assumed with probability 1
        - "incoming": Upstream maturation nodes given as naive probability that are assumed with probability 1
        - "connectedDown": Downstream consequence queries that are used to calculate likelihood of assignment
        - "downwardZero": Downstream nodes that evidence gives as zero, to contribute to likelihood of assignment to query
    We compute, given evidence E, for each valid assignment Q to query nodes and C to critical nodes:
        P("query" = Q, "critical" = C, "colourEvidence" = E, "lineage" = E, "downardZero" = E, collapsed("connectedDown") = E | "hardEvidence" = E), higher is better.
    "hidden" nodes are eliminated.
    "unknown" nodes are trivially eliminated (ignored).
    "incoming" nodes are eliminated.
    "downwardZero" evidence is preconditioned on the fact that it didn't get the plasmid from anywhere else.
    "connectedDown" evidence is rescaled so that if all incoming maturation nodes are 1, then it has maximum probability.
        - This is analogous to having this evidence as conditioning and then normalizing it out later, not quite the same, but much faster.
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
            return [block] + blocks, [None] + [(current if x == None else x) for x in parents]
        else:
            current = children[0]

def evaluateModel(modelFolder, dataFolder, modelName, modelExtension, save=True, debug=0, progressBar=False, loadedModel=None, loadedEvidence=None, loadedConjugateQueries=None, loadedNonConjugateQueries=None, loadedNaiveProbabilities=None):
    if progressBar:
        import tqdm

    if debug >= 1:
        print("Loading model.")
    if loadedModel:
        model = loadedModel
    else:
        with open(modelFolder + modelName + "_model" + modelExtension + "_contradictionsPruned_normalized.pickle", "rb") as f:
            model = pickle.load(f)

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

    if debug >= 1:
        print("Loading naive probabilities.")
    if loadedNaiveProbabilities:
        naiveProbabilities = loadedNaiveProbabilities
    else:
        with open(modelFolder + modelName + "_modeldata" + modelExtension + "_naiveProbabilities.pickle", "rb") as f:
            naiveProbabilities = pickle.load(f)

    #
    # Conjugate Queries:
    # 

    conQueryEvaluation = dict()

    if debug >= 1:
        print("Starting evaluation of conjugate queries.")

    if progressBar:
        iterator = tqdm.tqdm(sorted(conQueries.keys(), key=lambda x: int(x.split("_")[1]), reverse=True))
        #iterator = tqdm.tqdm(conQueries.keys())
    else:
        iterator = conQueries.keys()

    for query in iterator:
        #NoisyOrFactor._registry = dict()
        if debug >= 2:
            print(query + " - Initializing")
        if progressBar:
            iterator.set_description(desc=query + " - Initializing")

        # Should only come up in a child query of something split earlier.
        preconditionConstant = 1.0
        for node in conQueries[query]["zeroParent"]:
            parent = node
            while parent not in evidence:
                cpd = model.get_cpds(parent)
                for i in range(len(cpd.evidence)):
                    if cpd.evidence[i][0] == "g":
                        parent = cpd.evidence[i]
                    elif cpd.evidence[i][0] == "m":
                        preconditionConstant = preconditionConstant * (1 - (naiveProbabilities[cpd.evidence[i]] * cpd.evidence_noise[i]))
                    else:
                        print("Something weird happened. Shouldn't get here.")
        
        # Find all the possible assignments to query and critical nodes.
        firstNode = sorted(conQueries[query]["query"] + conQueries[query]["critical"], key=lambda x: int(x.split("_")[1]))[0]
        QCzeroAssignments = assignmentHelper(model, firstNode, conQueries[query]["query"], "g")

        if progressBar:
            iterator.set_description(desc=query + " - Constructing simplified downward queries")
        # Build "connectedDown" query blocks.
        downwardQueryCPDs = list()
        downwardQueryEvidence = list()
        downwardQueryConstant = 1.0
        for downwardQuery in conQueries[query]["connectedDown"].keys():
            fixed = True
            for hidden in conQueries[query]["connectedDown"][downwardQuery]:
                # Check if this query block only uses fixed maturation nodes.
                if fixed:
                    for parent in model.predecessors(hidden):
                        if parent[0] != "g":
                            continue
                        if parent not in conQueries[query]["lineage"]:
                            fixed = False
                            break
            # If so, it's not relevant to calculating the full query and we can skip. (It ends up trivial.)
            # This is maybe unjustified? Worth a discussion. Easy to calculate.
            if fixed:
                continue
            firstBlockNode = sorted(conQueries[downwardQuery]["query"] + conQueries[downwardQuery]["critical"], key=lambda x: int(x.split("_")[1]))[0]
            blocks, parents = blockBuilder(model, firstBlockNode, conQueries[downwardQuery]["query"])
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
                # Reduce to only query variables.
                values = dict()
                for variable in blockCPD.evidence:
                    if variable in conQueries[query]["hidden"]:
                        continue
                    if variable in conQueries[query]["incoming"]:
                        continue
                    if variable[0] == "g":
                        continue
                    values[variable] = naiveProbabilities[variable]
                blockCPD.reduce_with_uncertain_evidence([(key, value) for key, value in values.items()], inplace=True)
                blockCPDs.append(blockCPD)
            downwardQueryCPDs.append(blockCPDs)
            downwardQueryEvidence.append(conQueries[downwardQuery]["query"].copy())

        """
        downwardQueryFactors = []
        for i in range(len(downwardQueryCPDs)):
            toReduce = []
            toMarginalize = []
            for cpd in downwardQueryCPDs[i]:
                if cpd.variable in downwardQueryEvidence[i]:
                    toReduce.append((cpd.variable, 1))
                else:
                    toMarginalize.append(cpd.variable)
            factor = NoisyOrFactor("product", downwardQueryCPDs[i])
            downwardQueryFactors.append(NoisyOrFactor("marginalize", [NoisyOrFactor("reduce", [factor], toReduce)], argument=toMarginalize))
        """
        downwardQueryFactors = []
        downwardQueryToMarginalize = []
        for i in range(len(downwardQueryCPDs)):
            for cpd in downwardQueryCPDs[i]:
                if cpd.variable in downwardQueryEvidence[i]:
                    downwardQueryFactors.append(NoisyOrFactor("reduce", [cpd], argument=[(cpd.variable, 1)]))
                else:
                    downwardQueryToMarginalize.append(cpd.variable)
                    downwardQueryFactors.append(cpd)

        if progressBar:
            iterator.set_description(desc=query + " - Computing shared downward zeroes")

        # Calculate the downwardZero values that are shared for all assignments.
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
        if fixedDownwardZeroMultiplier == 0:
            breakpoint()

        # Precompute all variable downward zero values to be used as needed later.
        variableDownwardZeroPrecomputed = dict()
        for hidden in conQueries[query]["hidden"]:
            variableDownwardZeroPrecomputed[hidden] = 1.0
        for downwardZero in variableDownwardZero:
            cpd = model.get_cpds(downwardZero)
            for parent in model.predecessors(downwardZero):
                if parent in conQueries[query]["hidden"]:
                    variableDownwardZeroPrecomputed[hidden] = variableDownwardZeroPrecomputed[hidden] * (1.0 - cpd.evidence_noise[cpd.evidence.index(parent)])

        # Figure out which queries are incoming
        incomingQueries = dict()
        querylessIncoming = set()
        for incoming in conQueries[query]["incoming"]:
            for incomingQuery in conQueries.keys():
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
        #for incomingQuery in incomingQueries.keys():
        #    nodes = sorted(incomingQueries[incomingQuery], key=lambda x: int(x.split("_")[1]))
        # Instead, do it for everything, then use this to build "queries" for the queryless incoming
        for node in conQueries[query]["incoming"]:
            incomingChildren[node] = set()
            children = collections.deque(model.successors(node))
            while children:
                child = children.popleft()
                if child[0] != "m":
                    continue
                if child in evidence:
                    continue
                if child in conQueries[query]["incoming"]:
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
        for incoming in querylessIncoming:
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
            incomingQueries[incoming] = pseudoquery
        # All incoming nodes should be in incomingQueries now.

        # Build simplified CPDs for incoming hidden nodes.
        incomingCPDs = []
        for incoming in conQueries[query]["incoming"]:
            if incoming in incomingParents:
                alpha = naiveProbabilities[incomingParents[incoming]]
                beta = naiveProbabilities[incoming]
                #gamma = (beta - alpha)/(1.0 - alpha) #Floating point issues
                if (1.0 - alpha) <= 0.000000000001: # If 1-alpha is this close to zero, beta is almost 1 anyway.
                    gamma = beta
                else:
                    gamma = ((128 * beta) - (128 * alpha)) / (128 - (128 * alpha)) # Hopefully a scale up helps.
                incomingCPDs.append(BinaryNoisyOrCPD(incoming, gamma, evidence=[incomingParents[incoming]], evidence_noise=[1]))
            else:
                incomingCPDs.append(BinaryNoisyOrCPD(incoming, naiveProbabilities[incoming]))

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
                criticalNodeReducedCPDs.append(model.get_cpds(node).reduce([(key, value) for key, value in hardEvidence.items()], inplace=False))
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
                        switchPointCPDs.append(model.get_cpds(child).reduce([(key, value) for key, value in hardEvidence.items()], inplace=False))
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
                        switchPointCPDs.append(model.get_cpds(node).reduce([(key, value) for key, value in hardEvidence.items()], inplace=False))
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
            #fixedToMarginalize = set()
            """
            for cpd in criticalNodeReducedCPDs:
                #for parent in cpd.evidence:
                #    fixedToMarginalize.add(parent)
                #factor = cpd.to_factor()
                #factor.reduce([(cpd.variable, (0 if cpd.variable in assignment else 1))], inplace=True)
                factor = NoisyOrFactor("reduce", [cpd], argument=[(cpd.variable, (0 if cpd.variable in assignment else 1))])
                fixedFactors.append(factor)
            """

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
                workingFactors = []
                newFactors = []
                for factor in fixedFactors:
                    if len(set(incomingQueries[incomingQuery]) & set(factor.variables)) > 0:
                        workingFactors.append(factor)
                    else:
                        newFactors.append(factor)
                if len(workingFactors) == 0:
                    continue
                elif len(workingFactors) == 1:
                    productFactor = workingFactors[0]
                else:
                    productFactor = NoisyOrFactor("product", workingFactors)
                reducedFactors = []
                #print(list(itertools.product(*[assignmentHelperFromDictionary(incomingChildren, node) for node in incomingQueries[incomingQuery] if node not in incomingParents])))
                for m_assignment in itertools.product(*[assignmentHelperFromDictionary(incomingChildren, node) for node in incomingQueries[incomingQuery] if node not in incomingParents]):
                    zeroVariables = set()
                    for subAssignment in m_assignment:
                        for variable in subAssignment:
                            zeroVariables.add(variable)
                    constant = 1.0
                    # Add in constants from critical nodes set to zero, here.
                    for critical in assignment:
                        cpd = model.get_cpds(critical)
                        for i in range(len(cpd.evidence)):
                            if (cpd.evidence[i] in incomingQueries[incomingQuery]) and (cpd.evidence[i] not in zeroVariables):
                                constant = constant * (1 - cpd.evidence_noise[i])
                    if constant != 1:
                        reducedFactors.append(NoisyOrFactor("reduce", [NoisyOrFactor("product", [NoisyOrFactor("constant", [], argument=[constant]), productFactor])], argument=[(var, (0 if var in zeroVariables else 1)) for var in incomingQueries[incomingQuery]]))
                    else:
                        reducedFactors.append(NoisyOrFactor("reduce", [productFactor], argument=[(var, (0 if var in zeroVariables else 1)) for var in incomingQueries[incomingQuery]]))
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

            if progressBar:
                nestedIterator.set_description(desc="Doing final calculation                       ")
            
            finalFactorValue = finalFactor.get_value()
            """
            if finalFactorValue == 0:
                for factor in fixedFactors:
                    if len(factor.variables) == 0:
                         if factor.get_value() == 0:
                             break
                else:
                    breakpoint()
            """
            #print("Factor Value", finalFactorValue)
            #print("Colour Modifier", colourEvidenceProbability)
            
            totalAssignmentProbability = totalAssignmentProbability + (colourEvidenceProbability * finalFactorValue)
        #print("Total", totalAssignmentProbability)
        if len(conQueries[query]["parentQueries"]) == 0:
            if totalAssignmentProbability == 0:
                children = [x for x in conQueries.keys() if conQueries[x]["parentQueries"][0] == query]
                if (len(children) == 0) or (len([x for x in children if conQueryEvaluation[x] == 0]) > 0):
                    broken = False
                    for i in conQueries[query]["incoming"]:
                        if naiveProbabilities[i] > 0:
                            broken = True
                    for i in conQueries[query]["hardEvidence"]:
                        if i[0] == "m":
                            if evidence[i] == 1:
                                broken = True
                    if broken:
                        broken = False
                        for child in [x for x in children if conQueryEvaluation[x] == 0]:
                            for i in conQueries[query]["incoming"]:
                                if naiveProbabilities[i] > 0:
                                    broken = True
                            for i in conQueries[query]["hardEvidence"]:
                                if i[0] == "m":
                                    if evidence[i] == 1:
                                        broken = True
                        if broken:
                            print("Broken query!", query)
        totalAssignmentProbability = preconditionConstant * fixedDownwardZeroMultiplier * totalAssignmentProbability
        #print("Modified Total", totalAssignmentProbability)
        conQueryEvaluation[query] = totalAssignmentProbability

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

    if debug >= 1:
        print("All queries evaluated!")

    if save:
        if debug >= 1:
            print("Saving all query evaluations.")
        with open(modelFolder + modelName + "_modelevaluation" + modelExtension + "_completeConjugateQueries.pickle", "wb") as f:
            pickle.dump(conQueryEvaluation, f)
        with open(modelFolder + modelName + "_modelevaluation" + modelExtension + "_nonConjugateQueries.pickle", "wb") as f:
            pickle.dump(ncQueryEvaluation, f)

    return conQueryEvaluation, ncQueryEvaluation



    

   


