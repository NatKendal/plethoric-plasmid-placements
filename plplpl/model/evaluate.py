import collections
import itertools
import pickle

from pgmpy.factors import factor_product

from plplpl.NoisyOr import BinaryNoisyOrCPD
from plplpl.NoisyOr import NoisyOrBayesianNetwork

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
def assignmentHelper(model, node, endpoints, prefix, includeDeadEnds=False):
    if node in endpoints:
        return [[]]

    assignments = []
    children = [child for child in model.successors(node) if child[0] == prefix]
    subAssignments = []
    for child in children:
        subAssignments.append(assignmentHelper(model, child, endpoints, prefix))
    for subAssignment in itertools.product(*subAssignments):
        newAssignment = [node]
        for componentAssignment in subAssignment:
            for newNode in componentAssignment:
                newAssignment.append(newNode)
        assignments.append(newAssignment)
    if (len(assignments) > 0) or (includeDeadEnds):
        assignments.append([])
    return assignments

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
        iterator = tqdm.tqdm(conQueries.keys())
    else:
        iterator = conQueries.keys()

    for query in iterator:
        if debug >= 2:
            print(query + " - Initializing")
        if progressBar:
            iterator.set_description(desc=query + " - Initializing")
        
        # Find all the possible assignments to query and critical nodes.
        firstNode = sorted(conQueries[query]["query"] + conQueries[query]["critical"], key=lambda x: int(x.split("_")[1]))[0]
        QCzeroAssignments = assignmentHelper(model, firstNode, conQueries[query]["query"], "g")

        if progressBar:
            iterator.set_description(desc=query + " - Constructing simplified downward queries")
        # Build "connectedDown" query blocks.
        # We do some normalization now, so that it's used as a condition to check, instead of evidence to condition on.
        downwardQueryCPDs = list()
        downwardQueryEvidence = list()
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

        # For each assignment, build the query to evaluate.
        if progressBar:
            iterator.set_description(desc=query + " - Evaluating each assignment to gene nodes")
            nestedIterator = tqdm.tqdm([set(x) for x in QCzeroAssignments], leave=False)
        else:
            nestedIterator = [set(x) for x in QCzeroAssignments]

        totalAssignmentProbability = 0.0
        for assignment in nestedIterator:
            if progressBar:
                nestedIterator.set_description(desc="Calculating probability of colour evidence")
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
                    else:
                        colourAssignment[parent] = 1
                colourEvidenceProbability = colourEvidenceProbability * cpd.get_value(**colourAssignment)

            if progressBar:
                nestedIterator.set_description(desc="Calculating reduced cpds for critical region")
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
                    if (child[0] == "g") and (child not in assignment):
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
                            if int(grandChild.split("_")[1])-int(child.split("_")[1]) == model.constants["maturation_max"]:
                                maturationEndpoints.add(grandChild)
                        for successor in model.successors(node):
                            if successor[0] != "m":
                                continue
                            if int(successor.split("_")[1])-int(node.split("_")[1]) == model.constants["maturation_min"]:
                                maturationStartpoints.add(successor)

            if progressBar:
                nestedIterator.set_description(desc="Calculating each assignment to maturation nodes")

            # Get all fixed zero hidden nodes.
            fixedMaturationZeroes = set()
            queue = collections.deque(maturationStartpoints)
            while queue:
                node = queue.popleft()
                for parent in model.predecessors(node):
                    if (parent in conQueries[query]["hidden"]) and (parent not in fixedMaturationZeroes):
                        fixedMaturationZeroes.add(parent)
                        queue.append(parent)

            # Calculate possible assignments to the maturation nodes.
            """
            #DEBUGGING
            print(maturationStartpoints)
            print(maturationEndpoints)
            for x in maturationStartpoints:
                lower = 10000000
                upper = 0
                for y in assignmentHelper(model, x, maturationEndpoints, "m", includeDeadEnds=True):
                    for z in y:
                        if int(z.split("_")[1]) > upper:
                            upper = int(z.split("_")[1])
                        if int(z.split("_")[1]) < lower:
                            lower = int(z.split("_")[1])
                print(lower, upper)

            examine = [assignmentHelper(model, x, maturationEndpoints, "m", includeDeadEnds=True) for x in maturationStartpoints]
            breakpoint()
            """

            
            # Setup the factors that don't change between assignments.
            fixedToMarginalize = set()
            fixedFactors = list()
            for cpd in criticalNodeReducedCPDs:
                for parent in cpd.evidence:
                    toMarginalize.add(parent)
                factor = cpd.to_factor()
                factor.reduce([(cpd.variable, (0 if cpd.variable in assignment else 1))], inplace=True)
                fixedFactors.append(factor)
            for cpd in switchPointCPDs:
                for parent in cpd.evidence:
                    toMarginalize.add(parent)
                factor = cpd.to_factor()
                factor.reduce([(cpd.variable, 1)], inplace=True)
                fixedFactors.append(factor)

            maturationAssignments = list(itertools.product(*[assignmentHelper(model, x, maturationEndpoints, "m", includeDeadEnds=True) for x in maturationStartpoints]))

            if progressBar:
                nestedNestedIterator = tqdm.tqdm(maturationAssignments, leave=False)
            else:
                nestedNestedIterator = maturationAssignments

            overallProbability = 0.0
            nestedNestedCounter = 0
            for maturationAssignment in nestedNestedIterator:
                #if progressBar:
                #    nestedNestedIterator.set_description(desc="Calculating probability of assignment of maturation nodes")
                if progressBar:
                    nestedNestedCounter += 1
                    nestedNestedIterator.set_description(desc="(1) Working on " + str(nestedNestedCounter) + " of " + str(len(maturationAssignments)))
                # Build full list of hidden node zeroes for this maturation assignment.
                maturationZeroes = fixedMaturationZeroes.copy()
                for componentAssignment in maturationAssignment:
                    for node in componentAssignment:
                        maturationZeroes.add(node)
                maturationProbability = 1.0
                for maturationZero in maturationZeroes:
                    values = dict()
                    values[maturationZero] = 0
                    for parent in model.predecessors(maturationZero):
                        if parent[0] == "m":
                            values[parent] = 0
                        elif parent in assignment:
                            values[parent] = 0
                        else:
                            values[parent] = 1
                    maturationProbability = maturationProbability * model.get_cpds(maturationZero).get_value(**values)
                    # Check if this is right before a switchpoint.
                    # If so, calculate that too.
                    for child in model.successors(maturationZero):
                        if (child[0] == "m") and (child not in maturationZeroes):
                            values = dict()
                            values[child] = 1
                            for predecessor in model.predecessors(child):
                                if predecessor[0] == "m":
                                    values[predecessor] = 0
                                elif parent in assignment:
                                    values[predecessor] = 0
                                else:
                                    values[predecessor] = 1
                            maturationProbability = maturationProbability * model.get_cpds(child).get_value(**values)

                if progressBar:
                    nestedNestedIterator.set_description(desc="(2) Working on " + str(nestedNestedCounter) + " of " + str(len(maturationAssignments)))
                variableDownwardZeroProbability = 1.0
                for downwardZero in variableDownwardZero:
                    cpd = model.get_cpds(downwardZero)
                    for parent in model.predecessors(downwardZero):
                        if (parent in conQueries[query]["hidden"]) and (parent not in maturationZeroes):
                            variableDownwardZeroProbability = variableDownwardZeroProbability * (1.0 - cpd.evidence_noise[cpd.evidence.index(parent)])

                if progressBar:
                    nestedNestedIterator.set_description(desc="(3) Working on " + str(nestedNestedCounter) + " of " + str(len(maturationAssignments)))

                downwardQueryFactors = list()
                for i in range(len(downwardQueryCPDs)):
                    reducedBlockCPDs = list()
                    toReduce = list()
                    toMarginalize = list()
                    for cpd in downwardQueryCPDs[i]:
                        values = dict()
                        for j in range(len(cpd.evidence)):
                            if cpd.evidence[j] in conQueries[query]["hidden"]:
                                if cpd.evidence[j] in maturationZeroes:
                                    values[cpd.evidence[j]] = 0
                                else:
                                    values[cpd.evidence[j]] = 1
                        reducedBlockCPDs.append(cpd.reduce([(key, value) for key, value in values.items()], inplace=False))
                        if cpd.variable in downwardQueryEvidence[i]:
                            toReduce.append(cpd.variable)
                        else:
                            toMarginalize.append(cpd.variable)
                    factor = reducedBlockCPDs[0].to_factor()
                    for cpd in reducedBlockCPDs[1:]:
                        factor = factor * cpd.to_factor()
                    factor.reduce([(x, 1) for x in toReduce], inplace=True)
                    factor.marginalize(toMarginalize, inplace=True)
                    downwardQueryFactors.append(factor)

                if progressBar:
                    nestedNestedIterator.set_description(desc="(4) Working on " + str(nestedNestedCounter) + " of " + str(len(maturationAssignments)))

                # Compute the constant from all previously computed constants.
                finalProbability = fixedDownwardZeroMultiplier * colourEvidenceProbability * maturationProbability * variableDownwardZeroProbability

                allFactors = fixedFactors.copy()
                toMarginalize = fixedToMarginalize.copy()
               
                for factor in downwardQueryFactors:
                    for parent in factor.variables:
                        toMarginalize.add(parent)
                    allFactors.append(factor)

                # Eliminate all remaining variables.
                for variable in toMarginalize:
                    newFactors = list()
                    toMultiply = list()
                    for factor in allFactors:
                        if variable in factor.variables:
                            toMultiply.append(factor)
                        else:
                            newFactors.append(factor)
                    combinedFactor = factor_product(*toMultiply)
                    combinedFactor.marginalize([variable], inplace=True)
                    newFactors.append(combinedFactor)
                    allFactors = newFactors

                # All variables should be eliminated as this point.
                # All these factors should just be constants.
                for factor in allFactors:
                    finalProbability = finalProbability * factor.values

                overallProbability = overallProbability + finalProbability
            
            totalAssignmentProbability = totalAssignmentProbability + overallProbability
        
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



    

   


