import collections
import pickle

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
loadedEvidence: if we should use evidence already in memory instead of loading it.
loadedConjugateQueries: if we should use conjugate queries already in memory instead of loading them.
"""

def calculateNaiveProbabilities(modelFolder, dataFolder, modelName, modelExtension, save=True, debug=0, progressBar=False, loadedModel=None, loadedEvidence=None, loadedConjugateQueries=None):
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
        print("Loading evidence.")
    if loadedEvidence:
        evidence = loadedEvidence
    else:
        with open(modelFolder + modelName + "_modeldata" + modelExtension + "_evidence.pickle", "rb") as f:
            evidence = pickle.load(f)

    if debug >= 1:
        print("Loading Conjugate Queries.")
    if loadedConjugateQueries:
        conQueries = loadedConjugateQueries
    else:
        with open(modelFolder + modelName + "_modeldata" + modelExtension + "_completeConjugateQueries.pickle", "rb") as f:
            conQueries = pickle.load(f)

    # Final output dictionary.
    naiveProbabilities = dict()
    # Working dictionary for two stage approach.
    workingNaiveProbabilities = dict()
    # Stack of extra queries to compute in two stage approach.
    linkedQueries = set()
    # Queries examined so far in two stage approach.
    currentQueries = set()
    # Finished queries.
    finishedQueries = set()

    # Subprocess for recursion when calculating one node.
    # If force=True, then recompute even if the node is in naiveProbabilities or evidence.
    def doNode(node, destination, force=False):
        # If it's already done and we don't want to force, skip.
        if (node in naiveProbabilities) and (not force):
            return
        # If we have evidence and we don't want to force, skip.
        if (node in evidence) and (not force):
            naiveProbabilities[node] = evidence[node]
            return
        # If it's a maturation node downstream from a query and we don't want to force, find the recursive query and skip.
        if (node[0] == "m") and (not force):
            for parent in model.predecessors(node):
                if parent[0] == "g" and parent not in evidence:
                    for query in conQueries:
                        if parent in conQueries[query]["query"] or parent in conQueries[query]["critical"]:
                            linkedQueries.add(query)
                            return

        # Otherwise, compute it.
        arguments = dict()
        for parent in model.predecessors(node):
            if parent not in naiveProbabilities:
                doNode(parent, naiveProbabilities)
            if parent not in naiveProbabilities:
                arguments[parent] = 0
            else:
                arguments[parent] = naiveProbabilities[parent]
        destination[node] = model.get_cpds(node).get_self_values_with_uncertain_evidence(**arguments)[1][0]

    # Subprocess for handling a specific query.
    def doQuery(query, destination):
        currentQueries.add(query)
        if len(conQueries[query]["parentQueries"]) > 0:
            for parent in conQueries[query]["parentQueries"]:
                linkedQueries.add(query)
            return

        nodes = []
        networkedQueries = collections.deque()
        networkedQueries.append(query)
        allQueryNodes = []
        while networkedQueries:
            q = networkedQueries.popleft()
            currentQueries.add(q)
            for node in conQueries[q]["query"]:
                nodes.append(node)
            for node in conQueries[q]["critical"]:
                nodes.append(node)
            found = False
            for otherQuery in conQueries.keys():
                if q in conQueries[otherQuery]["parentQueries"]:
                    networkedQueries.append(otherQuery)
                    found = True
            if found == False:
                for node in conQueries[q]["query"]:
                    allQueryNodes.append(q)

        for node in sorted(nodes, key=lambda x: int(x.split("_")[1])):
            doNode(node, destination, force=True)

        # After we naively calculate, we normalize the query's naive probability to ensure that the end of the query is 1.
        finalProbability = min([destination[x] for x in allQueryNodes])
        # If the query is zero, then it can't have received the gene EXCEPT for from a linked query.
        # So we leave everything as zero for now and then recalculate later.
        if finalProbability == 0:
            return
        for node in nodes:
            destination[node] = min(destination[node]/finalProbability, 1)

    # First: Go through the Conjugate Queries, calculating them naively and all their dependences.
    if debug >= 1:
        print("Starting to naively calculate query nodes and dependences.")
    if progressBar:
        iterator = tqdm.tqdm(sorted(list(conQueries.keys()), key=lambda x: min([int(y.split("_")[1]) for y in conQueries[x]["query"]+conQueries[x]["critical"]])))
    else:
        iterator = sorted(list(conQueries.keys()), key=lambda x: min([int(y.split("_")[1]) for y in conQueries[x]["query"]+conQueries[x]["critical"]]))

    for query in iterator:
        if progressBar:
            iterator.set_description(desc="Working on query " + query)
        if debug >= 2:
            print("Working on query " + query)

        # Cleanup previously used dictionaries.
        linkedQueries = set()
        currentQueries = set()
        workingNaiveProbabilities = dict()

        # Calculate and normalize all linked queries independently.
        doQuery(query, naiveProbabilities)
        while len(linkedQueries) != 0:
            subQuery = linkedQueries.pop()
            if subQuery not in currentQueries:
                doQuery(subQuery, naiveProbabilities)

        # Calculate all the maturation node children from each linked query.
        for subQuery in currentQueries:
            for node in (conQueries[subQuery]["query"] + conQueries[subQuery]["critical"]):
                for child in model.successors(node):
                    if child[0] == "g":
                        continue
                    if child[0] == "c":
                        continue
                    doNode(child, naiveProbabilities, force=True)

        # Recalculate all linked queries using normalized results from independent case.
        for subQuery in currentQueries:
            # Recalculate critical region.
            doQuery(query, workingNaiveProbabilities)
        for subQuery in currentQueries:
            # Recalculate maturation node children.
            for node in (conQueries[query]["query"] + conQueries[query]["critical"]):
                for child in model.successors(node):
                    if child[0] == "g":
                        continue
                    if child[0] == "c":
                        continue
                    doNode(child, workingNaiveProbabilities, force=True)

        # Overwrite the independent naive query calculation with the recalculated one.
        for node in list(workingNaiveProbabilities.keys()):
            naiveProbabilities[node] = workingNaiveProbabilities.pop(node)

        # Note all the queries that have been handled.
        for subQuery in currentQueries:
            finishedQueries.add(subQuery)

    # Second: Calculate everything remaining.
    if debug >= 1:
        print("Finished naively calculating query nodes. Now calculating remaining nodes.")
    if progressBar:
        iterator = tqdm.tqdm(sorted(list(model.nodes), key=lambda x: int(x.split("_")[1])))
    else:
        iterator = sorted(model.nodes, key=lambda x: int(x.split("_")[1]))

    for node in iterator:
        if node not in naiveProbabilities:
            if progressBar:
                iterator.set_description(desc="Working on node " + node)
            if debug >= 2:
                print("Working on node " + node)

            doNode(node, naiveProbabilities)

    if save:
        if debug >= 1:
            print("Saving naive probabilities.")
            with open(modelFolder + modelName + "_modeldata" + modelExtension + "_naiveProbabilities.pickle", "wb") as f:
                pickle.dump(naiveProbabilities, f)

    return naiveProbabilities
