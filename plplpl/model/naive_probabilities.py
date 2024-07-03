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
loadedFullQueries: if we should use full queries already in memory instead of loading them.
"""

def calculateNaiveProbabilities(modelFolder, dataFolder, modelName, modelExtension, save=True, debug=0, progressBar=False, loadedModel=None, loadedEvidence=None, loadedConjugateQueries=None, loadedFullQueries=None):
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

    if debug >= 1:
        print("Loading Full Queries.")
    if loadedFullQueries:
        fullQueries = loadedFullQueries
    else:
        with open(modelFolder + modelName + "_modeldata" + modelExtension + "_fullQueries.pickle", "rb") as f:
            fullQueries = pickle.load(f)

    # Final output dictionary.
    naiveProbabilities = dict()
    # Stack of extra queries to compute in two stage approach.
    linkedQueries = set()
    # Queries examined so far in two stage approach.
    currentQueries = set()
    # Finished queries.
    finishedQueries = set()
    # Keeping track of maturation nodes.
    maturationNodes = dict()

    # Subprocess for recursion when calculating one node.
    # If force=True, then recompute even if the node is in naiveProbabilities or evidence.
    def doNode(node, destination, force=False, ignore={}):
        # If it's already done and we don't want to force, skip.
        if (node in naiveProbabilities) and (not force):
            return
        # If we have evidence and we don't want to force, just use the evidence value.
        if (node in evidence) and (not force):
            if node in conQueries: # Note that if it's in evidence because it's a query node, skip it.
                return
            naiveProbabilities[node] = evidence[node]
            return
        # If it's a maturation node downstream from a query and we don't want to force, find the recursive query and skip.
        if (node[0] == "m") and (not force):
            for parent in model.predecessors(node):
                if parent[0] == "g" and parent not in evidence:
                    for query in conQueries:
                        if parent in conQueries[query]["query"] or parent in conQueries[query]["critical"]:
                            parentQuery = query
                            while len(conQueries[parentQuery]["parentQueries"]) > 0:
                                parentQuery = conQueries[parentQuery]["parentQueries"][0]
                            linkedQueries.add(parentQuery)
                            return

        # Otherwise, compute it.
        arguments = dict()
        for parent in model.predecessors(node):
            if parent in ignore:
                arguments[parent] = 0
                continue

            if parent not in naiveProbabilities:
                doNode(parent, naiveProbabilities)
            if parent not in naiveProbabilities:
                arguments[parent] = 0
            else:
                arguments[parent] = naiveProbabilities[parent]
        destination[node] = model.get_cpds(node).get_self_values_with_uncertain_evidence(**arguments)[1][0]

    # Subprocess for handling a specific query.
    def doQuery(query, destination, ignore={}):
        currentQueries.add(query)

        nodes = set()
        allQueryNodes = set()
        for q in fullQueries[query]:
            for node in conQueries[q]["critical"]:
                nodes.add(node)
            for node in conQueries[q]["query"]:
                nodes.add(node)
                if len(conQueries[q]["childQueries"]) == 0:
                    allQueryNodes.add(node)

        # Include the gene nodes coming out of the query that disappear before evidence is given for them.
        for node in list(nodes):
            for child in model.successors(node):
                if (child[0] == "g") and (child not in evidence) and (child not in nodes):
                    queue = collections.deque()
                    queue.append(child)
                    while queue:
                        grandChild = queue.popleft()
                        nodes.add(grandChild)
                        for greatGrandChild in model.successors(grandChild):
                            if (greatGrandChild[0] == "g") and (greatGrandChild not in evidence) and (greatGrandChild not in nodes):
                                queue.append(greatGrandChild)


        for node in sorted(nodes, key=lambda x: int(x.split("_")[1])):
            doNode(node, destination, force=True, ignore=ignore)

        # After we naively calculate, we normalize the query's naive probability to ensure that the end of the query is 1.
        finalProbability = min([destination[x] for x in allQueryNodes])
        # If the query is zero, then it can't have received the gene EXCEPT for from a linked query.
        # So we leave everything as zero for now and then recalculate later.
        if finalProbability == 0:
            return
        for node in nodes:
            destination[node] = min(destination[node]/finalProbability, 1)

    # First: Go through the Conjugate Queries, calculating them naively and all their dependences.
    # Note: Going through the full queries and adding in all children to calculate at once.
    if debug >= 1:
        print("Starting to naively calculate query nodes and dependences.")
    if progressBar:
        iterator = tqdm.tqdm(sorted(list(fullQueries.keys()), key=lambda x: min([int(y.split("_")[1]) for y in conQueries[x]["query"]+conQueries[x]["critical"]])))
    else:
        iterator = sorted(list(fullQueries.keys()), key=lambda x: min([int(y.split("_")[1]) for y in conQueries[x]["query"]+conQueries[x]["critical"]]))

    for query in iterator:
        if progressBar:
            iterator.set_description(desc="Working on query " + query)
        if debug >= 2:
            print("Working on query " + query)

        if query in finishedQueries:
            continue

        # Cleanup previously used dictionaries.
        linkedQueries = set()
        currentQueries = set()

        # Calculate and normalize all linked queries independently.
        doQuery(query, naiveProbabilities)
        while len(linkedQueries) != 0:
            subQuery = linkedQueries.pop()
            if subQuery not in currentQueries:
                doQuery(subQuery, naiveProbabilities)

        # Calculate all the maturation node children from each linked query.
        for currentQuery in currentQueries:
            geneNodes = set()
            for subQuery in fullQueries[currentQuery]:
                for node in (conQueries[subQuery]["query"] + conQueries[subQuery]["critical"]):
                    geneNodes.add(node)

            # Include the gene nodes coming out of the query that disappear before evidence is given for them.
            for node in list(geneNodes):
                for child in model.successors(node):
                    if (child[0] == "g") and (child not in evidence) and (child not in geneNodes):
                        queue = collections.deque()
                        queue.append(child)
                        while queue:
                            grandChild = queue.popleft()
                            geneNodes.add(grandChild)
                            for greatGrandChild in model.successors(grandChild):
                                if (greatGrandChild[0] == "g") and (greatGrandChild not in evidence) and (greatGrandChild not in geneNodes):
                                    queue.append(greatGrandChild)

            nodes = set()
            for node in geneNodes:
                for child in model.successors(node):
                    if child[0] == "g":
                        continue
                    if child[0] == "c":
                        continue
                    nodes.add(child)

            maturationNodes[currentQuery] = nodes
        for currentQuery in currentQueries:
            for maturation in sorted(maturationNodes[currentQuery], key=lambda x: int(x.split("_")[1])):
                doNode(child, naiveProbabilities, force=True)

        # Recalculate all linked queries using normalized results from independent case.
        for subQuery in currentQueries:
            # Recalculate critical region.
            # We ignore the internal edges between two sibling nodes.
            doQuery(query, naiveProbabilities, ignore=maturationNodes[subQuery])

        for currentQuery in currentQueries:
            for maturation in sorted(maturationNodes[currentQuery], key=lambda x: int(x.split("_")[1])):
                doNode(child, naiveProbabilities, force=True)

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

            doNode(node, naiveProbabilities, force=True)


    # Third: Calculate incoming probabilities of queries.
    # i.e. probability given the normalized naive probability for every node and the precondition that your previous was 0.
    incomingProbabilities = dict()
    if debug >= 1:
        print("Finished naive probability calculation. Doing special incoming naive probability calculation for queries.")
    if progressBar:
        iterator = tqdm.tqdm(sorted(list(fullQueries.keys()), key=lambda x: min([int(y.split("_")[1]) for y in conQueries[x]["query"]+conQueries[x]["critical"]])))
    else:
        iterator = sorted(list(fullQueries.keys()), key=lambda x: min([int(y.split("_")[1]) for y in conQueries[x]["query"]+conQueries[x]["critical"]]))

    for query in iterator:
        if progressBar:
            iterator.set_description(desc="Working on query " + query)
        if debug >= 2:
            print("Working on query " + query)

        nodes = set()
        for subQuery in fullQueries[query]:
            for node in conQueries[subQuery]["critical"]:
                nodes.add(node)
            for node in conQueries[subQuery]["query"]:
                nodes.add(node)

        # Include the gene nodes coming out of the query that disappear before evidence is given for them.
        for node in list(nodes):
            for child in model.successors(node):
                if (child[0] == "g") and (child not in evidence) and (child not in nodes):
                    queue = collections.deque()
                    queue.append(child)
                    while queue:
                        grandChild = queue.popleft()
                        nodes.add(grandChild)
                        for greatGrandChild in model.successors(grandChild):
                            if (greatGrandChild[0] == "g") and (greatGrandChild not in evidence) and (greatGrandChild not in nodes):
                                queue.append(greatGrandChild)

        for node in nodes:
            arguments = dict()
            for parent in model.predecessors(node):
                if parent[0] == "g":
                    arguments[parent] = 0 # Assume previous node is zero.
                elif parent in maturationNodes[query]:
                    arguments[parent] = 0 # Discount internal edges of a query.
                else:
                    arguments[parent] = naiveProbabilities[parent]
            incomingProbabilities[node] = model.get_cpds(node).get_self_values_with_uncertain_evidence(**arguments)[1][0]

    # Fourth: Calculate incoming probabilities of everything else.
    if debug >= 1:
        print("Finished special incoming naive probability calculation for queries. Doing special incoming naive probability calculation for remaining nodes.")
    if progressBar:
        iterator = tqdm.tqdm(sorted(model.nodes, key=lambda x: int(x.split("_")[1])))
    else:
        iterator = sorted(model.nodes, key=lambda x: int(x.split("_")[1]))

    for node in iterator:
        if node in incomingProbabilities:
            continue

        if progressBar:
            iterator.set_description(desc="Working on node " + node)
        if debug >= 2:
            print("Working on node " + node)

        arguments = dict()
        for parent in model.predecessors(node):
            if parent[0] == node[0]:
                arguments[parent] = 0 # Assume previous node is zero.
            else:
                arguments[parent] = naiveProbabilities[parent]
        incomingProbabilities[node] = model.get_cpds(node).get_self_values_with_uncertain_evidence(**arguments)[1][0] 

    if save:
        if debug >= 1:
            print("Saving naive probabilities.")
        with open(modelFolder + modelName + "_modeldata" + modelExtension + "_naiveProbabilities.pickle", "wb") as f:
            pickle.dump(naiveProbabilities, f)
        with open(modelFolder + modelName + "_modeldata" + modelExtension + "_incomingProbabilities.pickle", "wb") as f:
            pickle.dump(incomingProbabilities, f)

    return naiveProbabilities, incomingProbabilities
