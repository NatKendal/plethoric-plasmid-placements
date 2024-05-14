from collections import deque
import pickle

from plplpl.NoisyOr import NoisyOrBayesianNetwork

"""
modelFolder: path to directory
dataFolder: path to directory
modelName: unique name for this model/computation
modelExtension: model extension in the form `_[conjugation function]_[colour function]_[maturation function]_contradictionsPruned` (don't include .pickle)
save: if we should save the model to a file (pickle)
debug: 0 = nothing, 1 = status, 2 = verbose
progressBar: if we should show a progress bar on long for loops
loadedModel: if we should use a model already in memory instead of loading one.
loadedEvidence: if we should use an evidence already in memory instead of loading one.
loadedForwardLinks: if we should use a forwardLinks already in memory instead of loading one.
loadedBackwardLinks: if we should use a backwardLinks already in memory instead of loading one.

Returns dictionary of queries
completeConjugationQueries[queryCoreVariable] = [list(queryVariables), list(criticalRegion), list(queryEvidence), list(queryHidden), list(queryForcedUnknown), list(connectedQueries), list(queryRequired)]
    - Query variables: all leaf unknowns to check probability of 1.
    - Critical Region: Core hidden variables to marginalize out.
    - Query Evidence: Nodes relevant to the query that are given as evidence.
    - Query Hidden: Non-Critical Hidden variables to marginalize out.
    - Query Forced Unknown: Nodes with evidence but the evidence is not given in order to make a counterfactual query.
    - Connected Queries: Query core variables to other queries not independent of this one.
    - Query Required: Variables that are not in a query but are unknown. Typically weird stuff near start or end.
"""

# NEXT TODO:
# Query Required shouldn't exist. Sit down and do the recursion until you find which query things are related to, or show that they aren't related.
def build_queries(modelFolder, dataFolder, modelName, modelExtension, save=True, debug=0, progressBar=False, loadedModel=None, loadedEvidence=None, loadedForwardLinks=None, loadedBackwardLinks=None):
    if progressBar:
        import tqdm

    if debug >= 1:
        print("Loading model.")
    if loadedModel:
        model = loadedModel
    else:
        with open(modelFolder + modelName + "_model" + modelExtension + ".pickle", "rb") as f:
            model = pickle.load(f)

    if debug >= 1:
        print("Loading cell names and unique ids.")
    with open(dataFolder + modelName + "_humanFriendlyNameLookup.pickle", "rb") as f:
        uid = pickle.load(f)
    with open(dataFolder + modelName + "_humanFriendlyName.pickle", "rb") as f:
        name = pickle.load(f)

    if debug >= 1:
        print("Loading backward links.")
    with open(dataFolder + modelName + "_backwardLinks.pickle", "rb") as f:
        backwardLinks = pickle.load(f)
    
    if debug >= 1:
        print("Loading forward links. (Pruned)")
    if loadedForwardLinks:
        forwardLinks = loadedForwardLinks
    else:
        with open(dataFolder + modelName + "_forwardLinksPostEvidence.pickle", "rb") as f:
            forwardLinks = pickle.load(f)

    if debug >= 1:
        print("Loading backward links. (Pruned)")
    if loadedBackwardLinks:
        backwardLinks = loadedBackwardLinks
    else:
        with open(dataFolder + modelName + "_backwardLinksPostEvidence.pickle", "rb") as f:
            backwardLinks = pickle.load(f)

    if debug >= 1:
        print("Loading first events.")
    with open(dataFolder + modelName + "_firsts.pickle", "rb") as f:
        firsts = pickle.load(f)

    if debug >= 1:
        print("Loading evidence.")
    if loadedEvidence:
        evidence = loadedEvidence
    else:
        with open(modelFolder + modelName + "_model" + modelExtension + "_evidence.pickle", "rb") as f:
            evidence = pickle.load(f)

    #
    # False Negative Queries
    #

    conjugationQueries = dict()
    if debug >= 1:
        print("Finding conjugation queries.")
    if progressBar:
        iterator = tqdm.tqdm(model)
    else:
        iterator = model
    for node in iterator:
        if node[0] != "g":
            continue
        if (node not in evidence) or (evidence[node] == 0) or (backwardLinks[uid[node[1:]]] == -1) or ("g" + name[backwardLinks[uid[node[1:]]]] in evidence):
            continue
        # If it passed the previous checks, then it is gene node that is 1 and its parent is unknown.
        # This is a query to examine.
        if debug >= 2:
            print("Building a query around " + node)
        if progressBar:
            iterator.set_description(desc="Building a query around " + node)

        hiddenNodes = []
        parent = uid[node[1:]]
        while (backwardLinks[parent] != -1) and ("g" + name[backwardLinks[parent]] not in evidence):
            parent = backwardLinks[parent]
            hiddenNodes.append("g" + name[parent])
        if (backwardLinks[parent] == -1):
            if debug >= 2:
                print("Skipping query since it leads to a dead end.")
            continue
        if evidence["g" + name[backwardLinks[parent]]] == 1:
            print("Something went wrong when calculating query on " + node)
            raise ValueError("Model leads to a malformed query around " + node)
        # node and hiddenNodes define a critical region for a query.
        # conjugationQueries[queryNode] = [[critical region], [evidence nodes], [interfering queries]]
        conjugationQueries[node] = hiddenNodes

    if debug >= 1:
        print("Found " + str(len(conjugationQueries)) + " conjugation queries. Calculating full conjugation queries now.")
    if progressBar:
        iterator = tqdm.tqdm(conjugationQueries.keys())
    else:
        iterator = conjugationQueries.keys()

    absorbedQueries = set()
    completeConjugateQueries = dict()
    for query in iterator:
        if query in absorbedQueries:
            continue
        if debug >= 2:
            print("Working on building query " + query)
        if progressBar:
            iterator.set_description("Working on building query " + query)

        queryVariables = set([query])
        criticalRegion = set(conjugationQueries[query])
        queryEvidence = set()
        queryForcedUnknown = set()
        queryHidden = set()
        connectedQueries = set()
        queryRequired = set()

        for child in model.successors(query):
            queryForcedUnknown.add(child)
        for parent in model.predecessors(query):
            if parent[0] == "g":
                continue
            elif parent in evidence:
                queryEvidence.add(parent)
            else:
                for grandParent in model.predecessors(parent):
                    if grandParent in evidence:
                        queryEvidence.add(grandParent)
                    else:
                        for otherQuery in conjugationQueries.keys():
                            if grandParent in conjugationQueries[otherQuery]:
                                if otherQuery != query:
                                    connectedQueries.add(otherQuery)
                                break
                        else:
                            queryRequired.add(grandParent)
                queryHidden.add(parent)

        queue = deque(conjugationQueries[query])
        while queue:
            critical = queue.popleft()
            for child in model.successors(critical):
                if child in queryForcedUnknown:
                    continue
                if child[0] == "g":
                    if (child not in criticalRegion) and (child not in queryVariables):
                        if child in evidence:
                            if evidence[child] == 0:
                                print("Model leads to a contradictory query around " + query + " " + critical + " " + child)
                                raise ValueError("Model leads to a contradictory query around " + query + " " + critical + " " + child)
                            # evidence[child] == 1
                            queryVariables.add(child)
                            absorbedQueries.add(child)
                            for grandChild in model.successors(child):
                                queryForcedUnknown.add(grandChild)
                            criticalRegion.difference(queryForcedUnknown)
                            queryEvidence.difference(queryForcedUnknown)
                            queryHidden.difference(queryForcedUnknown)
                            queryRequired.difference(queryForcedUnknown)
                            for predecessor in model.predecessors(child):
                                if predecessor[0] == "g":
                                    continue
                                elif predecessor in evidence:
                                    queryEvidence.add(predecessor)
                                else:
                                    for grandPredecessor in model.predecessors(predecessor):
                                        if grandPredecessor in evidence:
                                            queryEvidence.add(grandPredecessor)
                                        else:
                                            for otherQuery in conjugationQueries.keys():
                                                if grandPredecessor in conjugationQueries[otherQuery]:
                                                    if otherQuery != query:
                                                        connectedQueries.add(otherQuery)
                                                    break
                                            else:
                                                queryRequired.add(grandPredecessor)
                                    queryHidden.add(predecessor)
                        else:
                            criticalRegion.add(child)
                            queue.append(child)
                elif child[0] == "m":
                    if child in evidence:
                        print("Model expected no evidence but found it at query " + query + " " + critical + " " + child)
                        raise ValueError("Model expected no evidence but found it at query " + query + " " + critical + " " + child)
                    for grandChild in model.successors(child):
                        if grandChild[0] != "g":
                            continue
                        if grandChild in evidence:
                            queryEvidence.add(grandChild)
                            if evidence[grandChild] == 0:
                                # If it's a zero, then other incoming events aren't correlated.
                                continue
                            if backwardLinks[uid[grandChild[1:]]] == -1:
                                # This shouldn't happen, but in case it does, we can just ignore it.
                                continue
                            if "g" + name[backwardLinks[uid[grandChild[1:]]]] in evidence and evidence["g" + name[backwardLinks[uid[grandChild[1:]]]]]:
                                # If its parent had the gene, then incoming events are of no consequence.
                                continue
                            for predecessor in model.predecessors(grandChild):
                                if predecessor == child:
                                    continue
                                if predecessor in evidence:
                                    queryEvidence.add(predecessor)
                                else:
                                    for grandPredecessor in model.predecessors(predecessor):
                                        if grandPredecessor in evidence:
                                            queryEvidence.add(grandPredecessor)
                                        else:
                                            for otherQuery in conjugationQueries.keys():
                                                if grandPredecessor in conjugationQueries[otherQuery]:
                                                    if otherQuery != query:
                                                        connectedQueries.add(otherQuery)
                                                    break
                                            else:
                                                queryRequired.add(grandPredecessor)
                                    queryHidden.add(predecessor)
                        else:
                            for otherQuery in conjugationQueries.keys():
                                if grandChild in conjugationQueries[otherQuery]:
                                    if otherQuery != query:
                                        connectedQueries.add(otherQuery)
                                    break
                            else:
                                # If this path leads to the unknown, no dependence flows through here.
                                pass
                    queryHidden.add(child)
                else:
                    # child[0] == "c"
                    queryEvidence.add(child)
            for parent in model.predecessors(critical):
                if parent in criticalRegion:
                    continue
                elif parent in evidence:
                    queryEvidence.add(parent)
                elif parent[0] == "m":
                    for grandParent in model.predecessors(parent):
                        if grandParent in evidence:
                            queryEvidence.add(grandParent)
                        else:
                            for otherQuery in conjugationQueries.keys():
                                if grandParent in conjugationQueries[otherQuery]:
                                    if otherQuery != query:
                                        connectedQueries.add(otherQuery)
                                    break
                            else:
                                queryRequired.add(grandParent)
                else:
                    print("Shouldn't get here. Something went wrong when looking at " + query + " " + critical + " " + parent)
                    raise ValueError("Shouldn't get here. Something went wrong when looking at " + query + " " + critical + " " + parent)

        completeConjugateQueries[query] = [list(queryVariables), list(criticalRegion), list(queryEvidence), list(queryHidden), list(queryForcedUnknown), list(connectedQueries), list(queryRequired)]

    #
    # False Positive Queries:
    #

    if debug >= 1:
        print("Finished making complete conjugation queries. Making non-conjugation queries.")

    nonConjugateQueries = dict()
    if progressBar:
        iterator = tqdm.tqdm(model)
    else:
        iterator = model

    for node in iterator:
        if node[0] != "g":
            continue
        if node not in evidence:
            continue
        if evidence[node] == 1:
            continue
        # evidence[node] == 0
        if backwardLinks[uid[node[1:]]] == -1:
            continue
        
        if debug >= 2:
            print("Working on " + node)
        if progressBar:
            iterator.set_description(desc="Working on " + node)
        # We should check false positive chance of node.
        nonConjugateQueries[node] = [list(model.predecessors(node))]

    if save:
        with open(modelFolder + modelName + "_model" + modelExtension + "_completeConjugateQueries.pickle", "wb") as f:
            pickle.dump(completeConjugateQueries, f)
        with open(modelFolder + modelName + "_model" + modelExtension + "_nonConjugateQueries.pickle", "wb") as f:
            pickle.dump(nonConjugateQueries, f)

    return completeConjugateQueries, nonConjugateQueries
