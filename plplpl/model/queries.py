from collections import deque
import pickle

from plplpl.NoisyOr import NoisyOrBayesianNetwork

"""
modelFolder: path to directory
dataFolder: path to directory
modelName: unique name for this model/computation
modelExtension: model extension in the form `_[conjugation function]_[colour function]_[maturation function]` (don't include _contradictionsPruned.pickle)
save: if we should save the model to a file (pickle)
debug: 0 = nothing, 1 = status, 2 = verbose
progressBar: if we should show a progress bar on long for loops
loadedModel: if we should use a model already in memory instead of loading one.
loadedEvidence: if we should use an evidence already in memory instead of loading one.
loadedForwardLinks: if we should use a forwardLinks already in memory instead of loading one.
loadedBackwardLinks: if we should use a backwardLinks already in memory instead of loading one.

Returns dictionary of dictionaries
completeConjugateQueries[query] = {"query":queryVariables, "critical":criticalRegion, "evidence":queryEvidence, "unknown":queryForcedUnknown, "connectedDown":connectedDownwardQueries, "hidden":queryHidden,  "incoming":queryIncomingMature.difference(queryForcedUnknown).difference(queryHidden), "virtual":set(), "connectedUp":set(), "used":queryUsed, "required":set()}
    - query: all leaf unknowns to check probability of being 1.
    - critical: Core hidden variables to marginalize out.
    - evidence: Nodes relevant to the query that are given as evidence.
    - unknown: Nodes with evidence but the evidence is not given in order to make a counterfactual query.
    - connectedDown: Other queries which use a maturation node pointed to by a critical node in this query.
    - hidden: Non-Critical Hidden variables to marginalize out.
    - incoming: Incoming maturation nodes to critical variables to handle. Should be empty when dictionary is returned.
    - virtual: Incoming maturation nodes which can be given as virtual evidence once computed.
    - connectedUp: Other queries which have a maturation node pointed into a critical node of this query.
    - used: variables that are used in other queries in connectedDown
    - required: variables needed from queries in connectedUp


The goal of the conjugate queries is to compute the probability
P(Lineage receiving plasmid within critical range | not having plasmid before critical range AND observed status of rest of network) 
"""

def build_queries(modelFolder, dataFolder, modelName, modelExtension, save=True, debug=0, progressBar=False, loadedModel=None, loadedEvidence=None, loadedForwardLinks=None, loadedBackwardLinks=None):
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
        print("Loading forward links. (Pruned)")
    if loadedForwardLinks:
        forwardLinks = loadedForwardLinks
    else:
        with open(modelFolder + modelName + "_modeldata" + modelExtension + "_forwardLinksPostEvidence.pickle", "rb") as f:
            forwardLinks = pickle.load(f)

    if debug >= 1:
        print("Loading backward links. (Pruned)")
    if loadedBackwardLinks:
        backwardLinks = loadedBackwardLinks
    else:
        with open(modelFolder + modelName + "_modeldata" + modelExtension + "_backwardLinksPostEvidence.pickle", "rb") as f:
            backwardLinks = pickle.load(f)

    if debug >= 1:
        print("Loading evidence.")
    if loadedEvidence:
        evidence = loadedEvidence
    else:
        with open(modelFolder + modelName + "_modeldata" + modelExtension + "_evidence.pickle", "rb") as f:
            evidence = pickle.load(f)

    # Keep track of every edge in every query, for normalizing later.
    allEdges = set()

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
        if (node not in evidence) or (evidence[node] == 0) or (backwardLinks[uid[node[1:]]] == -1) or (("g" + name[backwardLinks[uid[node[1:]]]] in evidence) and (evidence["g" + name[backwardLinks[uid[node[1:]]]]] == 1)):
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
        iterator = tqdm.tqdm(list(conjugationQueries.keys()))
    else:
        iterator = list(conjugationQueries.keys())

    #print(list(conjugationQueries.keys()))

    absorbedQueries = set()
    completeConjugateQueries = dict()
    for query in iterator:
        if query in absorbedQueries:
            continue
        if debug >= 2:
            print("Working on building query " + query)
        if progressBar:
            iterator.set_description("Working on building query " + query)

        # The core query variables, the end of the critical range that evidence says must have the gene.
        queryVariables = set([query])
        # The remaining variables in the critical range.
        criticalRegion = set(conjugationQueries[query])
        # The colour nodes that constrain the possible assignments if they were included as evidence.
        queryForcedUnknown = set()
        # The colour nodes that are relevant evidence to the query.
        queryColourEvidence = set()
        # All gene nodes downwards from the query variables.
        queryLineage = set()
        # Hidden nodes in the query, maturation nodes that aren't fixed by evidence.
        queryHidden = set()
        # Nodes available in evidence for which the evidence is taken as fact. CPDs are implicitly changed.
        queryHardEvidence = set()
        # Maturation nodes not in evidence that point into gene nodes of the critical region/query variables.
        queryIncomingMature = set()
        # Other queries that this query's hidden nodes point into.
        # These will be collapsed later using naive priors.
        connectedDownwardQueries = set()
        # Gene nodes that have evidence = 0 but have a (potentially) active maturation node pointing into them.
        # These have CPD implicitly reformulated to be exclusively in terms of relevant nodes.
        queryDownwardGeneZero = set()
 
        # First: handle the main query node.
        for child in model.successors(query):
            if child[0] == "c":
                queryForcedUnknown.add(child)
            elif child[0] == "g":
                children = deque()
                children.append(uid[child[1:]])
                while children:
                    grandChild = children.popleft()
                    queryLineage.add("g" + name[grandChild])
                    for greatGrandChild in forwardLinks[grandChild]:
                        children.append(greatGrandChild)
            elif child[0] == "m":
                queryHidden.add(child)
            else:
                print("Something went wrong. Child wasn't a 'c', 'g', or 'm' node.")
                print(child, query)

        for parent in model.predecessors(query):
            if parent[0] == "g":
                continue
            elif parent in evidence:
                allEdges.add((parent, query))
                queryHardEvidence.add(parent)
            else:
                allEdges.add((parent, query))
                queryIncomingMature.add(parent)

        # Second: get all queries we will absorb.
        queue = deque(conjugationQueries[query])
        while queue:
            critical = queue.popleft()
            for child in model.successors(critical):
                if child in queryForcedUnknown:
                    continue
                elif child[0] != "g":
                    continue
                elif child in criticalRegion:
                    continue
                elif child in queryVariables:
                    continue
                elif child not in evidence:
                    criticalRegion.add(child)
                    queue.append(child)
                else:
                    # Found a query to absorb
                    if evidence[child] == 0:
                        print("Model leads to a contradictory query around " + query + " " + critical + " " + child)
                        raise AssertionError("Model leads to a contradictory query around " + query + " " + critical + " " + child)
                    queryVariables.add(child)

                    # Fix other queries that point to absorbed query.
                    absorbedQueries.add(child)
                    for otherQuery in completeConjugateQueries.keys():
                        if otherQuery == query:
                            pass
                        if child in completeConjugateQueries[otherQuery]["connectedDown"]:
                            completeConjugateQueries[otherQuery]["connectedDown"].remove(child)
                            completeConjugateQueries[otherQuery]["connectedDown"].add(query)

                    for grandChild in model.successors(child):
                        if grandChild[0] == "c":
                            queryForcedUnknown.add(grandChild)
                        elif grandChild[0] == "g":
                            children = deque()
                            children.append(uid[child[1:]])
                            while children:
                                successor = children.popleft()
                                queryLineage.add("g" + name[successor])
                                for grandSuccessor in forwardLinks[successor]:
                                    children.append(grandSuccessor)
                        elif grandChild[0] == "m":
                            queryHidden.add(child)
                        else:
                            print("Something went wrong. Child wasn't a 'c', 'g', or 'm' node.")
                    for predecessor in model.predecessors(child):
                        if predecessor[0] == "g":
                            continue
                        elif predecessor in evidence:
                            allEdges.add((predecessor, child))
                            queryHardEvidence.add(predecessor)
                        else:
                            allEdges.add((predecessor, child))
                            queryIncomingMature.add(predecessor)

        # Third: we handle the complete critical region.
        queue = deque(criticalRegion)
        while queue:
            critical = queue.popleft()
            for child in model.successors(critical):
                if child in queryForcedUnknown:
                    continue
                if child[0] == "g":
                    continue
                elif child[0] == "c":
                    queryColourEvidence.add(child)
                elif child[0] == "m":
                    if child in evidence:
                        print("Model expected no evidence but found it at query " + query + " " + critical + " " + child)
                        raise ValueError("Model expected no evidence but found it at query " + query + " " + critical + " " + child)
                    queryHidden.add(child)
            for parent in model.predecessors(critical):
                if parent in criticalRegion:
                    continue
                elif parent in evidence:
                    if parent[0] == "m":
                        allEdges.add((parent, critical))
                    queryHardEvidence.add(parent)
                elif parent[0] == "m":
                    allEdges.add((parent, critical))
                    queryIncomingMature.add(parent)
                else:
                    print("Shouldn't get here. Something went wrong when looking at " + query + " " + critical + " " + parent)
                    raise AssertionError("Shouldn't get here. Something went wrong when looking at " + query + " " + critical + " " + parent)

        # Fourth: get all relevant downstream consequences.
        for node in queryLineage:
            for child in model.successors(node):
                if child[0] == "g":
                    continue
                elif child[0] == "c":
                    queryForcedUnknown.add(child)
                elif child[0] == "m":
                    queryHidden.add(child)
                else:
                    print("Something went wrong, child " + child + " wasn't a 'c', 'g', or 'm' node.")

        # Fifth: We expand all the hidden nodes.
        for hidden in queryHidden:
            for grandChild in model.successors(hidden):
                if grandChild[0] != "g":
                    continue
                if grandChild in evidence:
                    if evidence[grandChild] == 0:
                        # If it's a zero, then other incoming events aren't correlated.
                        allEdges.add((hidden, grandChild))
                        queryDownwardGeneZero.add(grandChild)
                        continue
                    if backwardLinks[uid[grandChild[1:]]] == -1:
                        # This shouldn't happen, it occurred due to a bug that was thought to be resolved.
                        print("The problem, where we add a connection where it shouldn't, came up again.")
                    if ("g" + name[backwardLinks[uid[grandChild[1:]]]] in evidence) and (evidence["g" + name[backwardLinks[uid[grandChild[1:]]]]] == 1):
                        # If its parent had the gene, then incoming events are of no consequence.
                        continue
                    # evidence[grandChild] = 1 and parent doesn't. It's a query node exactly.
                    allEdges.add((hidden, grandChild))
                    for otherQuery in conjugationQueries.keys():
                        if grandChild in conjugationQueries[otherQuery] or grandChild == otherQuery:
                            target = otherQuery
                            if otherQuery in queryVariables:
                                target = query
                            elif otherQuery in absorbedQueries:
                                for key in completeConjugateQueries.keys():
                                    if otherQuery in completeConjugateQueries[key]["query"]:
                                        target = key
                                        break
                            if (target != query):
                                connectedDownwardQueries.add(target)
                            else:
                                # This is a self loop, but that's probably fine.
                                # Previously threw an error, but this case *does* occur.
                                # Need to ensure that such cases are handled properly later.
                                pass
                            break
                    else:
                        print("Node " + grandChild + " looks like a query node, but isn't.")
                        raise AssertionError("Node " + grandChild + " looks like a query node, but isn't.")
                else:
                    # This node is in the critical region of another query.
                    allEdges.add((hidden, grandChild))
                    for otherQuery in conjugationQueries.keys():
                        if grandChild in conjugationQueries[otherQuery]:
                            target = otherQuery
                            if otherQuery in queryVariables:
                                target = query
                            elif otherQuery in absorbedQueries:
                                for key in completeConjugateQueries.keys():
                                    if otherQuery in completeConjugateQueries[key]["query"]:
                                        target = key
                                        break
                            if (target != query):
                                connectedDownwardQueries.add(target)
                            else:
                                # This is a self loop, but that's probably fine.
                                # Previously threw an error but this case *does* occur.
                                # Need to ensure that such cases are handled properly later.
                                pass
                            break
                    else:
                        print("At " + hidden + " we have that " + grandChild + " should be in a critical region, but isn't.")
                        raise AssertionError("At " + hidden + " we have that " + grandChild + " should be in a critical region, but isn't.")
    
        completeConjugateQueries[query] = {"query":queryVariables, "critical":criticalRegion, "unknown":queryForcedUnknown, "colourEvidence":queryColourEvidence, "lineage":queryLineage, "hidden":queryHidden, "hardEvidence":queryHardEvidence, "incoming":queryIncomingMature, "connectedDown":connectedDownwardQueries, "downwardZero":queryDownwardGeneZero}

    for query in list(completeConjugateQueries.keys()):
        for key in completeConjugateQueries[query].keys():
            completeConjugateQueries[query][key] = list(completeConjugateQueries[query][key])

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
        nonConjugateQueries[node] = []
        for parent in model.predecessors(node):
            if parent[0] == "m":
                allEdges.add((parent, node))
            nonConjugateQueries[node].append(parent)

    if save:
        if debug >= 1:
            print("Saving conjugate queries, nonconjugate queries, and all relevant edges.")
        with open(modelFolder + modelName + "_modeldata" + modelExtension + "_completeConjugateQueries.pickle", "wb") as f:
            pickle.dump(completeConjugateQueries, f)
        with open(modelFolder + modelName + "_modeldata" + modelExtension + "_nonConjugateQueries.pickle", "wb") as f:
            pickle.dump(nonConjugateQueries, f)
        with open(modelFolder + modelName + "_modeldata" + modelExtension + "_queryEdgeList.pickle", "wb") as f:
            pickle.dump(sorted(list(allEdges)), f)

    return completeConjugateQueries, nonConjugateQueries, sorted(list(allEdges))
