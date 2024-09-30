import collections
import pickle

from plplpl.NoisyOr import BayesianNetwork
from plplpl.NoisyOr import NoisyOrCPD

"""
Find all evidence that could be supplied to the model from the data.

modelFolder: path to directory
dataFolder: path to directory
modelName: unique name for this model/computation
modelExtension: model extension in the form `_[conjugation function]_[colour function]_[maturation function]` (don't include .pickle)
colour_min: minimum number of timesteps for colour to appear
colour_max: maximum number of timesteps for colour to appear
maturation_min: minimum number of timesteps for maturation (if None, reads from model)
maturation_max: maximum number of timesteps for maturation (if None, reads from model)
save: if we should save the model to a file (pickle)
debug: 0 = nothing, 1 = status, 2 = verbose
progressBar: if we should show a progress bar on long for loops
loadedModel: if we should use a model already in memory instead of loading one.
"""
def get_evidence(modelFolder, dataFolder, modelName, modelExtension, colour_min, colour_max, save=True, debug=0, progressBar=False, loadedModel=None):
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
        print("Getting constants.")

    assert colour_min >= 1
    assert colour_max >= colour_min
    model.constants["colour_min"] = colour_min
    model.constants["colour_max"] = colour_max
    maturation_min = model.constants["maturation_min"]
    maturation_max = model.constants["maturation_max"]

    # Update model extension.
    modelExtension = modelExtension.split("_")
    modelExtension[2] = str(colour_min) + "-" + str(colour_max)
    modelExtension = "_".join(modelExtension)

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

    if debug >= 1:
        print("Starting main evidence calculation.")

    directParent = model.getDirectParents()
    directChildren = model.getDirectChildren()

    # Store all the computed evidence.
    evidence = dict()

    # Start with the gene node evidence.
    if progressBar:
        iterator = tqdm.tqdm(model)
    else:
        iterator = model
    for node in iterator:
        if node[0] != "g":
            continue

        if debug >= 2:
            print("Working on " + node)
        if progressBar:
            iterator.set_description(desc="Working on " + node)

        colour = colours[uid[node[1:]]]

        if colour == 0: # Red
            evidence["g"+node[1:]] = 1 # Red cells have the gene.
            evidence["m"+node[1:]] = 1 # Red cells are mature.
        elif colour == 1: # Green
            if not directParent[node]:
                # If this is the first frame, double check that we don't have any close children that would imply we have the plasmid.
                children = [node]
                depth = 0
                green = True
                while depth <= colour_min:
                    newChildren = []
                    for child in children:
                        for grandChild in directChildren[child]:
                            if colours[uid[grandChild[1:]]] != 1:
                                green = False
                            else:
                                newChildren.append(grandChild)
                    children = newChildren
                    depth += 1
                if green:
                    # Force gene and maturation to be zero. 
                    evidence["g" + node[1:]] = 0
                    evidence["m" + node[1:]] = 0
                    continue
                else:
                    # Force gene and maturation to be one.
                    evidence["g" + node[1:]] = 1
                    evidence["m" + node[1:]] = 1
                    continue
            depth = 0
            parent = node
            while depth < colour_max:
                if directParent[parent]:
                    parent = directParent[parent]
                    depth += 1
                else:
                    break
            if depth == colour_max:
                # Green implies not having the gene at some point in the past.
                if parent in evidence and evidence[parent] == 1:
                    # print("CONTRADICTION IN EVIDENCE CALCULATION AT " + name[parent] + "!!!")
                    # Resolved by link snapping later. Interpreted as parent having gene but not giving it to child.
                    pass
                else:
                    evidence[parent] = 0     
        elif colour == 2: # Yellow
            if not directParent[node]:
                # If this is the first frame, also force gene and maturation to be one.
                evidence["g" + node[1:]] = 1
                evidence["m" + node[1:]] = 1
                continue
            depth = 0
            parent = node
            while depth < colour_min:
                if directParent[parent]:
                    parent = directParent[parent]
                    depth += 1
                else:
                    break
            if depth == colour_min:
                # Yellow implies having the gene at some point in the past.
                #if "g" + name[parent] in evidence and evidence["g" + name[parent]] == 0:
                    #print("CONTRADICTION IN EVIDENCE CALCULATION AT " + name[parent] + "!!!")
                    # Resolved by bias towards having the gene and not transferring it to child.
                evidence[parent] = 1

    if debug >= 1:
        print("Completed basic evidence. Now handling contradictions in model.")

    if progressBar:
        iterator = tqdm.tqdm(model)
    else:
        iterator = model

    count = 0
   
    for node in iterator:
        if node[0] != "g":
            continue

        if node not in evidence:
            continue
        
        if evidence[node] == 0:
            continue

        if debug >= 2:
            print("Working on " + node)
        if progressBar:
            iterator.set_description(desc="Working on " + node)

        # evidence[node] == 1
        for child in list(model.children(node)):
            if child[0] != "g":
                continue
            # If you have the gene, but you have a child that doesn't.
            if (child in evidence) and (evidence[child] == 0):
                # Note that we found a contradiction.
                count += 1
                # Begin snapping links.
                # First, treat this as new a starting point. (Maturation forced to zero.)
                evidence["m"+child[1:]] = 0

                # Determine the timepoint from which all previous links need to be deleted.
                disconnectTime = int(node.split("_")[1])

                # Update the gene node.
                model.remove_edge((node, child), fromCPD=True)
                directParent[child] = None

                # Update the maturation nodes.
                queue = collections.deque()
                queue.append("m" + child[1:])
                while queue:
                    currentNode = queue.popleft()
                    for parent in list(model.parents(currentNode)):
                        if int(parent.split("_")[1]) <= disconnectTime:
                            model.remove_edge((parent, currentNode), fromCPD=True)
                            if parent == directParent[currentNode]:
                                directParent[currentNode] = None
                    for grandChild in model.children(currentNode):
                        if grandChild[0] == currentNode[0]:
                            queue.append(grandChild)


    if debug >= 1:
        print("Finished handling contradictions in model. Found " + str(count) + ". Cleaning up edge cases.")

    if progressBar:
        iterator = tqdm.tqdm(model)
    else:
        iterator = model

    for node in iterator:
        if node[0] != "g":
            continue
        
        if debug >= 2:
            print("Working on " + node)
        if progressBar:
            iterator.set_description(desc="Working on " + node)

        if node in evidence:
            # Push zeroes upward. If you don't have it, you didn't have it before.
            if (evidence[node] == 0) and (directParent[node] not in evidence):
                parent = node
                while (directParent[parent]) and (directParent[parent] not in evidence):
                    evidence[directParent[parent]] = 0
                    parent = directParent[parent]
                if directParent[parent] and (directParent[parent] in evidence) and (directParent[evidence] == 1):
                    print("WARNING! FOUND CONTRADICTION AT " + str(parent))
            elif (evidence[node] == 1):
                # evidence[node] == 1
                # Push ones downward, if you have it, your child has it.
                # (Unless they definitely don't, those are the contradictions handled above.)
                queue = collections.deque()
                for child in model.children(node):
                    if (child[0] == node[0]) and (child not in evidence):
                        queue.append(child)
                while queue:
                    child = queue.popleft()
                    if child not in evidence:
                        evidence[child] = 1
                        for grandChild in model.children(child):
                            if (grandChild[0] == child[0]) and (grandChild not in evidence):
                                queue.append(grandChild)
                    elif evidence[child] == 0:
                        print("WARNING! FOUND CONTRADICTION AT " + str(child))

    if debug >= 1:
        print("Handled edge cases. Doing double check for contradictions.")

    if progressBar:
        iterator = tqdm.tqdm(model)
    else:
        iterator = model

    for node in iterator:
        if node[0] != "g":
            continue
        
        if debug >= 2:
            print("Working on " + node)
        if progressBar:
            iterator.set_description(desc="Working on " + node)

        if directParent[node] and (node in evidence) and (evidence[node] == 0) and (evidence[directParent[node]] == 1):
            print("FOUND CONTRADICTION AT " + str(directParent[node]) + " -> " + str(node))

    if debug >= 1:
        print("Finished double check. Starting to add maturation evidence.")
    if progressBar:
        iterator = tqdm.tqdm(model)
    else:
        iterator = model

    for node in iterator:
        if node[0] != "m":
            continue

        if node in evidence:
            continue

        if debug >= 2:
            print("Working on " + node)
        if progressBar:
            iterator.set_description(desc="Working on " + node)

        allOne = True
        allZero = True

        for parent in model.parents(node):
            if parent[0] == "m":
                continue
            if parent not in evidence:
                allOne = False
                allZero = False
                break
            if evidence[parent] == 1:
                allZero = False
            else:
                allOne = False

        # Special case where all theoretical parent nodes aren't present in our graph.
        # i.e. this cell recently appeared.
        if allZero and allOne:
            parent = directParent[node]
            while parent not in evidence:
                parent = directParent[parent]
            evidence[node] = evidence[parent]
        elif allZero:
            evidence[node] = 0
        elif allOne:
            evidence[node] = 1

    if debug >= 1:
        print("Finished evidence calculation.")
    if save:
        if debug >= 1:
            print("Saving model.")
        with open(modelFolder + modelName + "_model" + modelExtension + "_contradictionsPruned.pickle", "wb") as f:
            pickle.dump(model, f)
        with open(modelFolder + modelName + "_modeldata" + modelExtension + "_evidence.pickle", "wb") as f:
            pickle.dump(evidence, f) 

    return model, evidence




