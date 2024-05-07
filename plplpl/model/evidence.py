import pickle

from plplpl.NoisyOr import NoisyOrBayesianNetwork

"""
Find all evidence that could be supplied to the model from the data.

modelFolder: path to directory
dataFolder: path to directory
modelName: unique name for this model/computation
modelExtension: model extension in the form `_[conjugation function]_[colour function]_[maturation function]` (don't include .pickle)
colour_min: minimum number of timesteps for colour to appear (if None, reads from model)
colour_max: maximum number of timesteps for colour to appear (if None, reads from model)
maturation_min: minimum number of timesteps for maturation (if None, reads from model)
maturation_max: maximum number of timesteps for maturation (if None, reads from model)
save: if we should save the model to a file (pickle)
debug: 0 = nothing, 1 = status, 2 = verbose
progressBar: if we should show a progress bar on long for loops
loadedModel: if we should use a model already in memory instead of loading one.
"""
#TESTING
# e = get_evidence("data/trap6test/", "data/trap6test/trap6test_data/", "trap6test", "_contactWeightsFixedNaive_colourDelayFunctionUniform_maturationDelayFunctionNormal", save=True, debug=1, progressBar=True)
# e = get_evidence("data/trap6test/", "data/trap6test/trap6test_data/", "trap6test", "_contactWeightsFixedNaive_colourDelayFunctionUniform_maturationDelayFunctionNormal", save=True, debug=1, progressBar=True, loadedModel=model)
def get_evidence(modelFolder, dataFolder, modelName, modelExtension, colour_min=None, colour_max=None, maturation_min=None, maturation_max=None, save=True, debug=0, progressBar=False, loadedModel=None):
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
    if not colour_min:
        colour_min = model.constants["colour_min"]
    if not colour_max:
        colour_max = model.constants["colour_max"]
    if not maturation_min:
        maturation_min = model.constants["maturation_min"]
    if not maturation_max:
        maturation_max = model.constants["maturation_max"]

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
        print("Loading backwards links.")
    with open(dataFolder + modelName + "_backwardLinks.pickle", "rb") as f:
        backwardLinks = pickle.load(f)

    if debug >= 1:
        print("Loading forward links.")
    with open(dataFolder + modelName + "_forwardLinks.pickle", "rb") as f:
        forwardLinks = pickle.load(f)

    if debug >= 1:
        print("Starting main evidence calculation.")

    evidence = dict()
    if progressBar:
        iterator = tqdm.tqdm(model)
    else:
        iterator = model
    for node in iterator:
        if node[0] != "c":
            continue

        if debug >= 2:
            print("Working on " + node)
        if progressBar:
            iterator.set_description(desc="Working on " + node)

        colour = colours[uid[node[1:]]]
        if colour == 0: # Red
            evidence["g"+node[1:]] = 1 # Red cells have the gene.
            evidence[node] = 1 # Red cells are coloured.
            evidence["m"+node[1:]] = 1 # Red cells are mature.
        if colour == 1: # Green
            evidence[node] = 0 # Green cells are not coloured.
            depth = 0
            parent = uid[node[1:]]
            seen = set()
            seen.add(parent)
            while depth < colour_max:
                if backwardLinks[parent] != -1:
                    parent = backwardLinks[parent]
                    seen.add(parent)
                    depth += 1
                else:
                    break
            if depth == colour_max:
                # Green implies not having the gene at some point in the past.
                if "g" + name[parent] in evidence and evidence["g" + name[parent]] == 1:
                    # print("CONTRADICTION IN EVIDENCE CALCULATION AT " + name[parent] + "!!!")
                    # Resolved by bias towards having the gene and not transferring it to child.
                    pass
                else:
                    evidence["g" + name[parent]] = 0
            start = colour_max - depth
            if start <= maturation_min - 1:
                children = [parent]
                for _ in range(start, maturation_min-1):
                    newChildren = []
                    skipCheck = (len(seen) == 0)
                    for child in children:
                        for grandChild in forwardLinks[child]:
                            if grandChild in seen or skipCheck:
                                seen.remove(grandChild)
                                newChildren.append(grandChild)
                    children = newChildren
                for child in children:
                    # Not having the gene implies not mature at some point in the relative future.
                    if "m" + name[child] in evidence and evidence["m" + name[child]] == 1:
                        # This contradiction case should be caught and fixed before we get here.
                        print("CONTRADICTION WHEN CALCULATING " + name[parent] + "!!!")
                        raise Exception("Shouldn't get here.")
                    else:
                        evidence["m" + name[child]] = 0
        elif colour == 2: # Yellow
            evidence[node] = 1 # Yellow cells are coloured.
            depth = 0
            parent = uid[node[1:]]
            seen = set()
            seen.add(parent)
            while depth < colour_min:
                if backwardLinks[parent] != -1:
                    parent = backwardLinks[parent]
                    seen.add(parent)
                    depth += 1
                else:
                    break
            if depth == colour_min:
                # Yellow implies having the gene at some point in the past.
                #if "g" + name[parent] in evidence and evidence["g" + name[parent]] == 0:
                    #print("CONTRADICTION IN EVIDENCE CALCULATION AT " + name[parent] + "!!!")
                    # Resolved by bias towards having the gene and not transferring it to child.
                evidence["g" + name[parent]] = 1
            start = colour_min - depth
            if start <= maturation_max:
                children = [parent]
                for _ in range(start, maturation_max):
                    newChildren = []
                    skipCheck = (len(seen) == 0)
                    for child in children:
                        for grandChild in forwardLinks[child]:
                            if grandChild in seen or skipCheck:
                                seen.remove(grandChild)
                                newChildren.append(grandChild)
                    children = newChildren
                for child in children:
                    # Having the gene implies mature at some point in the relative future.
                    evidence["m" + name[child]] = 1

    if debug >= 1:
        print("Completed basic evidence. Now handling contradictions in model.")

    if progressBar:
        iterator = tqdm.tqdm(model)
    else:
        iterator = model
   
    for node in iterator:
        if node[0] != "g":
            continue

        if (node not in evidence) or (evidence[node] == 0):
            continue

        if debug >= 2:
            print("Working on " + node)
        if progressBar:
            iterator.set_description(desc="Working on " + node)

        # evidence[node] == 1
        for child in list(forwardLinks[uid[node[1:]]]):
            if ("g" + name[child] in evidence) and (evidence["g"+name[child]] == 0):
                forwardLinks[uid[node[1:]]].remove(child)
                model.remove_edge(node, "g"+name[child])
                model.get_cpds("g"+name[child]).delete_evidence(node)
                model.remove_edge("c"+node[1:], "c"+name[child])
                model.get_cpds("c"+name[child]).delete_evidence("c"+node[1:])
                model.remove_edge("m"+node[1:], "m"+name[child])
                model.get_cpds("m"+name[child]).delete_evidence("m"+node[1:])
                # TODO recalculate evidence if we think child was red

    if debug >= 1:
        print("Finished handling contradictions in model. Adding in synchrony evidence.")

    with open(dataFolder + modelName + "_synchrony.pickle", "rb") as f:
        synchrony = pickle.load(f)

    if progressBar:
        iterator = tqdm.tqdm(synchrony)
    else:
        iterator = synchrony

    for cell in iterator:
        if ("g" + name[cell] in evidence) and (evidence["g" + name[cell]] == 1):
            continue

        if debug >= 2:
            print("Working on " + name[cell])
        if progressBar:
            iterator.set_description(desc="Working on " + name[cell])

        evidence["g" + name[cell]] = 1
        children = [cell]
        for _ in range(maturation_max):
            newChildren = []
            for child in children:
                newChildren = newChildren + forwardLinks[child]
            children = newChildren
        for child in children:
            # Having the gene implies mature at some point in the relative future.
            evidence["m" + name[child]] = 1

    if debug >= 1:
        print("Finished including synchrony evidence. Cleaning up edge cases.")

    if progressBar:
        iterator = tqdm.tqdm(model)
    else:
        iterator = model

    for node in iterator:
        if node[0] != "g":
            continue
        
        if node in evidence:
            continue

        if debug >= 2:
            print("Working on " + node)
        if progressBar:
            iterator.set_description(desc="Working on " + node)

        parent = uid[node[1:]]
        trace = []
        while (parent != -1) and ("g" + name[parent] not in evidence):
            trace.append(parent)
            parent = backwardLinks[parent]
        if parent == -1:
            continue
        if evidence["g" + name[parent]] == 1:
            for cell in trace:
                evidence["g" + name[cell]] = 1 # Colour implies we have the gene, but the experiment ended too early.
                children = [cell]
                for _ in range(maturation_max):
                    newChildren = []
                    for child in children:
                        newChildren = newChildren + forwardLinks[child]
                    children = newChildren
                for child in children:
                    # Having the gene implies mature at some point in the relative future.
                    evidence["m" + name[child]] = 1
    
    if debug >= 1:
        print("Finished cleaning up.")
    if save:
        if debug >= 1:
            print("Saving model.")
        with open(modelFolder + modelName + "_model" + modelExtension + "_contradictionsPruned.pickle", "wb") as f:
            pickle.dump(model, f)
        with open(modelFolder + modelName + "_model" + modelExtension + "_contradictionsPruned_evidence.pickle", "wb") as f:
            pickle.dump(evidence, f)
        with open(dataFolder + modelName + "_forwardLinksContradictionsPruned.pickle", "wb") as f:
            pickle.dump(forwardLinks, f)
        #TODO should consider also altering backwardlinks

    return model, evidence, forwardLinks




