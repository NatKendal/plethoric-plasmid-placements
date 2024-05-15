import pickle

from plplpl.NoisyOr import NoisyOrBayesianNetwork
from plplpl.NoisyOr import BinaryNoisyOrCPD

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
        elif colour == 1: # Green
            evidence[node] = 0 # Green cells are not coloured.
            if backwardLinks[uid[node[1:]]] == -1:
                # If this is the first frame, also force gene and maturation to match colour.
                evidence["g" + node[1:]] = 0
                evidence["m" + node[1:]] = 0
                continue
            depth = 0
            parent = uid[node[1:]]
            while depth < colour_max:
                if backwardLinks[parent] != -1:
                    parent = backwardLinks[parent]
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
        elif colour == 2: # Yellow
            evidence[node] = 1 # Yellow cells are coloured.
            if backwardLinks[uid[node[1:]]] == -1:
                # If this is the first frame, also force gene and activation to match colour.
                evidence["g" + node[1:]] = 1
                evidence["m" + node[1:]] = 1
                continue
            depth = 0
            parent = uid[node[1:]]
            while depth < colour_min:
                if backwardLinks[parent] != -1:
                    parent = backwardLinks[parent]
                    depth += 1
                else:
                    break
            if depth == colour_min:
                # Yellow implies having the gene at some point in the past.
                #if "g" + name[parent] in evidence and evidence["g" + name[parent]] == 0:
                    #print("CONTRADICTION IN EVIDENCE CALCULATION AT " + name[parent] + "!!!")
                    # Resolved by bias towards having the gene and not transferring it to child.
                evidence["g" + name[parent]] = 1

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
                evidence["m"+name[child]] = 0
                forwardLinks[uid[node[1:]]].remove(child)
                backwardLinks[child] = -1
                
                for parent in list(model.predecessors("g"+name[child])):
                    model.remove_edge(parent, "g"+name[child])
                model.remove_cpds("g"+name[child])
                model.add_cpds(BinaryNoisyOrCPD("g"+name[child], 0))

                for parent in list(model.predecessors("c"+name[child])):
                    model.remove_edge(parent, "c"+name[child])
                model.remove_cpds("c"+name[child])
                model.add_cpds(BinaryNoisyOrCPD("c"+name[child], 0))

                for parent in list(model.predecessors("m"+name[child])):
                    model.remove_edge(parent, "m"+name[child])
                model.remove_cpds("m"+name[child])
                model.add_cpds(BinaryNoisyOrCPD("m"+name[child], 0))

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

    if debug >= 1:
        print("Finished including synchrony evidence. Cleaning up edge cases.")

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
            if evidence[node] == 0:
                if backwardLinks[uid[node[1:]]] != -1 and ("g" + name[backwardLinks[uid[node[1:]]]] not in evidence):
                    parent = backwardLinks[uid[node[1:]]]
                    while parent != -1:
                        if "g" + name[parent] in evidence:
                            print("Errored in pushing back g=0 at " + node)
                            raise AssertionError("Errored in pushing back g=0 at " + node)
                        evidence["g" + name[parent]] = 0
                        parent = backwardLinks[parent]
            elif evidence[node] == 1:
                children = forwardLinks[uid[node[1:]]]
                while children:
                    child = children.pop()
                    if "g" + name[child] in evidence:
                        if evidence["g" + name[child]] == 0:
                            print("Errored in pushing forward g=1 at " + node)
                            raise AssertionError("Errored in pushing forward g=1 at " + node)
                    else:
                        evidence["g" + name[child]] = 1
                        for grandChild in forwardLinks[child]:
                            children.append(grandChild)
            else:
                print("Shouldn't get here.")
                raise AssertionError("Shouldn't get here.")
        else:
            parent = uid[node[1:]]
            seen = []
            while (parent != -1) and ("g" + name[parent] not in evidence):
                seen.append(parent)
                parent = backwardLinks[parent]
            if parent != -1:
                if evidence["g"+name[parent]] == 1:
                    # Push the 1 we found forward and continue.
                    for s in seen:
                        evidence["g"+name[s]] = 1
                    continue
                elif evidence["g"+name[parent]] == 0:
                    seen.append(parent)
            else: 
                print("Something weird happened. All initial cells should be flagged matching their colour, but " + node + " wasn't.")
                raise AssertionError("Something weird happened. All initial cells should be flagged matching their colour, but " + node + " wasn't.")
            children = forwardLinks[uid[node[1:]]]
            reached = None
            while children:
                child = children.pop()
                if "g" + name[child] in evidence:
                    if reached and reached != evidence["g" + name[child]]:
                        print("Found contradiction at " + node + " after it should have been fixed.")
                        raise AssertionError("Found contradiction at " + node + " after it should have been fixed.")
                    else:
                        reached = evidence["g" + name[child]]
                else:
                    seen.append(child)
                    for grandChild in forwardLinks[child]:
                        children.append(grandChild)
            if (reached == 0):
                # 0 pushes zero backward.
                for s in seen:
                    evidence["g"+name[s]] = 0
            elif reached == 1:
                # Found a critical region. Do nothing.
                continue
            elif reached == None:
                for s in seen:
                    evidence["g"+name[s]] = 0
                # This case assumes that cells lost before they change colour never got the plasmid.
                # Technically not guaranteed, but no evidence to the alternative exists.

    if debug >= 1:
        print("Finished cleaning up edge cases. Starting to add maturation evidence.")
    if progressBar:
        iterator = tqdm.tqdm(model)
    else:
        iterator = model

    for node in iterator:
        if node[0] != "g":
            continue

        if node not in evidence:
            continue

        if debug >= 2:
            print("Working on " + node)
        if progressBar:
            iterator.set_description(desc="Working on " + node)

        children = [uid[node[1:]]]
        for _ in range(maturation_min if evidence[node] == 0 else maturation_max):
            newChildren = []
            for child in children:
                newChildren = newChildren + forwardLinks[child]
            children = newChildren
        for child in children:
            evidence["m" + name[child]] = evidence[node]

    if debug >= 1:
        print("Basic maturation evidence added. Cleaning up maturation edge cases now.")
    if progressBar:
        iterator = tqdm.tqdm(model)
    else:
        iterator = model

    for node in iterator:
        if node[0] != "m":
            continue
        if node not in evidence:
            continue

        if debug >= 2:
            print("Working on " + node)
        if progressBar:
            iterator.set_description(desc="Working on " + node)

        if evidence[node] == 0:
            seen = []
            parent = backwardLinks[uid[node[1:]]]
            while (parent != -1) and ("m" + name[parent] not in evidence):
                seen.append(parent)
                parent = backwardLinks[parent]
            if (parent == -1) or (evidence["m" + name[parent]] == 0):
                for s in seen:
                    evidence["m" + name[s]] = 0
            else:
                print("Contradiction when working through maturation edge cases at " + node)
                raise AssertionError("Contradiction when working through maturation edge cases at " + node)
        elif evidence[node] == 1:
            children = forwardLinks[uid[node[1:]]]
            while children:
                child = children.pop()
                if "m" + name[child] in evidence:
                    if evidence["m" + name[child]] == 0:
                        print("Contradiction when working through maturation edge cases at " + node)
                        raise AssertionError("Contradiction when working through maturation edge cases at " + node)
                else:
                    evidence["m" + name[child]] = 0
                    for grandChild in forwardLinks[child]:
                        children.append(grandChild)

    if debug >= 1:
        print("Finished evidence calculation.")
    if save:
        if debug >= 1:
            print("Saving model.")
        with open(modelFolder + modelName + "_model" + modelExtension + "_contradictionsPruned.pickle", "wb") as f:
            pickle.dump(model, f)
        with open(modelFolder + modelName + "_modeldata" + modelExtension + "_evidence.pickle", "wb") as f:
            pickle.dump(evidence, f)
        with open(modelFolder + modelName + "_modeldata" + modelExtension + "_forwardLinksPostEvidence.pickle", "wb") as f:
            pickle.dump(forwardLinks, f)
        with open(modelFolder + modelName + "_modeldata" + modelExtension + "_backwardLinksPostEvidence.pickle", "wb") as f:
            pickle.dump(backwardLinks, f)

    return model, evidence, forwardLinks, backwardLinks




