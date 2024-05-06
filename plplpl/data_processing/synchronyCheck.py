# for checking if pairs of siblings / sets of cousins light up at the same time

import pickle

TOLERANCE = 2

# input: dictionary of cells per step, colours, forward links; list of first conjugation events
# output: dictionary of uid --> list representing how many descendants became transconjugants, and when
# [] = yellow and not first, -1 entry = descendant that never lit up
# non-negative entry = descendant lit up that many steps in future
def getSynchrony(byStep,colours,forwardLinks,firsts):

    # create the dictionary
    synchrony = dict()

    # sort the keys of byStep so we start from the end
    steps = list(byStep.keys())
    steps.sort()
    steps.reverse()

    # repeat for each time step, then for each cell
    for step in steps:
        for cell in byStep[step]:

            # if the cell is a donor, skip everything
            colour = colours[cell]

            if colour == 0:
                continue

            # if it has no forward links, initiate the list
            if forwardLinks[cell] == []:

                # recipient cells tagged with -1
                if colour == 1:
                    synchrony[cell] = [-1]

                # transconjugant - [0] for first, [] otherwise
                else:
                    if cell in firsts:
                        synchrony[cell] = [0]
                    else:
                        synchrony[cell] = []

            # if it does have forward links, combine future lists
            else:
                entry = []
                for child in forwardLinks[cell]:
                    entry.extend(synchrony[child])
                
                # if it is still empty, only yellow children
                # must check if first
                if entry == []:
                    if cell in firsts:
                        synchrony[cell] = [0]
                    else:
                        synchrony[cell] = []

                # if there is a -1, has an always green descendant
                # so we will always ignore it and can set to -1
                elif -1 in entry:
                    synchrony[cell] = [-1]

                # otherwise, update number of steps to each conjugation event
                else:
                    synchrony[cell] = [number+1 for number in entry]

    return synchrony


# input: list of all cells; dictionary of colours, forward links, synchrony;
# tolerance for disparity in siblings/cousins
# output: set of uid of green (recipient) cells which divide in next step and must have the gene
# may be more than one per lineage
def getFirstCertain(allCells,colours,forwardLinks,synchrony,tol=TOLERANCE):

    # create the set
    certain = set()

    # order is irrelevant, so go through all cells
    for cell in allCells:

        # red cells will not be in the dictionary, so skip them
        if colours[cell] == 0:
            continue

        # cells with <= 1 item in their list won't give certainty, skip
        if len(synchrony[cell]) <= 1:
            continue

        # only can be certain for the cell right before division
        # so we only want to check ones with at least two forward links
        if len(forwardLinks[cell]) > 1:

            # check the min and max, difference must be no more than tolerance
            min_val = min(synchrony[cell])
            max_val = max(synchrony[cell])
            diff = max_val - min_val

            if diff <= tol:
                certain.add(cell)

    return certain
            
# input: uid of certain cells; dictionary of colours and forward links
# output: list of uid of all green cells guaranteed to have the gene
def getAllGreenCertain(certain,colours,forwardLinks):

    allCertain = certain.copy()

    # for each certain cell, track down its lineage until the cell becomes yellow 
    # or splits into two cells (that parent will already be in certain)
    # or the lineage ends 
    for cell in certain:
        for child in forwardLinks[cell]:

            current = child
            nextFlag = True

            while nextFlag:
                # first check the colour - green added to allCertain, otherwise stop
                if colours[current] == 1:
                    allCertain.add(child)
                else:
                    nextFlag = False

                # if it has no forward links or too many, stop
                # otherwise take the one link as the new current cell
                if len(forwardLinks[current]) != 1:
                    nextFlag = False
                else: 
                    current = forwardLinks[current][0]

    return list(allCertain)

# input: uid of certain cells; dictionary byStep and backward links
# output: list of all cells guaranteed to have the gene
def getAllCertain(certain,byStep,backwardLinks):

    allCertain = certain

    # sort the keys of byStep so we start from the top
    steps = list(byStep.keys())
    steps.sort()

    # for each cell, if its backward link is in allCertain, add it
    for step in steps:
        for cell in byStep[step]:
            if backwardLinks[cell] in allCertain:
                allCertain.add(cell)

    return list(allCertain)

def saveSynchronyCertainty(dataFolder, modelName, save=True):

    # get all the relevant dictionaries
    with open(dataFolder + modelName + "_byStep.pickle", "rb") as f:
        byStep = pickle.load(f)

    with open(dataFolder + modelName + "_colours.pickle", "rb") as f:
        colours = pickle.load(f)

    with open(dataFolder + modelName + "_backwardLinks.pickle", "rb") as f:
        backwardLinks = pickle.load(f)

    with open(dataFolder + modelName + "_forwardLinks.pickle", "rb") as f:
        forwardLinks = pickle.load(f)

    with open(dataFolder + modelName + "_firsts.pickle", "rb") as f:
        firsts = pickle.load(f)

    # Save the list of unique cell ids for later.
    cells = list(colours.keys())

    synchrony = getSynchrony(byStep,colours,forwardLinks,firsts)

    certain = getFirstCertain(cells,colours,forwardLinks,synchrony,tol=TOLERANCE)

    allGreenCertain = getAllGreenCertain(certain,colours,forwardLinks)

    if save:
        with open(dataFolder + modelName + "_synchrony.pickle", "wb") as f:
            pickle.dump(allGreenCertain, f)

    return allGreenCertain




   
    
