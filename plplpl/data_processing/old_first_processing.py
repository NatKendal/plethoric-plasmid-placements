import csv
import math
import os
import pickle
import sys

MAX_DISTANCE = 10
# ------------------------------------------------------------------
### Getting the data organized
# ------------------------------------------------------------------

# input: file name
# output: dictionary storing all cells and info, dictionary of time step + uid/cell id
# makes python dictionaries to store the data from the csv while processing
def readFile(filename):
    # we will want a dictionary storing all the cells 
    # each will be its own dictionary of values
    raw = dict()

    # and a dictionary containing cells in each time step
    # each cell is a tuple (cell, unique)
    cellsByStep = dict()

    # keep track of the row -- starts at 2 because header
    # rows will be the new cell id values
    currentRow = 2

    #read in the csv file and store the information in the dictionaries 
    with open(filename, 'r') as file:
        
        reader = csv.reader(file, delimiter=',')

        # since there is a header, we skip line one
        next(reader)

        for row in reader:

            # first add relevant info to the cellsByStep
            step = int(row[0])
            cellId = int(row[2])

            if step in cellsByStep.keys():
                cellsByStep[step].append((cellId, currentRow))
            else:
                cellsByStep[step] = []
                cellsByStep[step].append((cellId,currentRow))

            # then, add everything to the main dictionary
            raw[currentRow] = {
                            'step': int(row[0]), 
                            #'objectNum': row[1], 
                            'cellId': int(row[2]),
                            'lineage': int(float(row[3])),
                            #'divideFlag': row[4],
                            'cellAge': row[5],
                            'growthRate': row[6],
                            'lifetime': row[7],
                            'startLength': row[8],
                            'endLength': row[9],
                            'parentCellId': int(row[10]),
                            'position': row[11],
                            #'time': row[12],
                            'width': float(row[13]),
                            'length': float(row[14]),
                            'ends': row[15],
                            'orientation': row[16],
                            'elongationRate': float(row[17]),
                            'avgElongationRate': row[18],
                            'gfp': float(row[19]),
                            'rfp': float(row[20]),
                            'flag': int(row[21])}

            # make things into ints/floats
            pos = raw[currentRow]['position'][1:-1].split(', ')
            raw[currentRow]['position'] = [float(pos[0]), float(pos[1])]

            ori = raw[currentRow]['orientation'][1:-1].split(', ')
            raw[currentRow]['orientation'] = [float(ori[0]), float(ori[1])]

            ends = raw[currentRow]['ends'][2:-2].split(', ')
            raw[currentRow]['ends'] = [[float(ends[0]), float(ends[1][:-1])], [float(ends[2][1:]), float(ends[3])]]

            currentRow += 1
            
    return raw, cellsByStep

# ------------------------------------------------------------------
### Functions for colour
# ------------------------------------------------------------------

# input: gfp and rfp values for one cell
# output: 0 for donor cell, 1 for recipient cell
# only run this if being a transconjugant is already ruled out
def redOrGreen(rfp, gfp):

    # fluorescence thresholds for red/green
    donorThreshold=1700/16383
    recipThreshold=2280/16383

    # check if it is above each threshold and assign a boolean value
    donor = (rfp >= donorThreshold)
    recip = (gfp >= recipThreshold)

    # if it is above only one threshold, return that one
    if donor and not recip:
        return 0
    
    elif recip and not donor:
        return 1
    
    # otherwise, normalize and take the larger one
    # bias is towards assuming it is a recipient
    elif gfp/recipThreshold >= rfp/donorThreshold:
        return 1
        
    else:
        return 0

# input: byStep, raw
# output: dictionary mapping uid to cell colour, list of transconjugants in frame 1
# 0 = red = donor, 1 = green = recipient, 2 = yellow = transconjugant
def getFirstFrameColour(byStep,raw):

    firstColours = dict()

    firstTransconjugants = []

    # first step 
    for cell in byStep[1]:

        if raw[cell]['flag'] == 1:
            firstColours[cell] = 2
            firstTransconjugants.append(cell)

        else:
            firstColours[cell] = redOrGreen(rfp=raw[cell]['rfp'],gfp=raw[cell]['gfp'])

    return firstColours, firstTransconjugants
    
# ------------------------------------------------------------------
### Actual Dictionary Outputs
# ------------------------------------------------------------------

# input: raw data
# output: dict of uid to time step
def dictSteps(raw):
    return {cell:raw[cell]["step"] for cell in raw}

# input: raw data
# output: dict of uid to cell id
def dictCellId(raw):
    return {cell:raw[cell]["cellId"] for cell in raw}

# input: raw data
# output: dict of uid to parent Cell Id
def dictParent(raw):
    return {cell:raw[cell]["parentCellId"] for cell in raw}

# input: raw data
# output: dict of uid to positions
def dictPosition(raw):
    return {cell:raw[cell]["position"] for cell in raw}

# input: raw data
# output: dict of uid to widths
def dictWidth(raw):
    return {cell:raw[cell]["width"] for cell in raw}

# input: raw data
# output: dict of uid to length
def dictLength(raw):
    return {cell:raw[cell]["length"] for cell in raw}

# input: raw data
# output: dict of uid to ends
def dictEnds(raw):
    return {cell:raw[cell]["ends"] for cell in raw}

# input: raw data
# output: dict of uid to ends
def dictOrientation(raw):
    return {cell:raw[cell]["orientation"] for cell in raw}

# input: raw data
# output: dict of uid to growth rates
# -1 for negative values, otherwise exactly as reported
def dictGrowth(raw):

    growth = dict()

    for cell in raw.keys():
        if raw[cell]['elongationRate'] < 0:
            growth[cell] = -1  
        else:
            growth[cell] = raw[cell]['elongationRate']

    return growth

# input: raw data, cellsByStep
# output: dict of uid to uid of past track/parent link
# no distinction between parent/self 
def dictBackwardLinks(raw,cellsByStep):

    # initalize to no backward link
    backwardLinks = {uid: -1 for uid in raw.keys()}

    for cell in raw.keys():

        # check the current time step
        step = raw[cell]['step']

        # if it is in the first time step, no backwards link
        if step == 1:
            backwardLinks[cell] = -1

        else:
            # first check if the cell id exists in the previous time frame
            found = False
            for item in cellsByStep[step-1]:
                if item[0] == raw[cell]['cellId']:
                    backwardLinks[cell] = item[1]
                    found = True
                    break
            
            # if we did not find its prior self, check for parents
            if not found:
                for item in cellsByStep[step-1]:
                    if item[0] == raw[cell]['parentCellId']:
                        backwardLinks[cell] = item[1]
                        break

    return backwardLinks

# input: raw data, cellsByStep
# output: dict of uid to uid of future track/child link
# no distinction between parent/self 
def dictForwardLinks(raw, cellsByStep):

    # initialize to no forward link
    forwardLinks = {uid: [] for uid in raw.keys()}

    for cell in raw.keys():

        # check the current time step
        step = raw[cell]['step']

        # if it is in the first time step, nothing to check 
        if step == 1:
            continue

        # check if the cell id exists in the previous time frame
        # or if the parent id exists
        # if it does, current cell is the forward link/child of that cell
        for item in cellsByStep[step-1]:
            if item[0] == raw[cell]['cellId']:
                forwardLinks[item[1]].append(cell)
                break

            elif item[0] == raw[cell]['parentCellId']:
                forwardLinks[item[1]].append(cell)
                break

    return forwardLinks

# input: dictSteps
# output: dictionary of time step to list of uid of cells in that step
def dictByStep(dictSteps):

    byStep = dict()

    for cell in dictSteps.keys():

        step = dictSteps[cell]

        # if the step is already in the dictionary, add the cell to the list
        if step in byStep.keys():
            byStep[step].append(cell)

        # otherwise, make the dictionary entry
        else:
            byStep[step] = [cell]

    return byStep

# input: raw data, byStep
# output: dict of uid to all potential neighbours 
def dictNeighbours(raw, byStep, maxDistance=MAX_DISTANCE):

    # initialize to no neighbours
    neighbours = {uid: [] for uid in raw.keys()}

    # we find neighbors for each time step
    for step in byStep.keys():

        # list of uid of cells not yet checked in step
        # second element in each item is the cell uid
        toCheck = [item for item in byStep[step]]

        # double loop - compare each cell against every other
        for cell1 in toCheck:

            # no need to compare against self / to ever check it later
            toCheck.remove(cell1)

            # check against all remaining cells
            for cell2 in toCheck:

                distance = math.dist(raw[cell1]['position'], raw[cell2]['position'])

                # potential neighbours must be within maxDistance
                # this is a quarter of the trap... maybe we should start smaller...
                if distance <= maxDistance:

                    # symmetric relation so add to both
                    neighbours[cell1].append(cell2)
                    neighbours[cell2].append(cell1)

    return neighbours

# input: raw data, backward links, cells in each step; dictionary of colours for first step
# output: dictionary mapping each uid to cell colour, list of first transconjugants
# 0 for donor cell, 1 for recipient cell, 2 is transconjugant
def dictColours(raw, backwardLinks, byStep, firstColours, firstTransconjugants):

    colours = firstColours
    firsts = firstTransconjugants

    # go forwards by step so that parents are always assigned colours first
    steps = list(byStep.keys())
    steps.sort()

    for step in steps:
        for cell in byStep[step]:

            previous = backwardLinks[cell]

            # see if the previous cell (self or parent) has a colour assigned
            if previous in colours.keys():

                # if it is red (0) or yellow (2), so is the cell
                if colours[previous] == 0:
                    colours[cell] = 0
                elif colours[previous] == 2:
                    colours[cell] = 2

                # if green (1) check first if it is flagged
                # if flagged, colour is 2 and is a first
                # otherwise stays green
                else:
                    if raw[cell]['flag'] == 1:
                        colours[cell] = 2
                        firsts.append(cell)
                    else:
                        colours[cell] = 1

            # if no previous cell, check the colour, first looking for flags
            # if transconjugant, necessarily the first in its lineage
            else:
                if raw[cell]['flag'] == 1:
                    colours[cell] = 2
                    firsts.append(cell)

                else:
                    colours[cell] = redOrGreen(rfp=raw[cell]['rfp'],gfp=raw[cell]['gfp'])

    return colours, firsts


# input: raw data, dictBackwardLinks
# output: dict of uid to lineages (uid of first cell in lineage)
def dictLineage(raw, dictBackwardLinks):

    lineage = {uid: -1 for uid in raw.keys()}

    for cell in raw.keys():

        current = cell

        # step back through parents until either the parent has a lineage
        # or you run out of parents to check
        while lineage[cell] == -1:

            if lineage[current] != -1:
                lineage[cell] = lineage[current]

            elif dictBackwardLinks[current] == -1:
                lineage[cell] = current

            else:
                current = dictBackwardLinks[current]

    return lineage

def dictCellAge(raw):
    return {cell:raw[cell]["cellAge"] for cell in raw}

def dictGrowthRate(raw):
    return {cell:raw[cell]["growthRate"] for cell in raw}

def dictLifetime(raw):
    return {cell:raw[cell]["lifetime"] for cell in raw}

def dictStartLength(raw):
    return {cell:raw[cell]["startLength"] for cell in raw}

def dictEndLength(raw):
    return {cell:raw[cell]["endLength"] for cell in raw}

def dictAvgElongationRate(raw):
    return {cell:raw[cell]["avgElongationRate"] for cell in raw}

def dictHumanFriendlyName(raw):
    return {cell:str(raw[cell]["cellId"])+"_"+str(raw[cell]["step"]) for cell in raw}

def saveAll(rawCSVfilename, saveDirectory, modelname):
    raw, cellsByStep = readFile(rawCSVfilename)
    backwardLinks = dictBackwardLinks(raw, cellsByStep)
    steps = dictSteps(raw)
    byStep = dictByStep(steps)
    firstColours, firstTransconjugants = getFirstFrameColour(byStep, raw)
    colours, firsts = dictColours(raw,backwardLinks,byStep,firstColours,firstTransconjugants)
    humanFriendlyName = dictHumanFriendlyName(raw)
    uid = {name: uid for uid, name in humanFriendlyName.items()}
    os.makedirs(saveDirectory, exist_ok=True)

    if not os.path.isfile(saveDirectory + modelname + "_firstFrameColours.pickle"):
        with open(saveDirectory + modelname + "_firstFrameColours.pickle", "wb") as f:
            pickle.dump(firstColours, f)

    if not os.path.isfile(saveDirectory + modelname + "_firstFrameTransconjugants.pickle"):
        with open(saveDirectory + modelname + "_firstFrameTransconjugants.pickle", "wb") as f:
            pickle.dump(firstTransconjugants, f)

    if not os.path.isfile(saveDirectory + modelname + "_cellId.pickle"):
        with open(saveDirectory + modelname + "_cellId.pickle", "wb") as f:
            pickle.dump(dictCellId(raw), f)

    if not os.path.isfile(saveDirectory + modelname + "_parent.pickle"):
        with open(saveDirectory + modelname + "_parent.pickle", "wb") as f:
            pickle.dump(dictParent(raw), f)

    if not os.path.isfile(saveDirectory + modelname + "_position.pickle"):
        with open(saveDirectory + modelname + "_position.pickle", "wb") as f:
            pickle.dump(dictPosition(raw), f)

    if not os.path.isfile(saveDirectory + modelname + "_width.pickle"):
        with open(saveDirectory + modelname + "_width.pickle", "wb") as f:
            pickle.dump(dictWidth(raw), f)

    if not os.path.isfile(saveDirectory + modelname + "_length.pickle"):
        with open(saveDirectory + modelname + "_length.pickle", "wb") as f:
            pickle.dump(dictLength(raw), f)

    if not os.path.isfile(saveDirectory + modelname + "_ends.pickle"):
        with open(saveDirectory + modelname + "_ends.pickle", "wb") as f:
            pickle.dump(dictEnds(raw), f)

    if not os.path.isfile(saveDirectory + modelname + "_orientaton.pickle"):
        with open(saveDirectory + modelname + "_orientation.pickle", "wb") as f:
            pickle.dump(dictOrientation(raw), f)

    if not os.path.isfile(saveDirectory + modelname + "_growth.pickle"):
        with open(saveDirectory + modelname + "_growth.pickle", "wb") as f:
            pickle.dump(dictGrowth(raw), f)

    if not os.path.isfile(saveDirectory + modelname + "_colours.pickle"):
        with open(saveDirectory + modelname + "_colours.pickle", "wb") as f:
            pickle.dump(colours, f)

    if not os.path.isfile(saveDirectory + modelname + "_forwardLinks.pickle"):
        with open(saveDirectory + modelname + "_forwardLinks.pickle", "wb") as f:
            pickle.dump(dictForwardLinks(raw, cellsByStep), f)

    if not os.path.isfile(saveDirectory + modelname + "_neighbours.pickle"):
        with open(saveDirectory + modelname + "_neighbours.pickle", "wb") as f:
            pickle.dump(dictNeighbours(raw, byStep), f)

    if not os.path.isfile(saveDirectory + modelname + "_byStep.pickle"):
        with open(saveDirectory + modelname + "_byStep.pickle", "wb") as f:
            pickle.dump(byStep, f)

    if not os.path.isfile(saveDirectory + modelname + "_lineage.pickle"):
        with open(saveDirectory + modelname + "_lineage.pickle", "wb") as f:
            pickle.dump(dictLineage(raw, backwardLinks), f)

    if not os.path.isfile(saveDirectory + modelname + "_cellAge.pickle"):
        with open(saveDirectory + modelname + "_cellAge.pickle", "wb") as f:
            pickle.dump(dictCellAge(raw), f)

    if not os.path.isfile(saveDirectory + modelname + "_growthRate.pickle"):
        with open(saveDirectory + modelname + "_growthRate.pickle", "wb") as f:
            pickle.dump(dictGrowthRate(raw), f)

    if not os.path.isfile(saveDirectory + modelname + "_lifetime.pickle"):
        with open(saveDirectory + modelname + "_lifetime.pickle", "wb") as f:
            pickle.dump(dictLifetime(raw), f)

    if not os.path.isfile(saveDirectory + modelname + "_startLength.pickle"):
        with open(saveDirectory + modelname + "_startLength.pickle", "wb") as f:
            pickle.dump(dictStartLength(raw), f)

    if not os.path.isfile(saveDirectory + modelname + "_endLength.pickle"):
        with open(saveDirectory + modelname + "_endLength.pickle", "wb") as f:
            pickle.dump(dictEndLength(raw), f)

    if not os.path.isfile(saveDirectory + modelname + "_avgElongationRate.pickle"):
        with open(saveDirectory + modelname + "_avgElongationRate.pickle", "wb") as f:
            pickle.dump(dictAvgElongationRate(raw), f)

    if not os.path.isfile(saveDirectory + modelname + "_humanFriendlyName.pickle"):
        with open(saveDirectory + modelname + "_humanFriendlyName.pickle", "wb") as f:
            pickle.dump(humanFriendlyName, f)

    if not os.path.isfile(saveDirectory + modelname + "_humanFriendlyNameLookup.pickle"):
        with open(saveDirectory + modelname + "_humanFriendlyNameLookup.pickle", "wb") as f:
            pickle.dump(uid, f)

    if not os.path.isfile(saveDirectory + modelname + "_raw.pickle"):
        with open(saveDirectory + modelname + "_raw.pickle", "wb") as f:
            pickle.dump(raw, f)

    if not os.path.isfile(saveDirectory + modelname + "_backwardLinks.pickle"):
        with open(saveDirectory + modelname + "_backwardLinks.pickle", "wb") as f:
            pickle.dump(backwardLinks, f)

    if not os.path.isfile(saveDirectory + modelname + "_step.pickle"):
        with open(saveDirectory + modelname + "_step.pickle", "wb") as f:
            pickle.dump(steps, f)

    if not os.path.isfile(saveDirectory + modelname + "_firsts.pickle"):
        with open(saveDirectory + modelname + "_firsts.pickle", "wb") as f:
            pickle.dump(firsts, f)

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Expects 3 arguments: [csvfile] [folder to save into] [modelname]")
    else:
        saveAll(sys.argv[1], sys.argv[2], sys.argv[3])

