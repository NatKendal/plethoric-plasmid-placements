import csv
import math 

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
    # each cell is a list [cell, unique]
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
                cellsByStep[step] = [(cellId,currentRow)]

            # then, add everything to the main dictionary
            raw[currentRow] = {
                            'step': int(row[0]), 
                            #'objectNum': row[1],
                            'cellId': int(row[2]),
                            'lineage': int(row[3]),
                            #'divideFlag': row[4],
                            #'cellAge': row[5],
                            #'growthRate': row[6],
                            #'lifetime': row[7],
                            #'startLength': row[8],
                            #'endLength': row[9],
                            'parentCellId': int(row[10]),
                            'position': row[11],
                            #'time': row[12],
                            'width': float(row[13]),
                            'length': float(row[14]),
                            'ends': row[15],
                            'orientation': row[16],
                            'elongationRate': float(row[17]),
                            #'avgElongationRate': row[18],
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
### Functions for processing the data 
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

# input: all data
# output: dict of all transconjugants, list of first transconjugants
# -1 for not conjugant, 2 for conjugant (not reusing 0/1 to avoid errors)
def getConj(raw):

    # stores all seen transconjugant cell id
    known = []

    # stores uid of all first transconjugants
    firsts = []

    # dictionary we want
    dictConj = dict()

    # sort the keys, just in case
    orderedCells = list(raw.keys())
    orderedCells.sort()

    # check each row to see if it is flagged
    # we sort it to make sure iwe go forward in time
    for cell in orderedCells:

        # checks
        status = -1

        # check if its cell Id is already there
        # in this case it is not a first, nor does it need to be added
        if raw[cell]['cellId'] in known:
            status = 2

        # check if parent cell Id is there
        # in this case it divided but was not a first
        # does need to be added
        elif raw[cell]['parentCellId'] in known:
            known.append(raw[cell]['cellId'])
            status = 2

        # check if flagged, and not yet known
        # this would be a first, add uid to that list
        elif raw[cell]['flag'] == 1:
            known.append(raw[cell]['cellId'])
            firsts.append(cell)
            status = 2

        # add the dictionary entry
        dictConj[cell] = status

    return dictConj, firsts
    
# ------------------------------------------------------------------
### Actual Dictionary Outputs
# ------------------------------------------------------------------

# input: raw data
# output: dict of uid to time step
# no changes to data
def dictSteps(raw):

    step = dict()

    for cell in raw.keys():

        step[cell] = raw[cell]['step']

    return step

# input: raw data
# output: dict of uid to cell id
# no changes to data
def dictCellId(raw):

    cellId = dict()

    for cell in raw.keys():

        cellId[cell] = raw[cell]['cellId']

    return cellId

# input: raw data
# output: dict of uid to parent Cell Id
# no changes to data
def dictParent(raw):

    parent = dict()

    for cell in raw.keys():

        parent[cell] = raw[cell]['parentCellId']

    return parent

# input: raw data
# output: dict of uid to positions
# no changes to data
def dictPosition(raw):

    position = dict()

    for cell in raw.keys():

        position[cell] = raw[cell]['position']

    return position

# input: raw data
# output: dict of uid to widths
# no changes to data
def dictWidth(raw):

    width = dict()

    for cell in raw.keys():

        width[cell] = raw[cell]['width']

    return width

# input: raw data
# output: dict of uid to length
# no changes to data
def dictLength(raw):

    length = dict()

    for cell in raw.keys():

        length[cell] = raw[cell]['length']

    return length

# input: raw data
# output: dict of uid to ends
# no changes to data
def dictEnds(raw):

    ends = dict()

    for cell in raw.keys():

        ends[cell] = raw[cell]['ends']

    return ends


# input: raw data
# output: dict of uid to ends
# no changes to data
def dictOrientation(raw):

    orientation = dict()

    for cell in raw.keys():

        orientation[cell] = raw[cell]['orientation']

    return orientation

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

# input: raw data, dictionary of transconjugants
# output: dictionary mapping each uid to cell colour
# 0 for donor cell, 1 for recipient cell, 2 is transconjugant
def dictColours(raw, dictConj):

    colours = dict()

    for cell in raw.keys():

        # check if it is a transconjugant first
        if dictConj[cell] == 2:
            colours[cell] = 2

        # otherwise check if red or green
        else:
            colours[cell] = redOrGreen(raw[cell]['rfp'], raw[cell]['gfp'])

    return colours

# input: raw data, cellsByStep
# output: dict of uid to uid of past track/parent link
# no distinction between parent/self 
def dictBackwardLinks(raw,cellsByStep):

    backwardLinks = dict()

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
                        found = True
                        break

            # if we also did not find a parent, error, set to -1
            if not found:
                backwardLinks[cell] = -1

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
        found = False
        for item in cellsByStep[step-1]:
            if item[0] == raw[cell]['cellId']:
                forwardLinks[item[1]].append(cell)
                break

            elif item[0] == raw[cell]['parentCellId']:
                forwardLinks[item[1]].append(cell)
                break

    return forwardLinks

# input: raw data, cellsByStep
# output: dict of uid to all potential neighbours 
def dictNeighbours(raw, cellsByStep, maxDistance=25):

    # initialize to no neighbours
    neighbours = {uid: [] for uid in raw.keys()}

    # we find neighbors for each time step
    for step in cellsByStep.keys():

        # list of uid of cells not yet checked in step
        # second element in each item is the cell uid
        toCheck = [item[1] for item in cellsByStep[step]]

        # double loop - compare each cell against every other
        for cell1 in toCheck:

            # no need to compare against self / to ever check it later
            toCheck.remove(cell1)

            # check against all remaining cells
            for cell2 in toCheck:

                distance = math.dist(raw[cell1]['position'], raw[cell2]['position'])

                # potential neighbours must be within 5 cell lengths (5 microns each)
                # this is a quarter of the trap... maybe we should start smaller...
                if distance <= maxDistance:

                    # symmetric relation so add to both
                    neighbours[cell1].append(cell2)
                    neighbours[cell2].append(cell1)

    return neighbours

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

# input: raw data
# output: dict of uid to lineages
# only updating id to uid
def dictLineage(raw):

    lineage = dict()

    ''' TO DO'''

    return lineage

