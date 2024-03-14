import csv

# ------------------------------------------------------------------
### Functions for processing the data 
# ------------------------------------------------------------------

# input: gfp and rfp values for one cell
# output: 0 for donor cell, 1 for recipient cell
# only run this if being a transconjugant is already ruled out
def redOrGreen(rfp,gfp):

    # fluorescence thresholds for red/green
    recipThreshold=2280/16383
    donorThreshold=1700/16383

    # check if it is above each threshold and assign a boolean value
    recip = (gfp >= recipThreshold)
    donor = (rfp >= donorThreshold)

    # if it is above only one threshold, return that one
    if recip and not donor:
        return 1
    
    elif donor and not recip:
        return 0
    
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
    known = set()

    # stores id of all first transconjugants
    firsts = []

    # dictionary we want
    dictConj = dict()

    # check each row to see if it is flagged
    # we sort it to make sure iwe go forward in time
    for cell in raw.keys().sort():

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
            known.add(raw[cell]['cellId'])
            status = 2

        # check if flagged, and not yet known
        # this would be a first, add uid to that list
        elif raw[cell]['flag'] == 1:
            known.add(raw[cell]['cellId'])
            firsts.append(cell)
            status = 2


        # add the dictionary entry
        dictConj[cell] = status

    return dictConj, firsts
    
# input: cellId, step number 
# output: link to prior frame (or note the lack thereof)
def getTrackLink(cellId, step):
    return 0

# input: 
def updateLineage(cellId):
    return 0
    
# ------------------------------------------------------------------
### Getting the data organized
# ------------------------------------------------------------------

# input: file name
# output: dictionary storing all cells and info, dictionary of time step + uid/cell id
# makes python dictionaries to store the data from the csv while processing
def readFile(file):

    # we will want a dictionary storing all the cells 
    # each will be its own dictionary of values
    raw = dict()

    # and a dictionary containing cells in each time step
    # each a dict with unique id: cell id
    cellsByStep = dict()

    # going to manually set the number of steps for now - can probably get this info?
    # create a new list for each time step
    numSteps = 225
    for i in range(1,numSteps+1):
        cellsByStep[i] = []

    # keep track of the row -- starts at 2 because header
    # rows will be the new cell id values
    currentRow = 2

    # read in the csv file and store the information in the dictionaries 
    with open('data.csv', 'r') as file:
        
        reader = csv.reader(file)

        # since there is a header, we skip line one
        next(reader)

        for row in reader:

            # first add relevant info to the cellsByStep
            step = int(row[0])
            cellId = int(row[2])
            cellsByStep[step].append({currentRow: cellId})

            # then, add everything to the main dictionary
            raw[currentRow] = {'step': step, 
                            #'objectNum': row[1],
                            'cellId': row[2],
                            'lineage': row[3],
                            #'divideFlag': row[4],
                            #'cellAge': row[5],
                            #'growthRate': row[6],
                            #'lifetime': row[7],
                            #'startLength': row[8],
                            #'endLength': row[9],
                            'parentCellId': row[10],
                            'position': row[11],
                            #'time': row[12],
                            'width': row[13],
                            'length': row[14],
                            'ends': row[15],
                            'orientation': row[16],
                            'elongationRate': row[17],
                            #'avgElongationRate': row[18],
                            'gfp': row[19],
                            'rfp': row[20],
                            'flag': row[21]}
            
            currentRow += 1
            
    return raw, cellsByStep

# ------------------------------------------------------------------
### Actual Processing Steps 
# ------------------------------------------------------------------

# input: raw data, dictionary of transconjugants
# output: dictionary mapping each uid to cell colour
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

