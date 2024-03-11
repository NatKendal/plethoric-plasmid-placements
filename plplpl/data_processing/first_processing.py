import csv

# ------------------------------------------------------------------
### Functions for processing the data 
# ------------------------------------------------------------------

# input: gfp and rfp values for one cell
# output: 0 for donor cell, 1 for recipient cell
# only run this if being a transconjugant is already ruled out
def donorOrRecip(gfp, rfp):

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
    
# input: cellId, step number 
# output: link to prior frame (or note the lack thereof)
def getTrackLink(cellId, step):

    return 0

# input: 
def updateLineage(cellId):
    return 0
    
# ------------------------------------------------------------------
### Storage to keep track of data during processing 
# ------------------------------------------------------------------

# we will want a dictionary storing all the cells 
# each will be its own dictionary of values
raw = dict()

# and a dictionary containing cells in each time step
# each a dict with cell id: new unique id
cellsByStep = dict()

# going to manually set the number of steps for now
# create a new list for each time step
numSteps = 225
for i in range(1,numSteps+1):
    cellsByStep[i] = []

# keep track of the row -- starts at 2 because header
# rows will be the new cell id values
currentRow = 2

# will add the processed data to a new dictionary
# indexed by new cell id, each one a dictionary of properties
updated = dict()

# ------------------------------------------------------------------
### Getting the data organized
# ------------------------------------------------------------------

# read in the csv file and store the information in the dictionaries 
with open('data.csv', 'r') as file:
    
    reader = csv.reader(file)

    # since there is a header, we skip line one
    next(reader)

    for row in reader:

        # first add relevant info to the cellsByStep
        step = int(row[0])
        cellId = int(row[2])
        cellsByStep[step].append({id:currentRow})

        # then, add everything to the main dictionary
        raw[currentRow] = {'step': step, 
                           'objectNum': row[1],
                           'cellId': row[2],
                           'lineage': row[3],
                           'divideFlag': row[4],
                           'cellAge': row[5],
                           'growthRate': row[6],
                           'lifetime': row[7],
                           'startLength': row[8],
                           'endLength': row[9],
                           'parent': row[10],
                           'position': row[11],
                           'time': row[12],
                           'width': row[13],
                           'length': row[14],
                           'ends': row[15],
                           'orientation': row[16],
                           'elongationRate': row[17],
                           'avgElongationRate': row[18],
                           'gfp': row[19],
                           'rfp': row[20],
                           'transconjugantFlag': row[21]}

# ------------------------------------------------------------------
### Actual Processing Step
# ------------------------------------------------------------------

# process will be repeated for each cell/row
for cell in raw.keys():

    return 0

'''
    # write a new file - this will overwrite the entire thing
    with open('formatted.csv', 'w', newline='') as new_file:
    
        # makes a writer object for the new file
        writer = csv.writer(new_file)

        # add the new headers
        writer.writerow(['step','id','lineage','growth rate','parent','first','type'])

        # repeat the formatting for each row, writing as we go
        for row in reader:

            # step stays the same
            step = row[0]

            # to get unique id, start with id and append the step 
            id = row[2] + '.' + str(step)

            # lineage stays the same, I think this is the "label"
            lineage = row[3]

            # make all the growths non-negative 
            if row[6] == '':
                growth = 0
            else:
                growth = max(float(row[6]),0)

            # get the parent, and modify the id to match new format
            parent = row[10]

            # first time step has no parent, leave as 0
            if int(step) == 1:
                new_parent = 0
            
            # if there is no parent means the cell did not divide??
            # set to own cell id in prior step
            elif int(parent) == 0:
                # old id + prior time step
                new_parent = row[2] + '.' + str(int(step)-1)
            
            # if there is a parent mentioned, it will be right id
            # so just add on the prior time step
            else:
                new_parent = parent + '.' + str(int(step)-1)

            # get the gfp and rfp and transconjugant flag
            gfp = float(row[19])
            rfp = float(row[20])
            flag = int(row[21])

            # non-lineage version
            elif new_parent in transconjugants:
                first = 0
                transconjugants.append(id)
                cell_type = 't' 

            # to determine the type of cell
            # if flagged as a first, it is a transconjugant
            if flag == 1 and lineage not in lineage:
                firsts += 1
            if flag == 1:
                # counting reported firsts
                transconjugants.append(id)
                lineages.add(lineage)
                cell_type = 't'
            
            elif lineage in lineages:
                transconjugants.append(id)
                cell_type = 't'

            # if not, it must be donor or recipient
            else:
                cell_type = donor_or_recip(gfp,rfp)

            # counting how many children of transconjugants are marked as first of lineage to light up 
            if (flag == 1 and new_parent in transconjugants):
                counter += 1

            # determining if it really is a first transconjugant
            if (flag == 1 and new_parent not in transconjugants):
                first = 1
            else:
                first = 0

            # write them
            writer.writerow([step,id,lineage,growth,new_parent,first,cell_type])

            '''

# the problem is with how im assigning parents -- cells are sometimes not tracked so the parent may not always exist???
# and also with the lineages - they dont update as they split.... have to decide based on parent after all....

if __name__ == "__main__":
    pass #do something
