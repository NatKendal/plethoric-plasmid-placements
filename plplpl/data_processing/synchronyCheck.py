# for checking if pairs of siblings / sets of cousins light up at the same time

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

