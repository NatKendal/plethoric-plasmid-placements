# a very first version of the colour delay function

import math

# input: minimum, maximum of range we care about, scaling of exponential distribution
# output: dictionary of exponential distribution values for all time steps in the range 

# IMPORTANT: this implementation gives a probability of 0 at the minimum -- may need to adjust this

def calcValues(scale=1/3, min=3, max=18):

    vals = dict()

    for i in range(min, max):

        vals[i] = 1 - math.exp(-1 * scale * (i-min))

    # force maximum value to be 1
    vals[max] = 1

    return vals

# input: time step of cell we are looking at and ancestor, dictionary of values in range
# output: probability the cell we are looking at is yellow, given that ancestor has gene
def probColour(tCell, tAncestor, dictValues):

    return ( dictValues.get(tCell-tAncestor, 0) )