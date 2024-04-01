# a very first version of the colour delay function

import math

# input: minimum, maximum of range we care about, 
# average and standard deviation of normal distribution 
# output: dictionary of normal distribution values for all time steps in the range 
def calcValues(avg=20, std=3, min=11, max=29):

    vals = dict()

    for i in range(min, max):

        vals[i] = 0.5 * (1 + math.erf( (i - avg) / (std * math.sqrt(2)) ))

    # force maximum value to be 1
    vals[max] = 1

    return vals

# input: time step of cell we are looking at and ancestor, dictionary of values in range
# output: probability the cell we are looking at is active, given that ancestor has gene
def probColour(tCell, tAncestor, dictValues):

    return ( dictValues.get(tCell-tAncestor, 0) )