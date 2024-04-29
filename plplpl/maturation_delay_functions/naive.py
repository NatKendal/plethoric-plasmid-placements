# a very first version of the maturation delay function

import math

# input: minimum, maximum of range we care about, 
# average and standard deviation of normal distribution
# output: dictionary of normal distribution values for all time steps in the range 
def calcValues(avg=16, std=5, min=2, max=31):

    vals = dict()

    for i in range(min, max):

        vals[i] = 0.5 * (1 + math.erf( (i - avg) / (std * math.sqrt(2)) ))

    # force maximum value to be 1
    vals[max+1] = 1

    return vals

# input: time step of cell we are looking at and ancestor, dictionary of values in range
# output: probability the cell we are looking at is mature, given that ancestor has gene
def probMature(tCell, tAncestor, dictValues):

    return ( dictValues.get(tCell-tAncestor, 0) )