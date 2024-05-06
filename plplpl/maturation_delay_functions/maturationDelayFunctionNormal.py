# specify the cdf, get the edge weights for maturation

import math
import numpy as np
import pathlib

from plplpl.base_functions import BaseDelayFunction

# input: average and standard deviation in minutes, step length
# output: dictionary of normal distribution values for all time steps in 3x standard deviation
def calcCdfValues(avg=80, std=20, step_length = 5):

    # convert minutes to steps
    # decrease min by one step so non-zero probability at the actual minimum
    avg = math.floor(avg/step_length) 
    std = math.ceil(std/step_length)

    min_val = avg - 3 * std 
    max_val = avg + 3 * std 

    assert min_val > 0, "minimum maturation time is less than 0"

    cdf_vals = dict()

    for i in range(min_val, max_val+1):

        cdf_vals[i] = 0.5 * (1 + math.erf( (i - avg) / (std * math.sqrt(2))))

    return cdf_vals, min_val, max_val

# input: dictionary of cdf values
# output: dictionary of edge weights
def calcMaturationWeights(cdf_vals):

    # get minimum and maximum
    min_val = min(cdf_vals.keys())
    max_val = max(cdf_vals.keys())

    # create dictionary for edge weights and enter first and last
    edge_vals = dict()

    # formula for weight_i is 1 - (1 - cdf(i))/(prod_k=1^i-1 (1-weight_k))
    for i in range(min_val, max_val+1):

        num = 1 - cdf_vals[i]
        denom = np.prod([1 - edge_vals[j] for j in range(min_val,i)])

        edge_vals[i] = 1 - num/denom

    return edge_vals

# input: time between cells in minutes, dictionary of edge values, step length
def getMaturationWeights(time, edge_vals, step_length = 5):

    assert time % step_length == 0, "time between cells is not a multiple of time between images"

    time = time / step_length

    # get minimum and maximum
    min_val = min(edge_vals.keys())
    max_val = max(edge_vals.keys())

    if time <= min_val:
        return 0
    
    elif time >= max_val:
        return 1
    
    else:
        return edge_vals[time]

class MaturationDelayFunctionNormal(BaseDelayFunction):
    def __init__(self):
        cdf_vals, min_val, max_val = calcCdfValues()
        edge_vals = calcMaturationWeights(cdf_vals)
        edge_vals[max_val] = 1.0
        super().__init__("maturationDelayFunctionNormal", 2, ["m"], {"maturation_min":min_val, "maturation_max":max_val}, min_val, max_val, edge_vals)
