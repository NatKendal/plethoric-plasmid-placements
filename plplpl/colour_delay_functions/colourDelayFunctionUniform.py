# specify the cdf, get the edge weights for colour

import math
import numpy as np
import pathlib

from plplpl.base_functions import BaseDelayFunction

# input: minimum and maximum of range in minutes, step length
# output: dictionary of uniform distribution values for all time steps in the range 
def calcCdfValues(min_val=30, max_val=150, step_length = 5):

    # convert minutes to steps
    # decrease min by one step so non-zero probability at the actual minimum
    # probability will be equal one at max
    min_val = math.floor(min_val/step_length) - 1
    max_val = math.ceil(max_val/step_length)

    cdf_vals = dict()

    for i in range(min_val + 1, max_val+1):

        cdf_vals[i] = (i - min_val) / (max_val - min_val)

    return cdf_vals, min_val + 1, max_val

# input: dictionary of cdf values
# output: dictionary of edge weights
def calcColourWeights(cdf_vals):

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
def getColourWeights(time, edge_vals, step_length = 5):

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

class ColourDelayFunctionUniform(BaseDelayFunction):
    def __init__(self):
        cdf_vals, min_val, max_val = calcCdfValues()
        edge_vals = calcColourWeights(cdf_vals)
        super().__init__("colourDelayFunctionUniform", 1, ["c"], {"colour_min":(min_val), "colour_max":(max_val)}, min_val, max_val, edge_vals)

