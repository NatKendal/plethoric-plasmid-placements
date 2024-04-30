import os
import pickle

from plplpl.base_functions import *
from plplpl.conjugation_functions import *
from plplpl.colour_delay_functions import *
from plplpl.maturation_delay_functions import *

def buildFunctions(location, force=False):
    functionList = [
        ContactWeightsFixedNaive(),
        MaturationDelayFunctionNormal(),
        ColourDelayFunctionUniform()
    ]

    for function in functionList:
        if force or (not os.path.isfile(location+function.name+".pickle")):
            function.save(location+function.name+".pickle")
