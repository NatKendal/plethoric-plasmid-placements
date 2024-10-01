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
        ColourDelayFunctionUniformV2(),
        ContactWeightsBaseline(),
        MaturationDelayFunctionNormalShort(),
        MaturationDelayFunctionUniform(),
        ContactWeightsArea(),
        ContactWeightsBoundary(),
        MaturationDelayFunctionNormalSkewed(),
        MaturationUniformV2L15U75(),
        MaturationUniformV2L30U90(),
        ColourUniformV2L30U120(),
        ColourUniformV2L30U150()
    ]

    for function in functionList:
        if force or (not os.path.isfile(location+function.name+".pickle")):
            function.save(location+function.name+".pickle")
