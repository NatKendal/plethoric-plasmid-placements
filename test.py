from plplpl.base_functions import buildFunctions

from plplpl.data_processing.first_processing import saveAll
from plplpl.data_processing.synchronyCheck import saveSynchronyCertainty

from plplpl.model import setupGraph
from plplpl.model import addDelayFunctionToModel
from plplpl.model import addConjugationFunctionToModel

from plplpl.model import get_evidence

from plplpl.model import build_queries

from plplpl.model import normalizeModel

from plplpl.model import calculateNaiveProbabilities

from plplpl.model import evaluateModel

import gc
import pickle
import sys

# When writing general test, use gc.collect to force garbage collector to clean up after each step.

if __name__ == "__main__":
    sys.settrace

    modelFolder = "data/trap6test/"
    dataFolder = "data/trap6test/trap6test_data/"
    modelName = "trap6test"

    #saveAll("data/trap6test/trap6test_raw.csv", dataFolder, modelName)
    #buildFunctions("functions/", force=True)
    #saveSynchronyCertainty(dataFolder, modelName)
    #model = setupGraph(modelFolder, dataFolder, modelName, 6, 30, 6, 18, save=True, debug=1, progressBar=True)
    gc.collect()
    #model = addDelayFunctionToModel(modelFolder, dataFolder, modelName, "_None_None_None", "functions/colourDelayFunctionUniform.pickle", save=True, debug=1, progressBar=True, safeMode=False, loadedModel=None)
    gc.collect()
    #model = addDelayFunctionToModel(modelFolder, dataFolder, modelName, "_None_colourDelayFunctionUniform_None", "functions/maturationDelayFunctionNormal.pickle", save=True, debug=1, progressBar=True, safeMode=False, loadedModel=None)
    gc.collect()
    #model = addConjugationFunctionToModel(modelFolder, dataFolder, modelName, "_None_colourDelayFunctionUniform_maturationDelayFunctionNormal", "functions/contactWeightsBaseline.pickle", save=True, debug=1, progressBar=True, safeMode=False, loadedModel=None, saveConjugationFunction=True)
    gc.collect()
    #model, evidence, forwardLinks, backwardLinks = get_evidence(modelFolder, dataFolder, modelName, "_contactWeightsBaseline_colourDelayFunctionUniform_maturationDelayFunctionNormal", save=True, debug=1, progressBar=True, loadedModel=None)
    gc.collect()
    #conQueries, nonConQueries, allEdges, fullQueries = build_queries(modelFolder, dataFolder, modelName, "_contactWeightsBaseline_colourDelayFunctionUniform_maturationDelayFunctionNormal", save=True, debug=1, progressBar=True, loadedModel=None, loadedEvidence=None, loadedForwardLinks=None, loadedBackwardLinks=None)
    gc.collect()
    #model = normalizeModel(modelFolder, dataFolder, modelName, "_contactWeightsBaseline_colourDelayFunctionUniform_maturationDelayFunctionNormal", "data/trap6test/trap6test_data/trap6test_contactWeightsBaseline.pickle", normalizeTo=200.0, save=True, debug=1, progressBar=True, loadedModel=None, loadedEdgeList=None, loadedBackwardLinks=None)
    gc.collect()
    #naiveProbabilities, incomingProbabilities = calculateNaiveProbabilities(modelFolder, dataFolder, modelName, "_contactWeightsBaseline_colourDelayFunctionUniform_maturationDelayFunctionNormal", save=True, debug=1, progressBar=True, loadedModel=None, loadedEvidence=None, loadedConjugateQueries=None, loadedFullQueries=None)
    gc.collect()
    conQueryEvaluation, ncQueryEvaluation, preconditonConstants, fullQueryEvaluation = evaluateModel(modelFolder, dataFolder, modelName, "_contactWeightsBaseline_colourDelayFunctionUniform_maturationDelayFunctionNormal", save=True, debug=1, progressBar=True, loadedModel=None, loadedEvidence=None, loadedConjugateQueries=None, loadedNonConjugateQueries=None, loadedNaiveProbabilities=None)
    gc.collect()

    if False:
        with open("data/trap6test/trap6test_data/trap6test_humanFriendlyNameLookup.pickle", "rb") as f:
            uid = pickle.load(f)

        with open("data/trap6test/trap6test_data/trap6test_humanFriendlyName.pickle", "rb") as f:
            name = pickle.load(f)
    
