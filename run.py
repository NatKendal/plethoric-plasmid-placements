import pickle

from plplpl.data_processing.first_processing import saveAll
from plplpl.base_functions import buildFunctions

from plplpl.model import setupGraph
from plplpl.model import addDelayFunctionToModel
from plplpl.model import addConjugationFunctionToModel

from plplpl.model import normalizeConjugation

from plplpl.model import get_evidence

from plplpl.model import find_queries

from plplpl.model import computeNaiveProbabilities

from plplpl.model import copyModelData

from plplpl.model import normalizeMaturation

from plplpl.model import evaluateModel

from plplpl.model import mergeQueries

def fullRun(rawFile, modelFolder, dataFolder, modelName, functionFolder, conjugationFunctionName, maturationFunctionName, maturation_min, maturation_max, colourFunctionName, colour_min, colour_max, timeout=10, debug=0, progressBar=False, rebuildFunctions=False):

    modelExtension = "_None_None_None"
    modelExtensionExtra = ""
    saveAll(rawFile, dataFolder, modelName)

    if rebuildFunctions:
        buildFunctions(functionFolder, force=True)

    model = setupGraph(modelFolder, dataFolder, modelName, maturation_min, maturation_max, save=True, debug=debug, progressBar=progressBar, doCheck=False)

    model = addDelayFunctionToModel(modelFolder, dataFolder, modelName, modelExtension, functionFolder+maturationFunctionName+".pickle", save=True, debug=1, progressBar=True, loadedModel=None)
    modelExtension = modelExtension.split("_")
    modelExtension[3] = maturationFunctionName
    modelExtension = "_".join(modelExtension)
    
    model = addConjugationFunctionToModel(modelFolder, dataFolder, modelName, modelExtension, functionFolder+conjugationFunctionName+".pickle", save=True, debug=1, progressBar=True, loadedModel=None, saveConjugationFunction=True)
    print("Checking CPDs.")
    print("Passed? " + str(model.checkCPDs(progressBar=True)))
    modelExtension = modelExtension.split("_")
    modelExtension[1] = conjugationFunctionName
    modelExtension = "_".join(modelExtension)

    model, evidence = get_evidence(modelFolder, dataFolder, modelName, modelExtension, colour_min, colour_max, save=True, debug=1, progressBar=True, loadedModel=None)
    print("Checking CPDs.")
    print("Passed? " + str(model.checkCPDs(progressBar=True)))
    modelExtension = modelExtension.split("_")
    modelExtension[2] = str(colour_min) + "-" + str(colour_max)
    modelExtension = "_".join(modelExtension)

    model = normalizeConjugation(modelFolder, dataFolder, modelName, modelExtension, "_contradictionsPruned", normalizeTo=1.0, save=True, debug=1, progressBar=True, loadedModel=None, edgeList=None)
    print("Checking CPDs.")
    print("Passed? " + str(model.checkCPDs(progressBar=True)))

    queries, criticalSegments = find_queries(modelFolder, dataFolder, modelName, modelExtension, "_conjugationNormalized", save=True, debug=1, progressBar=True, loadedModel=None)

    naiveProbabilities = computeNaiveProbabilities(modelFolder, dataFolder, modelName, modelExtension, "_conjugationNormalized", functionFolder + colourFunctionName + ".pickle", depth=2, save=True, debug=1, progressBar=True, loadedModel=None)
    preColourModelExtension = modelExtension
    modelExtension = modelExtension.split("_")
    modelExtension[2] = colourFunctionName
    modelExtension = "_".join(modelExtension)

    model = normalizeMaturation(modelFolder, dataFolder, modelName, preColourModelExtension, "_conjugationNormalized", modelExtension, normalizeTo=1.0, save=True, debug=1, progressBar=True, loadedModel=None, edgeList=None)
    print("Checking CPDs.")
    print("Passed? " + str(model.checkCPDs(progressBar=True)))

    copyModelData(modelFolder, modelName, preColourModelExtension, modelExtension)

    naiveProbabilities = computeNaiveProbabilities(modelFolder, dataFolder, modelName, modelExtension, "_maturationNormalized", functionFolder+colourFunctionName+".pickle", depth=2, save=True, debug=1, progressBar=True, loadedModel=None)

    queryPointValues = evaluateModel(modelFolder, dataFolder, modelName, modelExtension, "_maturationNormalized", timeout=10, save=True, debug=1, progressBar=True, loadedModel=None)

    queryEvaluation = mergeQueries(modelFolder, modelName, modelExtension, "_maturationNormalized", ".pickle", save=True, debug=1, progressBar=True)

if __name__ == "__main__":
    #def fullRun(rawFile, modelFolder, dataFolder, modelName, functionFolder, conjugationFunctionName, maturationFunctionName, maturation_min, maturation_max, colourFunctionName, colour_min, colour_max, timeout=10, debug=0, progressBar=False, rebuildFunctions=False)
    rawFile = "data/trap17v2test2/trap17raw.csv"
    modelFolder = "data/trap17v2test2/"
    dataFolder = "data/trap17v2test2/trap17v2test2_data/"
    modelName = "trap17v2test2"
    functionFolder = "functions/"
    conjugationFunctionName = "contactWeightsBaseline"
    maturationFunctionName = "maturationDelayFunctionUniform"
    maturation_min = 6
    maturation_max = 18
    colourFunctionName = "colourDelayFunctionUniformV2"
    colour_min = 6
    colour_max = 30
    timeout = 15
    debug = 1
    progressBar = True
    rebuildFunctions = True
    fullRun(rawFile, modelFolder, dataFolder, modelName, functionFolder, conjugationFunctionName, maturationFunctionName, maturation_min, maturation_max, colourFunctionName, colour_min, colour_max, timeout=timeout, debug=debug, progressBar=progressBar, rebuildFunctions=rebuildFunctions)
