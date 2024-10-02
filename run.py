import os
import pickle
import shutil

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

def multiRun(rawFile, modelFolder, modelName, functionFolder, conjugationList, maturationList, colourList, force=False, naiveDepth=2, evaluateTimeout=10, debug=0, progressBar=False, rebuildFunctions=False):
    os.makedirs(functionFolder, exist_ok=True)
    os.makedirs(modelFolder, exist_ok=True)
    os.makedirs(modelFolder+modelName+"_data/", exist_ok=True)
    os.makedirs(modelFolder+modelName+"_results/", exist_ok=True)
    if rebuildFunctions:
        buildFunctions(functionFolder, force=True)
    for conjugation in conjugationList:
        for maturation in maturationList:
            for colour in colourList:
                fullRun(rawFile, modelFolder, modelFolder+modelName+"_data/", modelName, functionFolder, conjugation, maturation[0], maturation[1], maturation[2], colour[0], colour[1], colour[2], force=force, naiveDepth=naiveDepth, evaluateTimeout=evaluateTimeout, debug=debug, progressBar=progressBar, rebuildFunctions=False)
                shutil.copyfile(modelFolder+modelName+"_modeldata_"+conjugation+"_"+colour[0]+"_"+maturation[0]+"_maturationNormalized_queryEvaluation.pickle", modelFolder+modelName+"_results/"+modelName+"_modeldata_"+conjugation+"_"+colour[0]+"_"+maturation[0]+"_maturationNormalized_queryEvaluation.pickle")

# If force=False, then we don't recalculate any model or model data that already exists.
def fullRun(rawFile, modelFolder, dataFolder, modelName, functionFolder, conjugationFunctionName, maturationFunctionName, maturation_min, maturation_max, colourFunctionName, colour_min, colour_max, force=True, naiveDepth=2, evaluateTimeout=10, debug=0, progressBar=False, rebuildFunctions=False):

    if debug >= 1:
        print("Starting full run on " + modelName + " using functions " + conjugationFunctionName + ", " + colourFunctionName + ", and " + maturationFunctionName)

    modelExtension = "_None_None_"+str(maturation_min)+"-"+str(maturation_max)
    modelExtensionExtra = ""

    model = None # Saved model to reduce saving and loading.

    if rebuildFunctions:
        buildFunctions(functionFolder, force=force)

    if force or (not os.path.exists(dataFolder+modelName+"_raw.pickle")):
        saveAll(rawFile, dataFolder, modelName)

    if force or (not os.path.exists(modelFolder+modelName+"_model"+modelExtension+".pickle")):
        model = setupGraph(modelFolder, dataFolder, modelName, maturation_min, maturation_max, save=True, debug=debug, progressBar=progressBar, doCheck=False)
    else:
        model = None

    nextModelExtension = modelExtension.split("_")
    nextModelExtension[3] = maturationFunctionName
    nextModelExtension = "_".join(nextModelExtension)
    if force or (not os.path.exists(modelFolder+modelName+"_model"+nextModelExtension+".pickle")):
        model = addDelayFunctionToModel(modelFolder, dataFolder, modelName, modelExtension, functionFolder+maturationFunctionName+".pickle", save=True, debug=debug, progressBar=progressBar, loadedModel=model)
    else:
        model = None
    modelExtension = nextModelExtension
    
    nextModelExtension = modelExtension.split("_")
    nextModelExtension[1] = conjugationFunctionName
    nextModelExtension = "_".join(nextModelExtension)
    if force or (not os.path.exists(modelFolder+modelName+"_model"+nextModelExtension+".pickle")):
        model = addConjugationFunctionToModel(modelFolder, dataFolder, modelName, modelExtension, functionFolder+conjugationFunctionName+".pickle", save=True, debug=debug, progressBar=progressBar, loadedModel=model, saveConjugationFunction=True)
        if debug >= 1:
            print("Checking CPDs.")
        check = model.checkCPDs(progressBar=progressBar)
        assert check
        if debug >= 1:
            print("Passed? " + str(check))
    else:
        model = None
    modelExtension = nextModelExtension

    nextModelExtension = modelExtension.split("_")
    nextModelExtension[2] = str(colour_min) + "-" + str(colour_max)
    nextModelExtension = "_".join(nextModelExtension)
    if force or (not os.path.exists(modelFolder+modelName+"_modeldata"+nextModelExtension+"_evidence.pickle")) or (not os.path.exists(modelFolder+modelName+"_model"+nextModelExtension+"_contradictionsPruned.pickle")):
        model, evidence = get_evidence(modelFolder, dataFolder, modelName, modelExtension, colour_min, colour_max, save=True, debug=debug, progressBar=progressBar, loadedModel=model)
        if debug >= 1:
            print("Checking CPDs.")
        check = model.checkCPDs(progressBar=progressBar)
        assert check
        if debug >= 1:
            print("Passed? " + str(check))
    else:
        model = None
    modelExtension = nextModelExtension

    if force or (not os.path.exists(modelFolder+modelName+"_model"+modelExtension+"_conjugationNormalized.pickle")):
        model = normalizeConjugation(modelFolder, dataFolder, modelName, modelExtension, "_contradictionsPruned", normalizeTo=1.0, save=True, debug=debug, progressBar=progressBar, loadedModel=model, edgeList=None)
        if debug >= 1:
            print("Checking CPDs.")
        check = model.checkCPDWeights(progressBar=progressBar)
        assert check
        if debug >= 1:
            print("Passed? " + str(check))
    else:
        model = None

    if force or (not os.path.exists(modelFolder+modelName+"_modeldata"+modelExtension+"_queries.pickle")) or (not os.path.exists(modelFolder+modelName+"_modeldata"+modelExtension+"_criticalSegments.pickle")):
        queries, criticalSegments = find_queries(modelFolder, dataFolder, modelName, modelExtension, "_conjugationNormalized", save=True, debug=debug, progressBar=progressBar, loadedModel=model)

    preColourModelExtension = modelExtension
    nextModelExtension = modelExtension.split("_")
    nextModelExtension[2] = colourFunctionName
    nextModelExtension = "_".join(nextModelExtension)
    if force or (not os.path.exists(modelFolder+modelName+"_modeldata"+nextModelExtension+"_conjugationNormalized_naiveProbabilities.pickle")):
        naiveProbabilities = computeNaiveProbabilities(modelFolder, dataFolder, modelName, modelExtension, "_conjugationNormalized", functionFolder + colourFunctionName + ".pickle", depth=naiveDepth, save=True, debug=debug, progressBar=progressBar, loadedModel=None)
    modelExtension = nextModelExtension
    
    if force or (not os.path.exists(modelFolder+modelName+"_model"+modelExtension+"_maturationNormalized.pickle")):
        model = normalizeMaturation(modelFolder, dataFolder, modelName, preColourModelExtension, "_conjugationNormalized", modelExtension, normalizeTo=1.0, save=True, debug=debug, progressBar=progressBar, loadedModel=model, edgeList=None)
        if debug >= 1:
            print("Checking CPDs.")
        check = model.checkCPDWeights(progressBar=progressBar)
        assert check
        if debug >= 1:
            print("Passed? " + str(check))
    else:
        model = None

    if force or (not os.path.exists(modelFolder+modelName+"_modeldata"+modelExtension+"_evidence.pickle")) or (not os.path.exists(modelFolder+modelName+"_modeldata"+modelExtension+"_criticalSegments.pickle")) or (not os.path.exists(modelFolder+modelName+"_modeldata"+modelExtension+"_queries.pickle")):
        copyModelData(modelFolder, modelName, preColourModelExtension, modelExtension)

    if force or (not os.path.exists(modelFolder+modelName+"_modeldata"+modelExtension+"_maturationNormalized_naiveProbabilities.pickle")):
        naiveProbabilities = computeNaiveProbabilities(modelFolder, dataFolder, modelName, modelExtension, "_maturationNormalized", functionFolder+colourFunctionName+".pickle", depth=naiveDepth, save=True, debug=debug, progressBar=progressBar, loadedModel=model)

    if force or (not os.path.exists(modelFolder+modelName+"_modeldata"+modelExtension+"_maturationNormalized_querypoints.pickle")):
        queryPointValues = evaluateModel(modelFolder, dataFolder, modelName, modelExtension, "_maturationNormalized", timeout=evaluateTimeout, save=True, debug=debug, progressBar=progressBar, loadedModel=model)

    if force or (not os.path.exists(modelFolder+modelName+"_modeldata"+modelExtension+"_maturationNormalized_queryEvaluation.pickle")):
        queryEvaluation = mergeQueries(modelFolder, modelName, modelExtension, "_maturationNormalized", functionFolder+colourFunctionName+".pickle", save=True, debug=debug, progressBar=progressBar)

def resultsToCSV(resultsFolder, outputFile):
    results = []
    resultNames = []
    for filename in os.listdir(resultsFolder):
        if filename[-7:] != ".pickle":
            continue
        print("Loading " + str(filename))
        names = filename.split("_")
        resultNames.append("_".join(names[2:5]))
        with open(os.path.join(resultsFolder, filename), "rb") as f:
            results.append(pickle.load(f))
    keys = results[0].keys()
    for result in results:
        if result.keys() != keys:
            print("Key mismatch, stopping result merge.")
            return
    safeKeys = set()
    zeroKeys = set()
    for key in keys:
        safe = True
        allZero = True
        for result in results:
            if isinstance(result[key], tuple):
                safe = False
                break
            elif result[key] != 0:
                allZero = False
        if safe:
            if allZero:
                zeroKeys.add(key)
            else:
                safeKeys.add(key)
    safeKeys = sorted(safeKeys, key=lambda x: (int(x.split("_")[1]), int(x[1:].split("_")[0])))
    with open(outputFile, "w") as f:
        f.write(",".join(["Query"]+resultNames)+"\n")
        for key in safeKeys:
            f.write(",".join([key] + [str(results[i][key]) for i in range(len(results))]) + "\n")

    print("Usable queries: " + str(len(safeKeys)))
    print("Zero queries: " + str(len(zeroKeys)))
    print("All queries: " + str(len(keys)))

if __name__ == "__main__":
    if False:
        #def fullRun(rawFile, modelFolder, dataFolder, modelName, functionFolder, conjugationFunctionName, maturationFunctionName, maturation_min, maturation_max, colourFunctionName, colour_min, colour_max, timeout=10, debug=0, progressBar=False, rebuildFunctions=False)
        rawFile = "data/trap17v2test4/trap17raw.csv"
        modelFolder = "data/trap17v2test4/"
        dataFolder = "data/trap17v2test4/trap17v2test4_data/"
        modelName = "trap17v2test4"
        functionFolder = "functions/"
        conjugationFunctionName = "contactWeightsBaseline"
        maturationFunctionName = "maturationDelayFunctionUniform"
        maturation_min = 6
        maturation_max = 18
        colourFunctionName = "colourDelayFunctionUniformV2"
        colour_min = 6
        colour_max = 30
        timeout = 20
        debug = 1
        progressBar = True
        rebuildFunctions = True
        fullRun(rawFile, modelFolder, dataFolder, modelName, functionFolder, conjugationFunctionName, maturationFunctionName, maturation_min, maturation_max, colourFunctionName, colour_min, colour_max, timeout=timeout, debug=debug, progressBar=progressBar, rebuildFunctions=rebuildFunctions)
    if True:
        #multiRun(rawFile, metaFolder, baseModelName, functionFolder, conjugationList, maturationList, colourList, timeout=10, debug=0, progressBar=False, rebuildFunctions=False)
        #multiRun("data/trap16bulk/trap16raw.csv", "data/trap16bulk/", "trap16bulk", "functions/", ["contactWeightsBoundary", "contactWeightsBaseline"], [("maturationUniformV2L15U75", 3, 15), ("maturationUniformV2L30U90", 6, 18)], [("colourUniformV2L30U120", 6, 24), ("colourUniformV2L30U150", 6, 30)], force=False, naiveDepth=2, evaluateTimeout=10, debug=1, progressBar=True, rebuildFunctions=True)
        multiRun("data/trap9bulk/trap9raw.csv", "data/trap9bulk/", "trap9bulk", "functions/", ["contactWeightsBoundary", "contactWeightsBaseline"], [("maturationUniformV2L15U75", 3, 15), ("maturationUniformV2L30U90", 6, 18)], [("colourUniformV2L30U120", 6, 24), ("colourUniformV2L30U150", 6, 30)], force=False, naiveDepth=2, evaluateTimeout=10, debug=1, progressBar=True, rebuildFunctions=True)
