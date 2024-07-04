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

from plplpl.NoisyOr import NoisyOrFactor

from multiprocessing import Pool
from pathlib import Path
from weakref import WeakValueDictionary

import gc
import pickle
import sys



def experiment1():
    metaFolder = "Experiment1/"
    functionFolder = "functions/"
    inputFiles = ["Experiment1/rawFiles/trap6.csv"]
    transmissionFunctions = ["contactWeightsBaseline", "contactWeightsBoundary", "contactWeightsArea"]
    colourFunctions = [("colourDelayFunctionNormal", 6, 30), ("colourDelayFunctionNormalSkewed", 6, 30)]
    maturationFunctions = [("maturationDelayFunctionUniform", 6, 18), ("maturationDelayFunctionNormal", 6, 18)]
    multiCoreMultiRun(metaFolder, functionFolder, inputFiles, transmissionFunctions, colourFunctions, maturationFunctions, save=True, cores=6)
    #multiRun(metaFolder, functionFolder, inputFiles, transmissionFunctions, colourFunctions, maturationFunctions, save=True, progressBar=True, debug=1)

def experiment2():
    metaFolder = "Experiment1/"
    functionFolder = "functions/"
    inputFiles = ["Experiment1/rawFiles/trap9.csv"]
    transmissionFunctions = ["contactWeightsBaseline", "contactWeightsBoundary", "contactWeightsArea"]
    colourFunctions = [("colourDelayFunctionNormal", 6, 30), ("colourDelayFunctionNormalSkewed", 6, 30)]
    maturationFunctions = [("maturationDelayFunctionUniform", 6, 18), ("maturationDelayFunctionNormal", 6, 18)]
    #multiCoreMultiRun(metaFolder, functionFolder, inputFiles, transmissionFunctions, colourFunctions, maturationFunctions, save=True, cores=6)
    multiRun(metaFolder, functionFolder, inputFiles, transmissionFunctions, colourFunctions, maturationFunctions, save=True, progressBar=True, debug=1)

def experiment3(trapFile):
    metaFolder = "Experiment3/"
    functionFolder = "functions/"
    inputFiles = ["Experiment3/rawFiles/" + trapFile]
    transmissionFunctions = ["contactWeightsBaseline", "contactWeightsBoundary", "contactWeightsArea"]
    colourFunctions = [("colourDelayFunctionNormal", 6, 30), ("colourDelayFunctionNormalSkewed", 6, 30)]
    maturationFunctions = [("maturationDelayFunctionNormal", 6, 18), ("maturationDelayFunctionNormalSkewed", 6, 18)]
    multiRun(metaFolder, functionFolder, inputFiles, transmissionFunctions, colourFunctions, maturationFunctions, save=True, progressBar=True, debug=1)

def experiment4():
    metaFolder = "Experiment4/"
    functionFolder = "functions/"
    inputFiles = ["Experiment4/rawFiles/trap1.csv", "Experiment4/rawFiles/trap2.csv", "Experiment4/rawFiles/trap3.csv", "Experiment4/rawFiles/trap4.csv", "Experiment4/rawFiles/trap5.csv", "Experiment4/rawFiles/trap8.csv", "Experiment4/rawFiles/trap10.csv", "Experiment4/rawFiles/trap11.csv", "Experiment4/rawFiles/trap13.csv", "Experiment4/rawFiles/trap14.csv", "Experiment4/rawFiles/trap18.csv"]
    transmissionFunctions = ["contactWeightsBaseline", "contactWeightsBoundary", "contactWeightsArea"]
    colourFunctions = [("colourDelayFunctionNormal", 6, 30), ("colourDelayFunctionNormalSkewed", 6, 30)]
    maturationFunctions = [("maturationDelayFunctionNormal", 6, 18), ("maturationDelayFunctionNormalSkewed", 6, 18)]
    multiRun(metaFolder, functionFolder, inputFiles, transmissionFunctions, colourFunctions, maturationFunctions, save=True, progressBar=True, debug=1)

# NOTE: This might not work properly with the NoisyOrFactor registry!
def multiCoreMultiRun(metaFolder, functionFolder, inputFiles, transmissionFunctions, colourFunctions, maturationFunctions, save=True, cores=1):
    metaFolder = metaFolder + "/"
    Path(metaFolder).mkdir(parents=True, exist_ok=True)
    run = 1
    with Pool(processes=cores) as pool:
        for inputFile in inputFiles:
            for transmissionFunction in transmissionFunctions:
                for colourFunction in colourFunctions:
                    for maturationFunction in maturationFunctions:
                        name = inputFile.split("/")[-1].split(".")[0] + "_" + transmissionFunction + "_" + colourFunction[0] + "_" + maturationFunction[0] + "_runID" + str(run)
                        Path(metaFolder+name+"/").mkdir(parents=True, exist_ok=True)
                        Path(metaFolder+name+"/"+name+"_data/").mkdir(parents=True, exist_ok=True)
                        print("Setup - " + name)
                        pool.apply_async(fullRun, args=(inputFile, name, metaFolder+name+"/", metaFolder+name+"/"+name+"_data/", functionFolder, transmissionFunction, colourFunction[0], colourFunction[1], colourFunction[2], maturationFunction[0], maturationFunction[1], maturationFunction[2]))
                        run += 1
        pool.close()
        pool.join()


# Expects inputFiles and tranmissionFunctions to be lists.
# Expects colourFunctions and maturationFunctions to be lists of tuples (functionName, min, max).
def multiRun(metaFolder, functionFolder, inputFiles, transmissionFunctions, colourFunctions, maturationFunctions, save=True, progressBar=False, debug=0, rebuildFunctions=False):
    metaFolder = metaFolder + "/"
    Path(metaFolder).mkdir(parents=True, exist_ok=True)
    run = 1
    total = len(inputFiles) * len(transmissionFunctions) * len(colourFunctions) * len(maturationFunctions)
    for inputFile in inputFiles:
        for transmissionFunction in transmissionFunctions:
            for colourFunction in colourFunctions:
                for maturationFunction in maturationFunctions:
                    name = inputFile.split("/")[-1].split(".")[0] + "_" + transmissionFunction + "_" + colourFunction[0] + "_" + maturationFunction[0] + "_runID" + str(run)
                    run += 1
                    Path(metaFolder+name+"/").mkdir(parents=True, exist_ok=True)
                    Path(metaFolder+name+"/"+name+"_data/").mkdir(parents=True, exist_ok=True)
                    evaluation = fullRun(inputFile, name, metaFolder+name+"/", metaFolder+name+"/"+name+"_data/", functionFolder, transmissionFunction, colourFunction[0], colourFunction[1], colourFunction[2], maturationFunction[0], maturationFunction[1], maturationFunction[2], save=save, progressBar=progressBar, debug=debug, rebuildFunctions=rebuildFunctions)
                    with open(metaFolder+name+"_evaluation.pickle", "wb") as f:
                        pickle.dump(evaluation, f)
                    print("Finished run " + str(run) + " of " + str(total))
                    rebuildFunctions = False

def fullRun(rawInputFile, modelName, modelFolder, dataFolder, functionFolder, transmissionFunctionName, colourDelayFunctionName, colour_min, colour_max, maturationDelayFunctionName, maturation_min, maturation_max, save=True, progressBar=False, debug=0, rebuildFunctions=False, statusUpdates=True, clearRegistry=True):

    if clearRegistry:
        NoisyOrFactor._registry = WeakValueDictionary()

    modelFolder = modelFolder + "/"
    dataFolder = dataFolder + "/"
    functionFolder = functionFolder + "/"

    if statusUpdates:
        print("Stage 0 - " + modelName + " - Initialized")

    if rebuildFunctions:
        buildFunctions(functionFolder, force=True)
    saveAll(rawInputFile, dataFolder, modelName)
    saveSynchronyCertainty(dataFolder, modelName)

    if statusUpdates:
        print("Stage 1 - " + modelName + " - Finished Calculating Data")

    gc.collect()

    model = setupGraph(modelFolder, dataFolder, modelName, colour_min, colour_max, maturation_min, maturation_max, save=save, debug=debug, progressBar=progressBar)

    gc.collect()

    if statusUpdates:
        print("Stage 2 - " + modelName + " - Finished Generating Directed Graph")

    model = addDelayFunctionToModel(modelFolder, dataFolder, modelName, "_None_None_None", functionFolder+colourDelayFunctionName+".pickle", save=save, debug=debug, progressBar=progressBar, safeMode=False, loadedModel=model)

    gc.collect()

    if statusUpdates:
        print("Stage 3 - " + modelName + " - Finished Adding Colour Function")

    model = addDelayFunctionToModel(modelFolder, dataFolder, modelName, "_None_"+colourDelayFunctionName+"_None", functionFolder+maturationDelayFunctionName+".pickle", save=save, debug=debug, progressBar=progressBar, safeMode=False, loadedModel=model)

    gc.collect()

    if statusUpdates:
        print("Stage 4 - " + modelName + " - Finished Adding Maturation Function")

    model = addConjugationFunctionToModel(modelFolder, dataFolder, modelName, "_None_"+colourDelayFunctionName+"_"+maturationDelayFunctionName, functionFolder+transmissionFunctionName+".pickle", save=save, debug=debug, progressBar=progressBar, safeMode=False, loadedModel=model, saveConjugationFunction=True)

    gc.collect()

    if statusUpdates:
        print("Stage 5 - " + modelName + " - Finished Adding Transmission Function")

    model, evidence, forwardLinks, backwardLinks = get_evidence(modelFolder, dataFolder, modelName, "_"+transmissionFunctionName+"_"+colourDelayFunctionName+"_"+maturationDelayFunctionName, save=save, debug=debug, progressBar=progressBar, loadedModel=model)

    gc.collect()

    if statusUpdates:
        print("Stage 6 - " + modelName + " - Finished Calculating Evidence")

    conQueries, nonConQueries, allEdges, fullQueries = build_queries(modelFolder, dataFolder, modelName, "_"+transmissionFunctionName+"_"+colourDelayFunctionName+"_"+maturationDelayFunctionName, save=save, debug=debug, progressBar=progressBar, loadedModel=model, loadedEvidence=evidence, loadedForwardLinks=forwardLinks, loadedBackwardLinks=backwardLinks)

    gc.collect()

    if statusUpdates:
        print("Stage 7 - " + modelName + " - Finished Calculating Queries")

    model = normalizeModel(modelFolder, dataFolder, modelName, "_"+transmissionFunctionName+"_"+colourDelayFunctionName+"_"+maturationDelayFunctionName, dataFolder+modelName+"_"+transmissionFunctionName+".pickle", normalizeTo=1.0, save=save, debug=debug, progressBar=progressBar, loadedModel=model, loadedEdgeList=allEdges, loadedBackwardLinks=backwardLinks)

    gc.collect()

    if statusUpdates:
        print("Stage 8 - " + modelName + " - Finished Normalizing Model")

    naiveProbabilities, incomingProbabilities = calculateNaiveProbabilities(modelFolder, dataFolder, modelName, "_"+transmissionFunctionName+"_"+colourDelayFunctionName+"_"+maturationDelayFunctionName, save=save, debug=debug, progressBar=progressBar, loadedModel=model, loadedEvidence=evidence, loadedConjugateQueries=conQueries, loadedFullQueries=fullQueries)

    gc.collect()

    if statusUpdates:
        print("Stage 9 - " + modelName + " - Finished Calculating Naive Probabilities")

    conQueryEvaluation, ncQueryEvaluation, preconditonConstants, fullQueryEvaluation = evaluateModel(modelFolder, dataFolder, modelName, "_"+transmissionFunctionName+"_"+colourDelayFunctionName+"_"+maturationDelayFunctionName, save=save, debug=debug, progressBar=progressBar, loadedModel=model, loadedEvidence=evidence, loadedConjugateQueries=conQueries, loadedNonConjugateQueries=nonConQueries, loadedNaiveProbabilities=naiveProbabilities, loadedIncomingProbabilities=incomingProbabilities)

    print("Finished - " + modelName)

    return fullQueryEvaluation

def reloadFunctions(functionFolder="functions/"):
    buildFunctions(functionFolder, force=True)

# Borrowed from stackoverflow/python before it was removed.
def strtobool(val):
    """Convert a string representation of truth to true (1) or false (0).
    True values are 'y', 'yes', 't', 'true', 'on', and '1'; false values
    are 'n', 'no', 'f', 'false', 'off', and '0'.  Raises ValueError if
    'val' is anything else.
    """
    val = val.lower()
    if val in ('y', 'yes', 't', 'true', 'on', '1'):
        return 1
    elif val in ('n', 'no', 'f', 'false', 'off', '0'):
        return 0
    else:
        raise ValueError("invalid truth value %r" % (val,))

if __name__ == "__main__":
    if True:
        if len(sys.argv) == 18:
            fullRun(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6], sys.argv[7], int(sys.argv[8]), int(sys.argv[9]), sys.argv[10], int(sys.argv[11]), int(sys.argv[12]), save=bool(strtobool(sys.argv[13])), progressBar=bool(strtobool(sys.argv[14])), debug=int(sys.argv[15]), rebuildFunctions=bool(strtobool(sys.argv[16])), statusUpdates=bool(strtobool(sys.argv[17])), clearRegistry=bool(strtobool(sys.argv[18])))
        else:
            print("Expected 17 arguments for:")
            print("fullRun(rawInputFile, modelName, modelFolder, dataFolder, functionFolder, transmissionFunctionName, colourDelayFunctionName, colour_min, colour_max, maturationDelayFunctionName, maturation_min, maturation_max, save=True, progressBar=False, debug=0, rebuildFunctions=False, statusUpdates=True, clearRegistry=True)")
