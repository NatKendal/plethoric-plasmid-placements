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

from multiprocessing import Pool


import gc
from pathlib import Path
import pickle
import sys

def experiment1():
    metaFolder = "/Users/josephmeleshko/Home/GitHub/plethoric-plasmid-placements/Experiment1/"
    functionFolder = "/Users/josephmeleshko/Home/GitHub/plethoric-plasmid-placements/functions"
    inputFiles = ["/Users/josephmeleshko/Home/GitHub/plethoric-plasmid-placements/Experiment1/rawFiles/trap6.csv"]
    transmissionFunctions = ["ContactWeightsBaseline", "ContactWeightsBoundary", "ContactWeightsArea"]
    colourFunctions = [("ColourDelayFunctionNormal", 6, 30), ("ColourDelayFunctionNormalSkewed", 6, 30)]
    maturationFunctions = [("MaturationDelayFunctionUniform", 6, 18), ("MaturationDelayFunctionNormal", 6, 18)]
    multiRun(metaFolder, "/Users/josephmeleshko/Home/GitHub/plethoric-plasmid-placements/functions/", inputFiles, transmissionFunctions, colourFunctions, maturationFunctions, save=True, progressBar=True, debug=1, rebuildFunctions=True)

def multiCoreMultiRun(metaFolder, functionFolder, inputFiles, transmissionFunctions, colourFunctions, maturationFunctions, save=True, cores=1):
    metaFolder = metaFolder + "/"
    Path(metaFolder).mkdir(parents=True, exist_ok=True)
    run = 1
    with Pool(processes=cores) as pool:
        for inputFile in inputFiles:
            for transmissionFunction in transmissionFunctions:
                for colourFunction in colourFunctions:
                    for maturationFunction in maturationFunctions:
                        name = inputFile.split(".")[0] + "_" + transmissionFunction + "_" + colourFunction[0] + "_" + maturationFunction[0] + "_runID" + str(run)
                        Path(metaFolder+name+"/").mkdir(parents=True, exist_ok=True)
                        Path(metaFolder+name+"/"+name+"_data/").mkdir(parents=True, exist_ok=True)
                        pool.apply_async(fullRun, (inputFile, name, metaFolder+name+"/", metaFolder+name+"/"+name+"_data/", functionFolder, transmissionFunction, colourFunction[0], colourFunction[1], colourFunction[2], maturationFunction[0], maturationFunction[1], maturationFunction[2]))
        pool.close()
        pool.join()


# Expects inputFiles and tranmissionFunctions to be lists.
# Expects colourFunctions and maturationFunctions to be lists of tuples (functionName, min, max).
def multiRun(metaFolder, functionFolder, inputFiles, transmissionFunctions, colourFunctions, maturationFunctions, save=True, progressBar=False, debug=0, rebuildFunctions=False):
    metaFolder = metaFolder + "/"
    Path(metaFolder).mkdir(parents=True, exist_ok=True)
    run = 1
    total = len(inputFiles) * len(transmissionFunction) * len(colourFunction) * len(maturationFunction)
    for inputFile in inputFiles:
        for transmissionFunction in transmissionFunctions:
            for colourFunction in colourFunctions:
                for maturationFunction in maturationFunctions:
                    name = inputFile.split(".")[0] + "_" + transmissionFunction + "_" + colourFunction[0] + "_" + maturationFunction[0] + "_runID" + str(run)
                    run += 1
                    Path(metaFolder+name+"/").mkdir(parents=True, exist_ok=True)
                    Path(metaFolder+name+"/"+name+"_data/").mkdir(parents=True, exist_ok=True)
                    evaluation = fullRun(inputFile, name, metaFolder+name+"/", metaFolder+name+"/"+name+"_data/", functionFolder, transmissionFunction, colourFunction[0], colourFunction[1], colourFunction[2], maturationFunction[0], maturationFunction[1], maturationFunction[2], save=save, progressBar=progressBar, debug=debug, rebuildFunctions=rebuildFunctions)
                    with open(metaFolder+name+"_evaluation.pickle", "wb") as f:
                        pickle.dump(evaluation, f)
                    print("Finished run " + str(run) + " of " + str(total))
                    rebuildFunctions = False

def fullRun(rawInputFile, modelName, modelFolder, dataFolder, functionFolder, transmissionFunctionName, colourDelayFunctionName, colour_min, colour_max, maturationDelayFunctionName, maturation_min, maturation_max, save=True, progressBar=False, debug=0, rebuildFunctions=False):

    modelFolder = modelFolder + "/"
    dataFolder = dataFolder + "/"
    functionFolder = functionFolder + "/"

    if rebuildFunctions:
        buildFunctions(functionFolder, force=True)
    saveAll(rawInputFile, dataFolder, modelName)
    saveSynchronyCertainty(dataFolder, modelName)

    gc.collect()

    model = setupGraph(modelFolder, dataFolder, modelName, colour_min, colour_max, maturation_min, maturation_max, save=save, debug=debug, progressBar=progressBar)

    gc.collect()

    model = addDelayFunctionToModel(modelFolder, dataFolder, modelName, "_None_None_None", functionFolder+colourDelayFunctionName+".pickle", save=save, debug=debug, progressBar=progressBar, safeMode=False, loadedModel=model)

    gc.collect()

    model = addDelayFunctionToModel(modelFolder, dataFolder, modelName, "_None_"+colourDelayFunctionName+"_None", functionFolder+colourDelayFunctionName+".pickle", save=save, debug=debug, progressBar=progressBar, safeMode=False, loadedModel=model)

    gc.collect()

    model = addConjugationFunctionToModel(modelFolder, dataFolder, modelName, "_None_"+colourDelayFunctionName+"_"+maturationDelayFunctionName, functionFolder+transmissionFunctionName+".pickle", save=save, debug=debug, progressBar=progressBar, safeMode=False, loadedModel=model, saveConjugationFunction=True)

    gc.collect()

    model, evidence, forwardLinks, backwardLinks = get_evidence(modelFolder, dataFolder, modelName, "_"+transmissionFunctionName+"_"+colourDelayFunctionName+"_"+maturationDelayFunctionName, save=save, debug=debug, progressBar=progressBar, loadedModel=model)

    gc.collect()

    conQueries, nonConQueries, allEdges, fullQueries = build_queries(modelFolder, dataFolder, modelName, "_"+transmissionFunctionName+"_"+colourDelayFunctionName+"_"+maturationDelayFunctionName, save=save, debug=debug, progressBar=progressBar, loadedModel=model, loadedEvidence=evidence, loadedForwardLinks=forwardLinks, loadedBackwardLinks=backwardLinks)

    gc.collect()

    model = normalizeModel(modelFolder, dataFolder, modelName, "_"+transmissionFunctionName+"_"+colourDelayFunctionName+"_"+maturationDelayFunctionName, dataFolder+modelName+"_"+transmissionFunctionName+".pickle", normalizeTo=1.0, save=save, debug=debug, progressBar=progressBar, loadedModel=model, loadedEdgeList=allEdges, loadedBackwardLinks=backwardLinks)

    gc.collect()

    naiveProbabilities, incomingProbabilities = calculateNaiveProbabilities(modelFolder, dataFolder, modelName, "_"+transmissionFunctionName+"_"+colourDelayFunctionName+"_"+maturationDelayFunctionName, save=save, debug=debug, progressBar=progressBar, loadedModel=model, loadedEvidence=evidence, loadedConjugateQueries=conQueries, loadedFullQueries=fullQueries)

    gc.collect()

    conQueryEvaluation, ncQueryEvaluation, preconditonConstants, fullQueryEvaluation = evaluateModel(modelFolder, dataFolder, modelName, "_"+transmissionFunctionName+"_"+colourDelayFunctionName+"_"+maturationDelayFunctionName, save=save, debug=debug, progressBar=progressBar, loadedModel=model, loadedEvidence=evidence, loadedConjugateQueries=conQueries, loadedNonConjugateQueries=nonConQueries, loadedNaiveProbabilities=naiveProbabilities, loadedIncomingProbabilities=incomingProbabilities)

    return fullQueryEvaluation


def strtobool (val):
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
    if len(sys.argv) == 17:
        fullRun(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6], sys.argv[7], int(sys.argv[8]), int(sys.argv[9]), sys.argv[10], int(sys.argv[11]), int(sys.argv[12]), save=bool(strtobool(sys.argv[13])), progressBar=bool(strtobool(sys.argv[14])), debug=int(sys.argv[15]), rebuildFunctions=bool(strtobool(sys.argv[16])))
