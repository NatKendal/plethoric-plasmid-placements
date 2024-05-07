from plplpl.data_processing.first_processing import saveAll
from plplpl.data_processing.synchronyCheck import saveSynchronyCertainty

from plplpl.model import setupGraph
from plplpl.model import addDelayFunctionToModel
from plplpl.model import addConjugationFunctionToModel

from plplpl.model import get_evidence

from plplpl.model import build_queries

if __name__ == "__main__":
    modelFolder = "data/trap6test/"
    dataFolder = "data/trap6test/trap6test_data/"
    modelName = "trap6test"
    #saveAll("data/trap6test/trap6test_raw.csv", dataFolder, modelName)
    #saveSynchronyCertainty(dataFolder, modelName)
    #model = setupGraph(modelFolder, dataFolder, modelName, 10, 36, 4, 29, save=True, debug=1, progressBar=True)
    #model = addDelayFunctionToModel(modelFolder, dataFolder, modelName, "_None_None_None", "functions/colourDelayFunctionUniform.pickle", save=True, debug=1, progressBar=True, safeMode=False, loadedModel=model)
    #model = addDelayFunctionToModel(modelFolder, dataFolder, modelName, "_None_colourDelayFunctionUniform_None", "functions/maturationDelayFunctionNormal.pickle", save=True, debug=1, progressBar=True, safeMode=False, loadedModel=model)
    #model = addConjugationFunctionToModel(modelFolder, dataFolder, modelName, "_None_colourDelayFunctionUniform_maturationDelayFunctionNormal", "functions/contactWeightsBaseline.pickle", save=True, debug=1, progressBar=True, safeMode=False, loadedModel=model, saveConjugationFunction=True)
    #model, evidence, forwardLinks = get_evidence(modelFolder, dataFolder, modelName, "_contactWeightsFixedNaive_colourDelayFunctionUniform_maturationDelayFunctionNormal", save=True, debug=1, progressBar=True, loadedModel=model)
    queries = build_queries(modelFolder, dataFolder, modelName, "_contactWeightsFixedNaive_colourDelayFunctionUniform_maturationDelayFunctionNormal_contradictionsPruned", save=True, debug=1, progressBar=True, loadedModel=None)
