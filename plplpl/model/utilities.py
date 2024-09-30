import shutil

def copyModelData(modelFolder, modelName, modelExtension1, modelExtension2):
    shutil.copyfile(modelFolder+modelName+"_modeldata"+modelExtension1+"_criticalSegments.pickle", modelFolder+modelName+"_modeldata"+modelExtension2+"_criticalSegments.pickle")
    shutil.copyfile(modelFolder+modelName+"_modeldata"+modelExtension1+"_evidence.pickle", modelFolder+modelName+"_modeldata"+modelExtension2+"_evidence.pickle")
    shutil.copyfile(modelFolder+modelName+"_modeldata"+modelExtension1+"_queries.pickle", modelFolder+modelName+"_modeldata"+modelExtension2+"_queries.pickle")
