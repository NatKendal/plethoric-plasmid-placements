import pickle

def mergeQueries(modelFolder, modelName, modelExtension, modelExtensionExtra, colourFunctionFile, save=True, debug=0, progressBar=False):
    if progressBar:
        import tqdm

    if debug >= 1:
        print("Loading queries.")
    with open(modelFolder + modelName + "_modeldata" + modelExtension + "_queries.pickle", "rb") as f:
        queries = pickle.load(f)

    if debug >= 1:
        print("Loading query point values.")
    with open(modelFolder + modelName + "_modeldata" + modelExtension + modelExtensionExtra + "_querypoints.pickle", "rb") as f:
        queryPointValues = pickle.load(f)

    if debug >= 1:
        print("Loading colour function.")
    with open(colourFunctionFile, "rb") as f:
        colourFunction = pickle.load(f)

    queryEvaluation = dict()
    
    if debug >= 1:
        print("Starting querypoint merge.")

    if progressBar:
        iterator = tqdm.tqdm(queries.keys())
    else:
        iterator = queries.keys()
    for query in iterator:
        if progressBar:
            iterator.set_description(desc="Working on query " + str(query))
        if debug >= 2:
            print("Working on query " + str(query))
        probability = 0.0
        queryTime = int(query.split("_")[1])
        for point in queries[query]["critical"]:
            if queryPointValues[point] == -1:
                probability = -1
                break
            else:
                probability += queryPointValues[point] * colourFunction.value(abs(queryTime - int(point.split("_")[1])))
        queryEvaluation[query] = probability

    if debug >= 1:
        print("Finished processing all queries.")

    if save:
        with open(modelFolder + modelName + "_modeldata" + modelExtension + modelExtensionExtra + "_queryEvaluation.pickle", "wb") as f:
            pickle.dump(queryEvaluation, f)

    return queryEvaluation
