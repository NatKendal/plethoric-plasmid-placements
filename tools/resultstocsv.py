import math
import pickle

if __name__ == "__main__":
    with open("filenames.txt", "r") as f:
        filenames = f.readlines()
    for i in range(len(filenames)):
        filenames[i] = filenames[i].strip()

    traps = dict()
    for filename in filenames:
        if filename.split("_")[0] not in traps:
            traps[filename.split("_")[0]] = []
        traps[filename.split("_")[0]].append(filename)

    mainFolder = "Experiment3Results/"

    for trap in traps:
        columnHeaders = ["query"]
        columns = []
        nonConQueries = None
        nonConjugateColumns = []
        queries = None
        for filename in sorted([f for f in traps[trap] if f.split("_")[-1] == "fullQueries.pickle"], key=lambda x: int(x.split("_")[4][5:])):
            columnHeaders.append("_".join(filename.split("_")[:5]))
            with open(mainFolder + filename, "rb") as f:
                fullQueryEval = pickle.load(f)
            with open(mainFolder + "_".join(filename.split("_")[:9]) + "_nonConjugateQueries.pickle", "rb") as f:
                ncQueryEval = pickle.load(f)

            if queries == None:
                queries = sorted(fullQueryEval.keys(), key=lambda x: int(x.split("_")[1]))
            else:
                if queries != sorted(fullQueryEval.keys(), key=lambda x: int(x.split("_")[1])):
                    print("Queries didn't match!")
                    print(trap, filename)
            newColumn = []
            for query in queries:
                newColumn.append(fullQueryEval[query])
            columns.append(newColumn)

            if nonConQueries == None:
                nonConQueries = sorted(ncQueryEval.keys(), key=lambda x: int(x.split("_")[1]))
            else:
                if nonConQueries != sorted(ncQueryEval.keys(), key=lambda x: int(x.split("_")[1])):
                    print("Non con queries didn't match!")
                    print(trap, filename)
            ncColumn = []
            for query in nonConQueries:
                ncColumn.append(ncQueryEval[query])
            nonConjugateColumns.append(ncColumn)

        with open(trap + "_conjugateResults.csv", "w") as f:
            f.write(",".join(columnHeaders) + "\n")
            for i in range(len(queries)):
                f.write(queries[i] + "," + ",".join([str(columns[j][i]) for j in range(len(columns))]) + "\n")
            f.flush()

        with open(trap + "_nonConjugateResults.csv", "w") as f:
            f.write(",".join(columnHeaders) + "\n")
            for i in range(len(nonConQueries)):
                f.write(nonConQueries[i] + "," + ",".join([str(nonConjugateColumns[j][i]) for j in range(len(nonConjugateColumns))]) + "\n")
            f.flush()






