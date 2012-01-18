paramOrder = ["zoneLength", "t", "tc", "beta", "x", "J"]
outOrder = ["muF", "muB", "D", "B", "A", "Dc", "Bc", "Ac"]

# Returns a list of dicts with values solving the equations associated with
# the parameters provided in paramList.
def executeRuns(paramList, paramFilePath, outputFilePath):
    if not validateParams(paramList):
        raise ValueError("paramList failed validation")
    numRuns = len(parameterList)
    paramFile = open(paramFilePath, 'w')
    paramFile.write(str(numRums) + "\n")
    for params in paramList:
        for key in paramOrder:
            paramFile.write(str(params[key]) + " ")
        paramFile.write("\n")
    paramFile.close()
    # --- TODO: run fortran code here ---
    outputList = []
    outputFile = open(outputFilePath, 'r')
    for line in outputFile.readlines():
        values = " ".split(line)
        lineDict = {}
        for i, key in enumerate(outOrder):
            lineDict[key] = values[i]
        outputList.append(lineDict)
    return outputList

# TODO
# All elements of paramList must be dicts containing each key in paramOrder.
# Return True if this holds; False if not.
def validateParams(paramList):
    return True
