import os

# Immutable interface "constants"
def paramOrder():
    return ["zoneLength", "t", "tc", "beta", "x", "J"]
def outOrder(): 
    return ["muF", "muB", "D", "B", "A", "Dc", "Bc", "Ac"]

# Returns a list of dicts with values solving the equations associated with
# the parameters provided in paramList.
def executeRuns(paramList, paramFilePath, outputFilePath):
    # is paramList ok to use?
    if not validateParams(paramList):
        raise ValueError("paramList failed validation")
    # write data for fortran program
    numRuns = len(paramList)
    paramFile = open(paramFilePath, 'w')
    paramFile.write(str(numRuns) + "\n")
    for params in paramList:
        for key in paramOrder():
            paramFile.write(str(params[key]) + " ")
        paramFile.write("\n")
    paramFile.close()
    # run fortran program and wait for it to finish
    os.spawnl(os.P_WAIT, "main.out", "main.out", paramFilePath, outputFilePath)
    # read output
    outputList = []
    outputFile = open(outputFilePath, 'r')
    for line in outputFile.readlines():
        if line[0] == "!":
            outputList.append({})
            continue
        values = " ".split(line)
        lineDict = {}
        for i, key in enumerate(outOrder()):
            lineDict[key] = values[i]
        outputList.append(lineDict)
    return outputList

# All elements of paramList must be dicts containing each key in paramOrder.
# Return True if this holds; False if not.
def validateParams(paramList):
    for params in paramList:
        for key in paramOrder():
            if key not in params:
                return False
    return True
