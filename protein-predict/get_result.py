resultFile = open("./predict-result","r")
cleanResult = open("./clean-result", "w")
next(resultFile)
for line in resultFile:
    cleanResult.write(line.split()[0] + ",")
resultFile.close()
cleanResult.close()