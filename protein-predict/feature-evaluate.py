testData = open("test-data","r")

posCount = 0
negCount = 0
for feature in range(1,21):
    print(str(feature))
    for line in testData:
        #print(line.split()[5][2:])
        if int(line.split()[5][2:])==feature and int(line.split()[0])==1:
            posCount += 1
        elif int(line.split()[5][2:])==feature and int(line.split()[0])==-1:
            negCount += 1
    print(posCount,negCount)
    

    trainData = open("train-data","r")
    posCount2 = 0
    negCount2 = 0
    for line in trainData:
        if int(line.split()[5][2:])==feature and int(line.split()[0])==1:
            posCount2 += 1
        elif int(line.split()[5][2:])==feature and int(line.split()[0])==-1:
            negCount2 += 1
    print(int(posCount2/5),negCount2)
    print("\n")
testData.close()
trainData.close()