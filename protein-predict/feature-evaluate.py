testData = open("test-data","r")

posCount = 0
negCount = 0
for line in testData:
    if float(line.split()[3][2:])>300 and int(line.split()[0])==1:
        posCount += 1
    elif float(line.split()[3][2:])>300 and int(line.split()[0])==-1:
        negCount += 1
print(posCount,negCount)

trainData = open("train-data","r")
posCount2 = 0
negCount2 = 0
for line in trainData:
    if float(line.split()[3][2:])>300 and int(line.split()[0])==1:
        posCount2 += 1
    elif float(line.split()[3][2:])>300 and int(line.split()[0])==-1:
        negCount2 += 1
print(int(posCount2/5),negCount2)