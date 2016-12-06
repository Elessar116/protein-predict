#project for predicting the protein deformation during protein-protein interaction
#co-conducted by Paul & Adam
import copy
import statistics

proteinDir = "./res/allprot.fa"
labelDir = "./res/label2prot2pair"
allProtein = open(proteinDir,"r")
labelFile = open(labelDir,"r")

proteinSerial = []
proteinSeq = []
hydroList = [("A",1.8),("R",-4.5),("N",-3.5),("D",-3.5),("C",2.5),("E",-3.5),("Q",-3.5),("G",-0.4),("H",-3.2),\
    ("I",4.5),("L",3.8),("K",-3.9),("M",1.9),("F",2.8),("P",-1.6),("S",-0.8),("T",-0.7),("W",-0.9),("Y",-1.3),("V",4.2),("X",0)]
hydroListTest = [("A",1.8),("R",-4.5),("N",-3.5),("D",-3.5),("C",2.5),("E",-3.5),("Q",-3.5),("G",0),("H",-3.2),\
    ("I",4.5),("L",3.8),("K",-3.9),("M",1.9),("F",2.8),("P",-1.6),("S",0),("T",0),("W",0),("Y",-1.3),("V",4.2),("X",0)]
hydroDict = (dict(hydroListTest))


def SeqToHydroList(str):    #function to transfrom amino acid seq into hydropathy index list
    temp = []
    for letter in str:
        if letter !="\n":
            temp.append(hydroDict[letter])
    return temp

for line in allProtein:
    if line[0] == ">":
        proteinSerial.append(line[1:])  #remove ">" and store serial number in list
    else:
        proteinSeq.append(SeqToHydroList(line))         #store sequence in list

noZeroSeq = [[x for x in seq if x!=0] for seq in proteinSeq] #remove 0s in sequence
#for seq in proteinSeq:
#    noZeroSeq.append([x for x in seq if x!=0]) 

#print(statistics.stdev(proteinSeq[0]))
min = 100
minIndex = 0

for index in range(len(noZeroSeq[0])-3):     #need at least two element to compute stdev
    stdevSum = (statistics.stdev(noZeroSeq[0][:index+2]))+(statistics.stdev(noZeroSeq[0][index+2:]))
    if stdevSum<min:
        min = stdevSum
        minIndex = index+2


#filter

isPositive = noZeroSeq[0][0]>0

oneFilterSeq = copy.deepcopy(noZeroSeq)              #deep copy of list

# first filter gets rid of single pulses
for i in range(len(noZeroSeq)):
    removeIndex = []
    seqLen = len(noZeroSeq[i])
    if noZeroSeq[i][0]<0 and noZeroSeq[i][1]>0:       #handle first element
        removeIndex.append(0)
    elif noZeroSeq[i][0]>0 and noZeroSeq[i][1]<0:
        removeIndex.append(0)

    for j in range(seqLen-2):
        if noZeroSeq[i][j+1] > 0 and noZeroSeq[i][j] < 0 and noZeroSeq[i][j+2] < 0:        #handle others
            removeIndex.append(j+1)
        elif noZeroSeq[i][j+1] < 0 and noZeroSeq[i][j] > 0 and noZeroSeq[i][j+2] > 0:
            removeIndex.append(j+1)

    if noZeroSeq[i][seqLen-2]<0 and noZeroSeq[i][seqLen-1]>0:       #handle last element
        removeIndex.append(seqLen-1)
    elif noZeroSeq[i][seqLen-2]>0 and noZeroSeq[i][seqLen-1]<0:
            removeIndex.append(seqLen-1)

    for item in list(reversed(removeIndex)):        #store first filtered sequence
        del oneFilterSeq[i][item]
######end of first filter

# second filter gets rid of double pulses
twoFilterSeq = copy.deepcopy(oneFilterSeq)          #deep copy of list
for i in range(len(oneFilterSeq)):
    removeIndex = []
    seqLen = len(oneFilterSeq[i])
    if oneFilterSeq[i][0]>0 and oneFilterSeq[i][1]>0 and oneFilterSeq[i][2]<0: #handle first two elements
        removeIndex.append(0)
    elif oneFilterSeq[i][0]<0 and oneFilterSeq[i][1]<0 and oneFilterSeq[i][2]>0:
        removeIndex.append(0)
    for j in range(seqLen-3):
        if oneFilterSeq[i][j+1]>0 and oneFilterSeq[i][j+2]>0 and oneFilterSeq[i][j]<0 and oneFilterSeq[i][j+3]<0:   #handle other elements
            removeIndex.append(j+1)
        elif oneFilterSeq[i][j+1]<0 and oneFilterSeq[i][j+2]<0 and oneFilterSeq[i][j]>0 and oneFilterSeq[i][j+3]>0:
            removeIndex.append(j+1)
    if oneFilterSeq[i][seqLen-1]>0 and oneFilterSeq[i][seqLen-2]>0 and oneFilterSeq[i][seqLen-3]<0: #handle last two elements
        removeIndex.append(seqLen-2)
    elif oneFilterSeq[i][seqLen-1]<0 and oneFilterSeq[i][seqLen-2]<0 and oneFilterSeq[i][seqLen-3]>0:
        removeIndex.append(seqLen-2)

    for item in list(reversed(removeIndex)):    #store second filterd sequence
        del twoFilterSeq[i][item+1]
        del twoFilterSeq[i][item]

#####end of second filter


#export file for excel analysis
excelFile1 = open("proteinOneFilter.txt","w")   
excelFile2 = open("proteinTwoFilter.txt","w")

for i in oneFilterSeq[0]:
    excelFile1.write(str(i)+"\n")
excelFile1.close()
for i in twoFilterSeq[0]:
    excelFile2.write(str(i)+"\n")
excelFile2.close()

#helllo101010101101010