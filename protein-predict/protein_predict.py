﻿#project for predicting the protein deformation during protein-protein interaction
#co-conducted by Paul & Adam
import copy
import statistics
import os 
import time
import pprint
import numpy
proteinDir = "./res/allprot.fa"
labelDir = "./res/label2prot2pair"
allProtein = open(proteinDir,"r")
labelFile = open(labelDir,"r")

proteinSerial = []
proteinSeq = []
proteinLen=[]
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
        proteinLen.append(len(line)-1)
     

noZeroSeq = [[x for x in seq if x!=0] for seq in proteinSeq] #remove 0s in sequence

min = 100
minIndex = 0

for index in range(len(noZeroSeq[0])-3):     #need at least two element to compute stdev
    stdevSum = (statistics.stdev(noZeroSeq[0][:index+2]))+(statistics.stdev(noZeroSeq[0][index+2:]))
    if stdevSum<min:
        min = stdevSum
        minIndex = index+2


#filter

oneFilterSeq = copy.deepcopy(noZeroSeq)              #deep copy of list

# first filter gets rid of single pulses
for i in range(len(noZeroSeq)):
    removeIndex = []
    seqLen = len(noZeroSeq[i])
    if noZeroSeq[i][0]*noZeroSeq[i][1]<0:       #handle first element
        removeIndex.append(0)
    
    for j in range(seqLen-2):
        if noZeroSeq[i][j+1]*noZeroSeq[i][j] < 0 and noZeroSeq[i][j+1]*noZeroSeq[i][j+2] < 0:        #handle others
            removeIndex.append(j+1)

    if noZeroSeq[i][seqLen-2]*noZeroSeq[i][seqLen-1]<0:       #handle last element
        removeIndex.append(seqLen-1)
    

    for item in list(reversed(removeIndex)):        #store first filtered sequence
        del oneFilterSeq[i][item]
######end of first filter

# second filter gets rid of double pulses
twoFilterSeq = copy.deepcopy(oneFilterSeq)          #deep copy of list
for i in range(len(oneFilterSeq)):
    removeIndex = []
    seqLen = len(oneFilterSeq[i])
    if oneFilterSeq[i][0]*oneFilterSeq[i][1]>0 and oneFilterSeq[i][1]*oneFilterSeq[i][2]<0: #handle first two elements
        removeIndex.append(0)
   
    for j in range(seqLen-3):
        if oneFilterSeq[i][j+1]>0 and oneFilterSeq[i][j+2]>0 and oneFilterSeq[i][j]<0 and oneFilterSeq[i][j+3]<0:   #handle other elements
            removeIndex.append(j+1)
        elif oneFilterSeq[i][j+1]<0 and oneFilterSeq[i][j+2]<0 and oneFilterSeq[i][j]>0 and oneFilterSeq[i][j+3]>0:
            removeIndex.append(j+1)

    if oneFilterSeq[i][seqLen-1] * oneFilterSeq[i][seqLen-2]>0 and oneFilterSeq[i][seqLen-2]*oneFilterSeq[i][seqLen-3]<0: #handle last two elements
        removeIndex.append(seqLen-2)

    for item in list(reversed(removeIndex)):    #store second filterd sequence
        del twoFilterSeq[i][item+1]
        del twoFilterSeq[i][item]

#####end of second filter


#export file for excel analysis
excelFile1 = open("./data/proteinOneFilter.txt","w")   
excelFile2 = open("./data/proteinTwoFilter.txt","w")


for i in oneFilterSeq[0]:
    excelFile1.write(str(i)+"\n")
excelFile1.close()
for i in twoFilterSeq[0]:
    excelFile2.write(str(i)+"\n")
excelFile2.close()


temp=0      #存取每一個index前一個變數
eachAccum=[]  #每一列相同符號的數字個數相加
maxFourRatio=[] #前四個所佔濾波兩次後長度的比例
hydroFeature=[] #最後輸出二維list
for i in range(len(twoFilterSeq)):
    count=0 #算個數
    for j,item in  enumerate(twoFilterSeq[i]): #每一行的每一個胺基酸親水指數       
        if item*temp<0:  #切換數字正負號時進入           
            eachAccum.append(count*int(temp/abs(temp)))  #前面的count都是正，乘上數字正負號                              
            count=1
        else:
            count+=1
        temp=item
    eachAccum.append(count*int(temp/abs(temp))) #最後一組不會切換正負號所以直接加入count的結果沒關係

    sortedAccum=sorted(range(len(eachAccum)), key=lambda k: abs(eachAccum[k]), reverse=True)   #將所有數字轉正號排序方便取下前四名
    if len(sortedAccum)>=4:
        maxFourRatio.append((abs(eachAccum[sortedAccum[0]])+abs(eachAccum[sortedAccum[1]])+abs(eachAccum[sortedAccum[2]])+abs(eachAccum[sortedAccum[3]]))/len(twoFilterSeq[i]))
    else: 
        maxFourRatio.append(0)
    
    tempList=sorted(sortedAccum[0:4]) #將前四名順序再次排序
    hydroFeature.append([int(eachAccum[y]/abs(eachAccum[y])) for y in tempList]) #透過四個名次索引值取出後判斷正負號存入
    temp=0    #每行結束將temp初始化為0，才不會影響下一組正負號切換
    
    del eachAccum[:] #清空每行累加結果


accumLen=[]       #store each the amount of each segment            
for i in range(50,7000,50):    
    accumLen.append(len([ x for x in proteinLen if x<i and x>=i-50 ])) #count the amount of length in region i-50~i

excelFile3 = open("./data/sortedProteinLength.txt","w")   
for i in accumLen:
    excelFile3.write(str(i)+"\n")
excelFile3.close()





featureTable = {"1 1 1 1":1,"1 1 1 -1":2,"-1 1 1 1":2,"1 1 -1 1":3,"1 -1 1 1":3,"1 1 -1 -1":4,"-1 -1 1 1":4,"1 -1 1 -1":5,"-1 1 -1 1":5,\
                "1 -1 -1 1":6,"1 -1 -1 -1":7,"-1 -1 -1 1":7,"-1 1 1 -1":8,"-1 1 -1 -1":9,"-1 -1 1 -1":9,"-1 -1 -1 -1":10}


trainData = open("train-data","w")      #open train data
testData = open("test-data","w")        #open test data
specialData = open("./data/special-feature.txt","w")        #特別存取6 8 的結果
twoFilterLengthData = open("./data/twoFilterLengthData.txt","w")        #特別存取6 8 的結果
ratioData = open("./data/ratioData.txt","w") #特別存取6 8 的結果
def SepLen(x):     #transform length to feature
        if x<150:
            return 1
        elif x<250:
            return 2
        elif x<350:
            return 3
        else: 
            return 4


for i,line in enumerate(labelFile):     #挑出長的觀察
    index = proteinSerial.index(line.split()[1]+"\n")   #find index of protein in protein serial
    if len(twoFilterSeq[index])>400: #這挑出大於400的看
        twoFilterLengthData.write(line.split()[1]+" ")
        twoFilterLengthData.write(line.split()[0])
        twoFilterLengthData.write(" "+str(len(twoFilterSeq[index])))
        twoFilterLengthData.write("\n")

for i,line in enumerate(labelFile):     #write train data
    if i< 500:
        index = proteinSerial.index(line.split()[1]+"\n")   #find index of protein in protein serial
        if len(hydroFeature[index]) == 4:  
           
            temp = " ".join(map(str,hydroFeature[index]))   #transform list to str for search in dict
            feature2=str(SepLen(proteinLen[index]))
            feature = str(featureTable[temp])               #search feature table for feature
           
            if feature=="6": #將feature1特徵為6的紀錄                
                specialData.write(line.split()[1]+" ")
                specialData.write(line.split()[0])
                for i in twoFilterSeq[index]:
                    specialData.write(" "+str(i))
                specialData.write("\n")  
              
            if line.split()[0]=="1":
                trainData.write(line.split()[0] + " 1:" + feature + " 2:" + feature2+"\n")  #write four times if label equals to 1
                trainData.write(line.split()[0] + " 1:" + feature + " 2:" + feature2+"\n")
                trainData.write(line.split()[0] + " 1:" + feature + " 2:" + feature2+"\n")
                trainData.write(line.split()[0] + " 1:" + feature + " 2:" + feature2+"\n")
            else:
                trainData.write(line.split()[0] + " 1:" + feature + " 2:" + feature2+"\n")
        else:
            if line.split()[0]=="1":
                trainData.write(line.split()[0] + " 1:11" +" 2:" + feature2+"\n")       #deal with short proteins
                trainData.write(line.split()[0] + " 1:11" +" 2:" + feature2+"\n")
                trainData.write(line.split()[0] + " 1:11" +" 2:" + feature2+"\n")
                trainData.write(line.split()[0] + " 1:11" +" 2:" + feature2+"\n")
            else:
                trainData.write(line.split()[0] + " 1:11" +" 2:" + feature2+"\n")
        

labelFile.close()
labelFile = open(labelDir,"r")          #re-open file
for i,line in enumerate(labelFile):     #write test data
    if i>= 500:
        index = proteinSerial.index(line.split()[1]+"\n")
        if len(hydroFeature[index]) == 4:
            temp = " ".join(map(str,hydroFeature[index]))
            feature = str(featureTable[temp])
            
            if feature=="6":                
                specialData.write(line.split()[1]+" ")
                specialData.write(line.split()[0])
                for i in twoFilterSeq[index]:
                    specialData.write(" "+str(i))
                specialData.write("\n")  

            testData.write(line.split()[0] + " 1:" + feature)
        else:
            testData.write(line.split()[0] + " 1:11")

        feature2=str(SepLen(proteinLen[index]))
        testData.write(" 2:" + feature2+"\n")

labelFile.close()
labelFile = open(labelDir,"r")          #re-open file
for i,line in enumerate(labelFile):     #吐出所有6 8 在整個長度的比例
    index = proteinSerial.index(line.split()[1]+"\n")
    if len(hydroFeature[index]) == 4:
        temp = " ".join(map(str,hydroFeature[index]))
        feature = str(featureTable[temp])            
        if feature=="6" or feature=="8":                
           ratioData.write(line.split()[1]+" ")
           ratioData.write(line.split()[0]+" ")
           ratioData.write(str(feature)+" ")
           ratioData.write(str(maxFourRatio[index])+" ")           
           ratioData.write("\n")  

checkva100_6_8_maxFourRatio = open("./data/va100_6_8_maxFourRatio.txt","w")
check50_6_8_maxFourRatio = open("./data/check50_6_8_maxFourRatio.txt","w")
check50_6 = open("check50-6","w")
check50_8 = open("check50-8","w")
vaFile = open("./res/vate100-lst","r")
predictData = open("predict-data","w")
for line in vaFile:     #take VA50 and write predict-data
    index = proteinSerial.index(line.split()[0] + "\n")
    if len(hydroFeature[index]) == 4:
        temp = " ".join(map(str,hydroFeature[index]))
        feature = str(featureTable[temp])
        if feature=="6" or feature=="8":                
           checkva100_6_8_maxFourRatio.write(line.split()[0]+" ")
           checkva100_6_8_maxFourRatio.write(str(feature)+" ")
           checkva100_6_8_maxFourRatio.write(str(maxFourRatio[index])+" ")           
           for i in twoFilterSeq[index]:
              checkva100_6_8_maxFourRatio.write(" "+str(i))
           checkva100_6_8_maxFourRatio.write("\n")    
        if feature =="6":
            tempp = " ".join(map(str,twoFilterSeq[index]))
            check50_6.write(tempp + "\n")
        if feature =="8":
            tempp = " ".join(map(str,twoFilterSeq[index]))
            check50_8.write(tempp + "\n")
        predictData.write("0" + " 1:" + feature)
    else:
        predictData.write("0" + " 1:11")

    feature2=str(SepLen(proteinLen[index]))
    predictData.write(" 2:" + feature2+"\n")
check50_6.close()
check50_8.close()
forGraphic = open("./data/forGraphic.txt","w")
labelFile.close()
labelFile = open(labelDir,"r")          #re-open file
numOfPos = 0
numOfNeg = 0
for line in labelFile:
    index = proteinSerial.index(line.split()[1]+"\n")
    if numOfPos != 10 and line.split()[0] == "1":
        numOfPos += 1
        temp = " ".join(map(str,twoFilterSeq[index]))
        forGraphic.write("1 " + temp + "\n")
    elif numOfNeg != 10 and line.split()[0] == "-1":
        numOfNeg += 1
        temp = " ".join(map(str,twoFilterSeq[index]))
        forGraphic.write("-1 " + temp + "\n")

predictData.close()
vaFile.close()
trainData.close()
testData.close()
labelFile.close()
allProtein.close()
specialData.close()
twoFilterLengthData.close()
ratioData.close()
check50_6_8_maxFourRatio.close()
checkva100_6_8_maxFourRatio.close()