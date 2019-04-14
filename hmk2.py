
import math
from operator import itemgetter 

def transpose(m):
    return [*zip(*m)]

def pearson(input1, input2):
    #Computing Pearson similarity coefficient
    # The vectors are in the correct format

    samplesize = len(input2)

    xsum = 0
    ysum = 0
    xysum = 0
    xsquaresum = 0
    ysquaresum = 0
    for i in range(len(input1)):
        x = float(input1[i])
        y = float(input2[i])
        xy = x * y
        xsquare = x * x
        ysquare = y * y
        xsum += x
        ysum += y
        xysum += xy
        xsquaresum += xsquare
        ysquaresum += ysquare

    numerator = (samplesize * xysum) - (xsum * ysum)
    denominator = math.sqrt((samplesize * xsquaresum - xsum * xsum) * (samplesize * ysquaresum - ysum * ysum))
    return (numerator / denominator)
    
#   You should create a classifier for each drug that takes a cell line expression 
# profile as input and produces a score that predicts whether it is sensitive or resistant to the given drug.
def knncluster(numOfClusters, input, drugMatrix):
    #input should be newMatrix, which is the one that was transposed 

    
    ##Computes the counts of sensitivity for top k most similiar cell lines
    finalDict = {}
    i = 1
    while i < len(input):
        #Compute pearson coefficients for all cell lines except for the last
        pearsonDict = {} # should be a dictionary associated with j
        j = 1
        
        while j < len(input):
            if(i != j):
                pearsonDict[str(j)] = pearson(input[j], input[i])
            j += 1
         
        index = 0
        count = 0
        #Sort and get top five similar cell lines
        #Compute number of cell lines that predicted sensitivity divided by k and i (cell line num) and add to dictionary
        for key, value in sorted(pearsonDict.items(), key = itemgetter(1), reverse = True):
            if(index < 5 and drugMatrix[int(key)] != "NA"):
                
                if(int(drugMatrix[int(key)]) == 1):
                    count += 1
                
            else:
                break
            index += 1

        finalDict[str(i)] = count
        print("This is count", count)
        i += 1

    ##Sort dictionary by predicted sensitivity nums
    print("Sorted final dictionary counts") 
    for key, value in sorted(finalDict.items(), key = itemgetter(1), reverse = True):
        print(key, value)

    
    return finalDict 

def sensitivity(counts):

    return 0

def one_specificity(counts):

def main():

    #Reading in data from txt file
    with open('DREAM_data.txt','r') as f:
        data = f.readlines()

    count = 0
    cellidsline = []
    sensitivity = []
    geneExpression = []
    for line in data:
        line = line.strip()
        if (count == 0):
            cellidsline = line.split('\t')
        if(count >= 1 and count <= 5):
            sensVals = line.split('\t')
            sensitivity.append(sensVals)
        if(count > 5):
            vals = line.split('\t')
            geneExpression.append(vals)

        count += 1



    newMatrix = transpose(geneExpression)
    print("Second row of geneExpression")
    for i in range(0,5):
        print(newMatrix[1][i], end =" ")

    #QUESTION 2
    #Everolimus(mTOR): first drug
    first5 = knncluster(5, newMatrix, sensitivity[0])
    #Compute sensitivity and 1-specificity for each threshold
    sensNum = sensitivity(first5)
    specNum = one_specificity(first5)
    #Plot those points as well as y= x

    #Disulfiram(ALDH2)
    second5 = knncluster(5, newMatrix, sensitivity[1])

    #Methylglyoxol(Pyruvate)
    third5 = knncluster(5, newMatrix, sensitivity[2])

    #Mebendazole(Tubulin)
    fourth5 = knncluster(5, newMatrix, sensitivity[3])

    #4-HC(DNA alkylator)
    fifth5 = knncluster(5, newMatrix, sensitivity[4])



    #QUESTION 3

if __name__ == '__main__':
    main()
