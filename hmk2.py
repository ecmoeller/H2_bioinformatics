
import math
from operator import itemgetter 

import plotly
plotly.tools.set_credentials_file(username='moell162', api_key='XWeDQVjHzfvEN01WsE7p')

import plotly.plotly as py
import plotly.graph_objs as go

# Create random data with numpy
import numpy as np

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

# Compute sensitivity for each threshold
def sensitivityCompute(counts, k, drugMatrix):
    #TP / (TP + FN)
    
    #starting at k
    threshold = k
    array = []
    while threshold > 0:
        tp = 0
        fn = 0
        for key, value in counts.items():
            if(drugMatrix[int(key)] != "NA"):
                if(int(drugMatrix[int(key)]) == 1):
                    #True positive
                    if(int(value) >= threshold):
                        tp += 1
                    #False negative
                    else:
                        fn += 1
        sens = tp / (tp + fn)
        array.append(sens)
        threshold -= 1

    for i in array:
        print(i)

    return array

# Compute 1 - specificty for each threshold
def one_specificity(counts, k, drugMatrix):
    # 1 - (TN / (TN + FP))

    #starting at k
    threshold = k
    array = []
    while threshold > 0:
        tn = 0
        fp = 0
        for key, value in counts.items():
            if(drugMatrix[int(key)] != "NA"):
                if(int(drugMatrix[int(key)]) == 0):
                    #False positive
                    if(int(value) >= threshold):
                        fp += 1
                    #True negative
                    else:
                        tn += 1
        sens = 1 - (tn / (tn + fp))
        array.append(sens)
        threshold -= 1

    for i in array:
        print(i)

    return array

def linegraph(sensNum, specNum, name):
    # Create a trace
    knn = go.Scatter(
        x = specNum,
        y = sensNum,
        name = "knn"
    )
    xaxis = [0, .1, .2, .3, .4, .5, .6]
    yaxis = [0, .1, .2, .3, .4, .5, .6]
    random = go.Scatter(
        x = xaxis,
        y = yaxis,
        name = "random"
    )


    layout = go.Layout(
        title=go.layout.Title(
            text=name
            
        ),
        xaxis=go.layout.XAxis(
            title=go.layout.xaxis.Title(
                text='1-Specificity'
            )
        ),
        yaxis=go.layout.YAxis(
            title=go.layout.yaxis.Title(
                text='Sensitivity'
            )
        )
    )

    data = [knn, random]

    fig = go.Figure(data=data, layout=layout)
    py.plot(fig, filename=name)


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
    sensNum = sensitivityCompute(first5, 5, sensitivity[0])
    specNum = one_specificity(first5, 5, sensitivity[0])
    #Plot those points as well as y= x
    linegraph(sensNum, specNum, "Everolimus(mTOR) k=5")

    #Disulfiram(ALDH2)
    second5 = knncluster(5, newMatrix, sensitivity[1])
    sensNum = sensitivityCompute(second5, 5, sensitivity[1])
    specNum = one_specificity(second5, 5, sensitivity[1])
    linegraph(sensNum, specNum, "Disulfiram(ALDH2) k=5")

    #Methylglyoxol(Pyruvate)
    third5 = knncluster(5, newMatrix, sensitivity[2])
    sensNum = sensitivityCompute(third5, 5, sensitivity[2])
    specNum = one_specificity(third5, 5, sensitivity[2])
    linegraph(sensNum, specNum, "Methylglyoxol(Pyruvate) k=5")

    #Mebendazole(Tubulin)
    fourth5 = knncluster(5, newMatrix, sensitivity[3])
    sensNum = sensitivityCompute(fourth5, 5, sensitivity[3])
    specNum = one_specificity(fourth5, 5, sensitivity[3])
    linegraph(sensNum, specNum, "Mebendazole(Tubulin) k=5")

    #4-HC(DNA alkylator)
    fifth5 = knncluster(5, newMatrix, sensitivity[4])
    sensNum = sensitivityCompute(fifth5, 5, sensitivity[4])
    specNum = one_specificity(fifth5, 5, sensitivity[4])
    linegraph(sensNum, specNum, "4-HC(DNA alkylator) k=5")

    #QUESTION 3

if __name__ == '__main__':
    main()
