# -*- coding: utf-8 -*-
import numpy as np
import math
def isTarget(ID):
    if int(ID)<1e10 or ID in ["9999999999"]:
        return True
    else :
        return False
    
def antiprojection(starksi, starita, Tan_A, Tan_D):
    import numpy
    anti_reflect_right = starksi / (numpy.cos(numpy.radians(Tan_D)) - starita * numpy.sin(numpy.radians(Tan_D)))
    alpha11 = numpy.arctan(anti_reflect_right) + numpy.radians(Tan_A)
    alpha22 = alpha11 * 180.0 / numpy.pi

    delta11_rightleft = (numpy.sin(numpy.radians(Tan_D)) + starita * numpy.cos(numpy.radians(Tan_D))) / (
                numpy.cos(numpy.radians(Tan_D)) - starita * numpy.sin(numpy.radians(Tan_D)))
    delta11_rightright = numpy.cos(alpha11 - numpy.radians(Tan_A))
    delta11 = numpy.arctan(np.array(delta11_rightleft) * np.array(delta11_rightright)) * 180.0 / numpy.pi
    return alpha22, delta11
 
    
class plsq:
    def __init__(self):
        self.x = []
        self.fun = [] 
# In[3]:

'Generate a star exclusion dictionary based on the sorted statListSorted, and return meanStdList as the median standard deviation within each magnitude range'
def genExcludeStarDict(statListSorted,excludeStarDict):
    import statistics as stat
    tempStdList=[]
    meanStdList=[]
    if len(statListSorted)==0:
        return []
    magThreshold=int(float(statListSorted[0][0]))+1
    for starStd in statListSorted:
        if float(starStd[0])<magThreshold :
            tempStdList.append([float(starStd[1]),starStd[2],starStd[0]])
        if float(starStd[0])>=magThreshold or starStd==statListSorted[-1]:
            
            inliersIndexs,deleteCount=DeleteOutliers(np.array(tempStdList,dtype=float)[:,0],3)
            
            for deleteItem in np.array(tempStdList)[~inliersIndexs]:
                if not isTarget(deleteItem[1]):
                    excludeStarDict[deleteItem[1]]=[deleteItem[2],deleteItem[0]] 
            if len(tempStdList)>0:
                meanStdList.append([magThreshold-0.5,stat.median(np.array(tempStdList,dtype=float)[inliersIndexs][:,0])])
                tempStdList=[[float(starStd[1]),starStd[2],starStd[0]]]
            else:
                tempStdList.append([float(starStd[1]),starStd[2],starStd[0]])
            magThreshold+=1
        
    print(str(len(excludeStarDict))+" Stars that will be discarded in the next iteration")
    return meanStdList    
  
 
def genLeftMatrix(xArray,yArray,model,weights,magArray=[],magOrder=0) :
    if model==1:
        return np.mat(weights*[np.append(np.ones(int(len(xArray)/2)),np.zeros(int(len(xArray)/2)))
                               ,xArray,yArray,np.append(np.zeros(int(len(xArray)/2)),np.ones(int(len(xArray)/2)))])    
    else:
        generatedLeftList=genPolynomialLeftMatrix(xArray,yArray,model-1)
        if magArray!=[]:
            for magTermNum in range(magOrder):
                generatedLeftList.append(np.array(magArray)**(magTermNum+1)) 
        return np.mat(weights*generatedLeftList)
 
def genPolynomialLeftMatrix(xArray,yArray,order):
    if order>1:
        lowOrderMatrix=genPolynomialLeftMatrix(xArray,yArray,order-1)
        for itemI in range(order+1):
            lowOrderMatrix.append(xArray**itemI*yArray**(order-itemI))
            # print("+x^",itemI,"*y^",(order-itemI),end='')
        # print("")
        return lowOrderMatrix
    else:
        # print("1+x+y",end='')
        return [np.ones(len(xArray)),xArray,yArray]

def model_Correct_GD(allFrameStarList,xCoeff,yCoeff):
    import numpy as np 
    import mylib.ioFile_ReadMoreCols as ioFile
    import re 
    order=int((2*len(xCoeff))**0.5) -1
    tempAllFrameStarList=[]   
    for curFrameStarInfos in allFrameStarList: 
        newCurFrameStars=[]
        for num,starInfo in enumerate(curFrameStarInfos): 
            distoredX=float(starInfo.x)
            distoredY=float(starInfo.y)
      
            newX=np.array(np.mat(genPolynomialLeftMatrix(np.array([distoredX]),np.array([distoredY]),order)).T *np.mat(xCoeff).T)[0][0]
            newY=np.array(np.mat(genPolynomialLeftMatrix(np.array([distoredX]),np.array([distoredY]),order)).T *np.mat(yCoeff).T)[0][0]
         
            newCurFrameStars.append(ioFile.StarInfo(*re.split(r'\s+',  '{ID} {newX} {newY} {RA} {DE} {Mag} {pmRA} {pmDE}  {plx} {rv}'.format(**starInfo._asdict(),
                                                             newX= newX,newY= newY ))))
        tempAllFrameStarList.append(newCurFrameStars)
    
    allFrameStarList[:]=tempAllFrameStarList 
    return allFrameStarList 
   
def julian_date( year, month, day, hour):
#    jd12h =  day - 32075  + 1461  * ( year + 4800  + ( month - 14 ) // 12 ) // 4 \
#		+ 367  * ( month - 2  - ( month - 14 ) // 12  * 12 ) // 12  \
#        - 3  * (( year + 4900  + ( month - 14 ) // 12 )// 100 ) // 4   
    year=float(year)
    month=float(month)
    day=float(day)
    hour=float(hour)
    jd12h=367*year-7*( year+(month+9)//12)//4-3*((year+(month-9)//7)//100+1)//4+275*month//9+day+1721029
    tjd =  jd12h - 0.5 + hour / 24 
    return tjd 
   

def DeleteOutliers(inputArray,timesSigma,isDebug=False): 
    inliersIndexs=np.ones(len(inputArray),dtype=bool)
    deleteCount=0
    if len(inputArray)<2:
        return inliersIndexs,deleteCount
    meanVal=inputArray.mean()
    stdVal=inputArray.std()
    listLen=len(inputArray)
    for index in range(listLen):
        if abs(inputArray[listLen-1-index]-meanVal)>timesSigma*stdVal:
            inliersIndexs[listLen-1-index]=False            
            deleteCount+=1 
    if deleteCount>0 and isDebug: 
        print("The number of outliers thrown by DeleteOutliers:{}".format(deleteCount))
        print("listLen=",listLen) 
        print("stdVal=",stdVal)
        print("meanVal=",meanVal) 
        print("inputArray=",inputArray)
    return inliersIndexs,deleteCount
 