# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
""" 

import sys
sys.path.append("..")
import re,os  
import math
import GDSolver.mylib.myLeastSquare as mylsq
import numpy as np 
import GDSolver.mylib.ioFile_ReadMoreCols  as ioFile 
import statistics as stat 
import scipy
from scipy import linalg 
from matplotlib import pyplot as plt  
from  GDSolver.mylib.utils import model_Correct_GD,genLeftMatrix,genPolynomialLeftMatrix,plsq  
import GDSolver.mylib.precisionStatistic as MuApp
import GDSolver.mylib.plotStatistics as pltS



inputFileName='as_obs_20110228_gsc2038_G2_J2000_paper.txt' 
# inputFileName='test.txt'  

'Matched file names that need to correct geometric distortions'
correctedFileName=inputFileName 


'image size'
frameWidth=2048 
frameHeight=2048     

'The order of geometric distortions contained in the image'
model=4

ITERATE_TIMES= 3 #usually 2-4
BinCorrect=1 
x_partition=16
y_partition=16
GDCModel=model+2

perrx=0
perry=0

plt.rcParams['figure.figsize'] = (6.0, 6)
plt.rcParams['figure.dpi'] = 80
fontScale=0.8

"scale of the plotted quiver"
quiverScale=  550
'order of the distortion after this GDC process' 
correctOrder=0
 
 
NO_STAR_MODIFY=True
NO_STAR_REJECT=False
plotPrecisionStat=True

 

PAPER_FIG_TITLE="" 


'Do we use RANSAC to remove outliers for additional reduction operations'
USE_RANSAC_IN_REDUCTION=False
Residual_Threshold=10  

'The magnitude range read and processed by the current process'
MIN_PROCCESSING_MAG=1
MAX_PROCCESSING_MAG=30
  

bWriteFittedResult=True
bWriteStatResult=True
bRefineResult=False
           
if bRefineResult:
    outputFileName="refined_model"+str(model)+"_"+inputFileName
else:
    outputFileName="model"+str(model)+"_"+inputFileName


fittedFileName="iterT"+str(ITERATE_TIMES)+"_"+inputFileName[:inputFileName.find('.txt')]+"_"+str(model*(model+1))+"C.txt"

bSeparateStars=False

bOutputHead=True 

def DeleteLargeStd(inputArray,timesMean): 
    inliersIndexs=np.ones(len(inputArray),dtype=bool)
    deleteCount=0
    if len(inputArray)<2:
        return inliersIndexs,deleteCount
    meanVal=inputArray.mean()
    listLen=len(inputArray)
    for index in range(listLen):
        if inputArray[listLen-1-index]>timesMean*meanVal:
            inliersIndexs[listLen-1-index]=False            
            deleteCount+=1 
    return inliersIndexs,deleteCount

def EquatorialToStandard_Xi(RA, DE,meanRA,meanDE):
    RA=RA*math.pi/180
    DE=DE*math.pi/180
    meanRA=meanRA*math.pi/180
    meanDE=meanDE*math.pi/180
    return math.cos(DE) * math.sin(RA - meanRA) /(math.sin(DE) * math.sin(meanDE) + math.cos(DE) * math.cos(meanDE) * math.cos(RA - meanRA)) 

def EquatorialToStandard_Eta(RA, DE,meanRA,meanDE):
    RA=RA*math.pi/180
    DE=DE*math.pi/180
    meanRA=meanRA*math.pi/180
    meanDE=meanDE*math.pi/180 
    return (math.sin(DE) * math.cos(meanDE) - math.cos(DE) * math.sin(meanDE) * math.cos(RA - meanRA)) /\
         (math.sin(DE) * math.sin(meanDE) + math.cos(DE) * math.cos(meanDE) * math.cos(RA - meanRA)) 

def findBadIndexs(dataArray): 
    badDataIndexList=[] 
    deleteCount=1
    originArray=np.ones(len(dataArray))
    originArray[:]=dataArray
    while deleteCount>0:
        deleteCount=0
        dataStd=np.std(dataArray)
        dataMean=np.mean(dataArray)
        for num in range(len(dataArray)):
            if abs(dataArray[num]-dataMean)-3*dataStd>0:
                for oriNum in range(len(originArray)):
                    
                    if abs(dataArray[num]-originArray[oriNum])<0.0001:
                        badDataIndexList.append(oriNum)
                deleteCount+=1
        
        
        tempDataArray=np.ones(len(dataArray)-deleteCount)
        newArrayCount=0
        for num in range(len(dataArray)):
            if num not in badDataIndexList:
                tempDataArray[newArrayCount]= originArray[num] 
                newArrayCount+=1
        dataArray=tempDataArray
    return badDataIndexList

def deleteFrameByIndex( allLists,badFrameIndex):
    
    for curList in allLists:
        list.sort(badFrameIndex,reverse=True) 
        for index in badFrameIndex:
            curList.pop(index)  
        
        
'y = A2 + (A1-A2)/(1 + exp((x-x0)/dx))'
def fittingFunc(p,x):
    return   p[1] + (p[0]-p[1])/(1 + math.e**((x-p[2])/p[3]))

def fittingFuncResidual(p,x,obj):
    return   fittingFunc(p,x)- obj           
    
 
def refineFrameByCosts(allLists,model):
    fittingResultList=allLists[2]
    
    costArray=np.ones(len(fittingResultList))
    
    for num,plsqs in enumerate(fittingResultList):
        if model==1:
            costArray[num]=plsqs.cost
        else:
            
            costArray[num]=(plsqs[0].cost**2+plsqs[1].cost**2)**0.5
    badFrameIndex=findBadIndexs(costArray)
    deleteFrameByIndex(allLists,badFrameIndex)
    
    
    xLinearItemList=[]
    yLinearItemList=[]
    
    for plsqs in fittingResultList:
        if model==1:
            xLinearItemList.append(plsqs.x[1])
            yLinearItemList.append(plsqs.x[2])
        else: 
            xLinearItemList.append(plsqs[0].x[1])
            yLinearItemList.append(plsqs[0].x[2])
    
    badFrameIndex=findBadIndexs(xLinearItemList)
    
    deleteFrameByIndex(allLists,badFrameIndex)
    badFrameIndex=findBadIndexs(yLinearItemList)
    deleteFrameByIndex(allLists,badFrameIndex)


def solvePlateConstants(inputFileName, model , plsqPrecisionCurve, inputFileNamePrefix,isFitMagTerm=False,isChangeTransform=False
                         ,minMag=-9999,maxMag=9999):  
    allFrameStarList=[]
    allHeadInfoList=[] 
    if inputFileNamePrefix=="": 
        if(os.path.exists(os.path.dirname(__file__)+'/matchedFiles/'+inputFileName)):
            filePath=os.path.dirname(__file__)+"/matchedFiles/"
        elif(os.path.exists(os.path.dirname(__file__)+'/modifiedFiles/'+inputFileName)):
            filePath=os.path.dirname(__file__)+"/modifiedFiles/"
        else:
            print("The specified matching file was not found",inputFileName)
            return
    else:
        filePath=os.path.dirname(__file__)+"/tempFiles/"
                
    excludeReadStarDict=dict()    
    if(os.path.exists('HomogenizeData/excluedeStar/'+"exclude_"+inputFileName)):
        with open('HomogenizeData/excluedeStar/'+"exclude_"+inputFileName) as excludeStars:
            for line in excludeStars:
                line=line.strip() 
                fields= re.split(r'\s+',line) 
                excludeReadStarDict[fields[0]]=1
    
    ioFile.read_mathcedStar(allHeadInfoList,allFrameStarList, inputFileNamePrefix+inputFileName,filePath,minMag=minMag,maxMag=maxMag
                            ,excludeReadStarDict=excludeReadStarDict) 
    fittingResultList=[]
    allKsiList=[]
    allEtaList=[]
    deleteFrameList=[]
    picNum=1 
    allMagList=[]
    allWeightedResRAList=[]
    allWeightedResDEList=[] 
    allME1RAList=[]
    allME1DEList=[] 
    allFittingErrRAList=[]
    allFittingErrDEList=[] 
    allGDCoeffX=[]
    allGDCoeffY=[]
    
    allMeanGDCCoeffX=np.zeros(int(GDCModel*(GDCModel+1)/2))
    allMeanGDCCoeffY=np.zeros(int(GDCModel*(GDCModel+1)/2))
    
    allMeanGDCoeffX=np.zeros(int(model*(model+1)/2))
    allMeanGDCoeffY=np.zeros(int(model*(model+1)/2))
    newAllFrameStarList=[]
    ransacDelStarsCounts=0
    reverseFrameList=[] 
    unReverseFrameList=[]
    reverseRotateFrameList=[] 
    unReverseRotateFrameList=[]
    perrxList=[]
    perryList=[]
    
    for headInfo,curFrameStarList in zip(allHeadInfoList,allFrameStarList): 
        meanRA=0
        meanDE=0
        for starInfo in curFrameStarList:
            meanRA+=float(starInfo.RA)/float(headInfo.STARCOUNT)
            meanDE+=float(starInfo.DE)/float(headInfo.STARCOUNT)
            
            
        
        xArray=np.ones(headInfo.STARCOUNT)
        yArray=np.ones(headInfo.STARCOUNT) 
        ksiArray=np.ones(headInfo.STARCOUNT)
        etaArray=np.ones(headInfo.STARCOUNT)
        magArray=np.ones(headInfo.STARCOUNT)
        
        targetX=[]
        targetY=[]
        targetKsi=[]
        targetEta=[]
        for num,starInfo in enumerate(curFrameStarList):
            'ID of gaia source is larger than 10e10'
            if int(starInfo.ID)>10e10:
                xArray[num]=float(starInfo.x)
                yArray[num]=float(starInfo.y)  
                ksiArray[num]=EquatorialToStandard_Xi(float(starInfo.RA) , float(starInfo.DE) ,meanRA,meanDE)*180/np.pi*3600
                etaArray[num]=EquatorialToStandard_Eta(float(starInfo.RA) , float(starInfo.DE) ,meanRA,meanDE)*180/np.pi*3600
                magArray[num]=float(starInfo.Mag)
            else:
                xArray=np.delete(xArray,-1)
                yArray=np.delete(yArray,-1)
                ksiArray=np.delete(ksiArray,-1)
                etaArray=np.delete(etaArray,-1)
                magArray=np.delete(magArray,-1)
                targetX.append(float(starInfo.x))
                targetY.append(float(starInfo.y))
                targetKsi.append(EquatorialToStandard_Xi(float(starInfo.RA) , float(starInfo.DE) ,meanRA,meanDE)*180/np.pi*3600)
                targetEta.append(EquatorialToStandard_Eta(float(starInfo.RA) , float(starInfo.DE) ,meanRA,meanDE)*180/np.pi*3600)
        targetX=np.array(targetX)
        targetY=np.array(targetY)
        targetKsi=np.array(targetKsi)
        targetEta=np.array(targetEta)
        
        if isChangeTransform:
            tempArray=xArray
            xArray=ksiArray
            ksiArray=tempArray
            tempArray=yArray
            yArray=etaArray
            etaArray=tempArray
            
            tempArray=targetX
            targetX=targetKsi
            targetKsi=tempArray
            tempArray=targetY
            targetY=targetEta
            targetEta=tempArray
        
        
        if len(plsqPrecisionCurve.x)!=0:
            weights= 1/fittingFunc(plsqPrecisionCurve.x,magArray) [:len(curFrameStarList)-len(targetX)]
        else:
            weights=np.ones(len(curFrameStarList)-len(targetX))
                     
        if model==1: 
            isReverse=False
            weights=np.append(weights,weights)  
            xCArray=np.append(xArray,-yArray)
            yCArray=np.append(yArray, xArray)
            leftMatrix=genLeftMatrix(xCArray,yCArray,model,weights) 
            leftMatrix=scipy.matrix(leftMatrix)
            allCoeff=np.array(linalg.inv(leftMatrix*leftMatrix.T)*(leftMatrix*scipy.matrix(weights*np.append(ksiArray,etaArray)).T))
            resAll=mylsq.residuals(allCoeff,xCArray,yCArray,np.append(ksiArray,etaArray),model,np.ones(len(xCArray))) 
             
            distantObjA=0
            distantObjB=1
            if len(xArray)>2:
                
                findOutDistantObjs=False
                
                for distantObjA in range(len(xArray)):
                    for distantObjB in range(len(xArray)):    
                        if abs(ksiArray[distantObjA]-ksiArray[distantObjB])>100 and abs(etaArray[distantObjA]-etaArray[distantObjB])>100 :
                            findOutDistantObjs=True
                        if findOutDistantObjs:
                            break
                    if findOutDistantObjs:
                        break
                    
            # Determine whether the image is flipped
            if (xArray[distantObjA]-xArray[distantObjB])/(ksiArray[distantObjA]-ksiArray[distantObjB])/\
            ((yArray[distantObjA]-yArray[distantObjB])/(etaArray[distantObjA]-etaArray[distantObjB]))>0 or resAll.max()>10: 
                isReverse=True
                
                xCArray=np.append( xArray,yArray)
                yCArray=np.append(-yArray, xArray)
                leftMatrix=genLeftMatrix(xCArray,yCArray,model,weights) 
                leftMatrix=scipy.matrix(leftMatrix)
                allCoeff=np.array(linalg.inv(leftMatrix*leftMatrix.T)*(leftMatrix*scipy.matrix(weights*np.append(ksiArray,etaArray)).T))
                resAll=mylsq.residuals(allCoeff,xCArray,yCArray,np.append(ksiArray,etaArray),model,np.ones(len(xCArray))) 
 
                
            tempCoeffs=np.ones(6)
            tempCoeffs[0]=allCoeff[1]
            tempCoeffs[1]=allCoeff[2]
            tempCoeffs[2]=allCoeff[2]
            tempCoeffs[3]=-allCoeff[1]
            tempCoeffs[4]=allCoeff[0]
            tempCoeffs[5]=allCoeff[3]
            if isReverse:
                tempCoeffs[1]=-allCoeff[2]
                tempCoeffs[3]=allCoeff[1]   
            allCoeff=tempCoeffs      
        else:
            if isFitMagTerm:
                magList=list(magArray)
            else:
                magList=[]
            if model*(model+1)/2>headInfo.STARCOUNT-len(targetX):
                deleteFrameList.append(picNum-1)             
                print("The number of stars in frame",picNum, 
                      " is too low to fit the current order model, and will be deleted in the output result")
            leftMatrix=genLeftMatrix(xArray,yArray,model,weights,magList)  
            leftMatrix=scipy.matrix(leftMatrix)
            
            xCoeff=np.array(linalg.inv(leftMatrix*leftMatrix.T)*(leftMatrix*scipy.matrix(weights*ksiArray).T))
            yCoeff=np.array(linalg.inv(leftMatrix*leftMatrix.T)*(leftMatrix*scipy.matrix(weights*etaArray).T))
            resRA=mylsq.residuals(xCoeff,xArray,yArray,ksiArray,model,np.ones(len(curFrameStarList)-len(targetX)),magList)
            resDE=mylsq.residuals(yCoeff,xArray,yArray,etaArray,model,np.ones(len(curFrameStarList)-len(targetX)),magList)
            targetResRA=[]
            targetResDE=[]
            if len(targetX)!=0:
                targetResRA=mylsq.residuals(xCoeff,targetX,targetY,targetKsi,model,np.ones(len(targetX)),magList)
                targetResDE=mylsq.residuals(yCoeff,targetX,targetY,targetEta,model,np.ones(len(targetX)),magList)
             
            allMagList.extend(magArray)
            allWeightedResRAList.extend( weights*resRA )
            allWeightedResDEList.extend( weights*resDE )
            
            
            allME1RAList.append((((weights*resRA)**2).sum()/(len(resRA)-model*(model+1)/2))**0.5)
            allME1DEList.append((((weights*resDE)**2).sum()/(len(resRA)-model*(model+1)/2))**0.5)
            
            
            gridImgPos=np.array([[x,y] for x in range(int(xArray.min()),int(xArray.max()),int((xArray.max()-xArray.min()) /20)) \
                              for y in range(int(yArray.min()),int(yArray.max()),int((yArray.max()-yArray.min()) /20))]).astype(float)
                
            gridXArray=gridImgPos[:,0]
            gridYArray=gridImgPos[:,1]
            order=model-1
            
            if abs(xCoeff[1])>abs(xCoeff[2]):
                isRotate=False
                
                if xCoeff[1]/yCoeff[2]<0:
                    isReverse=False
                    unReverseFrameList.append(headInfo.FILENAME)
                else:
                    isReverse=True
                    reverseFrameList.append(headInfo.FILENAME)
            else:
                isRotate=True
                
                if xCoeff[2]/yCoeff[1]>0:
                    isReverse=False
                    unReverseRotateFrameList.append(headInfo.FILENAME)
                else:
                    isReverse=True
                    reverseRotateFrameList.append(headInfo.FILENAME)
            
            
            
            gridXCoeff,gridYCoeff=alignModel(order,gridXArray,gridYArray,xCoeff,yCoeff,correctOrder,isReverse)
    
    
            if abs(gridXCoeff[1])>abs(gridXCoeff[2]):
                scaleXAppr=gridXCoeff[1]/abs(gridXCoeff[1])
            else:
                scaleXAppr=gridXCoeff[2]/abs(gridXCoeff[2])
            if abs(gridYCoeff[1])>abs(gridYCoeff[2]):
                scaleYAppr=gridYCoeff[1]/abs(gridYCoeff[1])
            else:
                scaleYAppr=gridYCoeff[2]/abs(gridYCoeff[2])  
                
                
            ratioStdToImgae= (gridXCoeff[1]**2+gridXCoeff[2]**2)**0.5
            stdWidth=frameWidth/ratioStdToImgae
            stdHeight=frameHeight/ratioStdToImgae
            boxLenX= ( stdWidth/x_partition)
            boxLenY= (stdHeight /y_partition)
            X = scaleXAppr*np.arange(-stdWidth/2+ boxLenX/2, stdWidth/2-boxLenX/2+1, boxLenX).astype(float)
            Y = scaleYAppr*np.arange(-stdHeight/2+ boxLenY/2, stdHeight/2-boxLenY/2+1, boxLenY) .astype(float) 
            gridValueX=np.zeros((len(X),len(Y)))
            gridValueY=np.zeros((len(X),len(Y)))    
            gridPosX=np.mat(genPolynomialLeftMatrix(X,Y,order)).T  *np.mat( gridXCoeff).T
            gridPosY=np.mat(genPolynomialLeftMatrix(X,Y,order)).T  *np.mat( gridYCoeff).T
            for xNum,xGrid in enumerate(X):
                for yNum,yGrid in enumerate(Y):
                    xIndex=xNum
                    yIndex=yNum
                    
                    if isRotate:
                        xIndex=yNum
                        yIndex=xNum
                        
                        
                    gridValueX[yIndex][xIndex]=np.mat(genPolynomialLeftMatrix(np.array([xGrid]),np.array([yGrid]),order)).T \
                    *np.mat(xCoeff.T-gridXCoeff).T
                    gridValueY[yIndex][xIndex]=np.mat(genPolynomialLeftMatrix(np.array([xGrid]),np.array([yGrid]),order)).T\
                    *np.mat(yCoeff.T-gridYCoeff).T
            
            fig,ax=plt.subplots()  
            ax.quiver(np.array(gridPosX).flatten(), np.array(gridPosY).flatten(), gridValueX,gridValueY,angles='xy', scale_units='xy'
                          , scale=1/quiverScale,        width =0.005 ,headlength =1,headwidth =1  )  
            sX,sY=np.meshgrid(np.array(gridPosX), np.array(gridPosY))
            
            extendFig=0
            if abs(gridValueX.max())>abs(gridValueY.max()):
                extendFig=abs(gridValueX.max())*quiverScale
            else:
                extendFig=abs(gridValueY.max())*quiverScale
            ax.set_xlim(-extendFig,frameWidth+extendFig)
            ax.set_ylim(-extendFig,frameHeight+extendFig)
            
            ax.scatter(sX, sY, color='k', s=10) 
            ax.set_title(headInfo.FILENAME+"     "+str(round(((gridValueX**2+gridValueY**2)**0.5 ).max(),4))+" Pixel"+
                     "  med: "+str(round(stat.median(np.array((gridValueX**2+gridValueY**2)**0.5 ).flatten()),4))+" Pixel"+
                     "  m.e.1: "+str(round((float(allME1RAList[-1])**2+float(allME1DEList[-1])**2)**0.5,4)))
            plt.show()
            xGDCCoeff,yGDCCoeff,xGDCoeff,yGDCoeff=getGDCoeffs(xCoeff,yCoeff,X,Y,gridXCoeff,gridYCoeff,order,GDCModel)
             
            
            poptx=popty=pcovx=pcovy=[]
            xCorrectedArray=np.array(np.mat(genPolynomialLeftMatrix(xArray,yArray,order)).T *np.mat(gridXCoeff).T).flatten()
            yCorrectedArray=np.array(np.mat(genPolynomialLeftMatrix(xArray,yArray,order)).T *np.mat(gridYCoeff).T ).flatten()
            xDistortedArray=ksiArray-xCorrectedArray
            yDistortedArray=etaArray-yCorrectedArray 
            xCorrectedArray=(xCorrectedArray-frameWidth/2)/(frameWidth/2)
            yCorrectedArray=(yCorrectedArray-frameHeight/2)/(frameHeight/2)
            from scipy.optimize import curve_fit
            # only the fitting errors at the second and third orders have been considered now
            if order==3: 
                p0=np.zeros(10)
                poptx, pcovx= curve_fit( polynomialFunc3, xdata=np.array([xCorrectedArray,yCorrectedArray]) 
                                        , ydata=xDistortedArray,p0=p0,maxfev=1000000)
                popty, pcovy= curve_fit( polynomialFunc3, xdata=np.array([xCorrectedArray,yCorrectedArray]) 
                                        , ydata=yDistortedArray,p0=p0,maxfev=1000000)
            elif order==2:  
                p0=np.zeros(6)
                poptx, pcovx= curve_fit( polynomialFunc2, xdata=np.array([xCorrectedArray,yCorrectedArray]) 
                                        , ydata=xDistortedArray,p0=p0,maxfev=1000000)
                popty, pcovy= curve_fit( polynomialFunc2, xdata=np.array([xCorrectedArray,yCorrectedArray]) 
                                        , ydata=yDistortedArray,p0=p0,maxfev=1000000)
            
            perrxList.append(np.sqrt(np.diag(pcovx)))
            perryList.append(np.sqrt(np.diag(pcovy))) 
        
            allFittingErrRAList.append((np.sqrt(np.diag(pcovx))).sum())
            allFittingErrDEList.append((np.sqrt(np.diag(pcovy))).sum())

            
            allGDCoeffX.append(xGDCoeff)
            allGDCoeffY.append(yGDCoeff)
            
            allMeanGDCoeffX+=xGDCoeff
            allMeanGDCoeffY+=yGDCoeff
        
        
        if USE_RANSAC_IN_REDUCTION:
            ransacInliers,delCounts=ransacForTransform(xArray,yArray,ksiArray,etaArray,Residual_Threshold)
            newCurFrameStarList=[]
            for num,val in enumerate(ransacInliers):
                if val :
                    newCurFrameStarList.append(curFrameStarList[num])   
            for targetOrder in range(num+1,len(curFrameStarList)):
                newCurFrameStarList.append(curFrameStarList[targetOrder]) 
            newAllFrameStarList.append(newCurFrameStarList)
            ransacDelStarsCounts+=delCounts
            ksiArray=ksiArray[ransacInliers]
            etaArray=etaArray[ransacInliers]
        else:
            newAllFrameStarList.append(curFrameStarList)
            ransacInliers=np.ones(len(curFrameStarList)-len(targetX),dtype=bool)
        if model==1:
            plsqA=plsq()
            plsqA.x=allCoeff
            plsqA.fun=resAll[np.append(ransacInliers,ransacInliers)]
            fittingResultList.append(plsqA)
        else:
            plsqA=plsq()
            plsqB=plsq()
            plsqA.x=xCoeff
            plsqA.fun=np.append( resRA[ransacInliers],targetResRA)
            plsqB.x=yCoeff
            plsqB.fun=np.append( resDE[ransacInliers],targetResDE)
            fittingResultList.append((plsqA,plsqB)) 
        allKsiList.append(np.append(ksiArray,targetKsi))
        allEtaList.append(np.append(etaArray,targetEta))
        picNum+=1
    if plsqPrecisionCurve!=[]:
        titleStr="weighted residuals RA"
    else:
        titleStr="residuals RA"
    
    
    
    allMeanGDCoeffX/=len(allHeadInfoList) 
    allMeanGDCoeffY/=len(allHeadInfoList)       
    if len(allFittingErrRAList)==0:
        # Do not use weights for GD solution
        for frameNO in range(len(allME1RAList)):
            allMeanGDCoeffX+=allGDCoeffX[frameNO]/len(allME1RAList)
            allMeanGDCoeffY+=allGDCoeffY[frameNO]/len(allME1RAList) 
    else:
       
        # Calculate weights using the sum of fitting errors of normalized parameters
        allMeanGDCoeffX=0
        allMeanGDCoeffY=0
        for frameNO in range(len(allFittingErrRAList)):
            allMeanGDCoeffX+=allGDCoeffX[frameNO]/allFittingErrRAList[frameNO]**2
            allMeanGDCoeffY+=allGDCoeffY[frameNO]/allFittingErrDEList[frameNO]**2
        allMeanGDCoeffX/=(1/np.array(allFittingErrRAList)**2).sum()
        allMeanGDCoeffY/=(1/np.array(allFittingErrDEList)**2).sum()
    
    boxLenX=frameWidth/x_partition
    boxLenY=frameHeight/y_partition
    X = np.arange(boxLenX/2, frameWidth-boxLenX/2+1, boxLenX).astype(float)
    Y = np.arange(boxLenY/2, frameHeight-boxLenY/2+1, boxLenY) .astype(float)
     
    
    xx,yy=np.meshgrid(X,Y) 
    order=model-1   
    for alignTimes in range(3): 
        alignXCoeff,alignYCoeff=alignModel(order,xx.flatten(),yy.flatten(),allMeanGDCoeffX.reshape(len(allMeanGDCoeffX),1)\
                                            ,allMeanGDCoeffY.reshape(len(allMeanGDCoeffY),1),correctOrder,isReverse=True) 
  
        allMeanGDCCoeffX,allMeanGDCCoeffY,allMeanGDCoeffX,allMeanGDCoeffY=\
        getGDCoeffs(allMeanGDCoeffX,allMeanGDCoeffY,xx[0],yy[:,0],alignXCoeff,alignYCoeff,order,GDCModel)
     
    
    global perrx,perry 
    perrx=np.sqrt(np.square(perrxList).sum(axis=0))/len(perrxList)
    perry=np.sqrt(np.square(perryList).sum(axis=0))/len(perryList)
    
    gridValueX=np.zeros((len(X),len(Y)))
    gridValueY=np.zeros((len(X),len(Y)))  
    
    with open("distortionMatrix/"+inputFileName[:-4]+"_GD.txt","w+") as GDFile:
        for yNum,yGrid in enumerate(Y):
            for xNum,xGrid in enumerate(X):
                
    #            gridValueX[yNum][xNum]=xGrid-np.mat(genPolynomialLeftMatrix(np.array([xGrid]),np.array([yGrid]),order)).T \
    #            *np.mat(allMeanGDCCoeffX).T
    #            gridValueY[yNum][xNum]=yGrid-np.mat(genPolynomialLeftMatrix(np.array([xGrid]),np.array([yGrid]),order)).T\
    #            *np.mat(allMeanGDCCoeffY).T 
                'Note that the indexes of matrices and arrays are arranged by row, column, and corresponding coordinates are y, x'
                gridValueX[yNum][xNum]=np.mat(genPolynomialLeftMatrix(np.array([xGrid]),np.array([yGrid]),order)).T \
                *np.mat(allMeanGDCoeffX).T-xGrid
                gridValueY[yNum][xNum]=np.mat(genPolynomialLeftMatrix(np.array([xGrid]),np.array([yGrid]),order)).T\
                *np.mat(allMeanGDCoeffY).T -yGrid
                GDFile.write("{}    {}    {:.4f}    {:.4f}\n".format(xGrid,yGrid,gridValueX[yNum][xNum],gridValueY[yNum][xNum]))
            
    with open("distortionMatrix/"+inputFileName[:-4]+"_reverseGD.txt","w+") as GDFile:
        GDFile.write("{}  {}\n".format(len(X), len(Y)))
        for yNum,yGrid in enumerate(Y):
            for xNum,xGrid in enumerate(X):  
                GDFile.write("{:10.4f}".format(gridValueX[yNum][xNum]))   
            GDFile.write("\n")
                
        GDFile.write("\n")
        for yNum,yGrid in enumerate(Y):
            for xNum,xGrid in enumerate(X):   
                GDFile.write("{:10.4f}".format(gridValueY[yNum][xNum]))   
            GDFile.write("\n")
            
   
    plt.rcParams['figure.figsize'] = (12.0, 12)
    plt.rcParams['figure.dpi'] = 160
    fontScale=1.5
            
    fig,ax=plt.subplots() 
    ax.quiver(X, Y, gridValueX,gridValueY,angles='xy', scale_units='xy'
                  , scale=1/quiverScale,        width =0.005 ,headlength =1,headwidth =1  )  
    extendFig=0
    if abs(gridValueX.max())>abs(gridValueY.max()):
        extendFig=abs(gridValueX.max())*quiverScale
    else:
        extendFig=abs(gridValueY.max())*quiverScale
    ax.set_xlim(-extendFig,frameWidth+extendFig)
    ax.set_ylim(-extendFig,frameHeight+extendFig)
    
    sX,sY=np.meshgrid(X, Y)
    ax.scatter(sX, sY, color='k', s=40*fontScale) 
    titleStr="GD     Max: "+str(round(((gridValueX**2+gridValueY**2)**0.5 ).max(),3))+" Pixel"+\
             "   Median: "+str(round(stat.median(np.array((gridValueX**2+gridValueY**2)**0.5 ).flatten()),3))+" Pixel"

    
    ax.set_title(titleStr,  fontsize=20*fontScale)
    ax.set_xlabel('X (pixel)',size =20*fontScale)
    ax.set_ylabel('Y (pixel)',size =20*fontScale)
    ax.tick_params(labelsize=20*fontScale,length=14*fontScale,width=2.5*fontScale) 
    for sps in ax.spines.values():
        sps.set_linewidth(2.5*fontScale)
    plt.show()
    plotResiduals(allMagList,allWeightedResDEList,titleStr,fontScale )
    allLists=[allHeadInfoList,newAllFrameStarList,fittingResultList,allKsiList,allEtaList]
    deleteFrameByIndex(allLists,deleteFrameList)
    
    
    allHeadInfoList2=[]
    allFrameStarList2=[]
    if(os.path.exists('matchedFiles/'+correctedFileName)):
        inputFilePath="matchedFiles/"
    elif(os.path.exists('modifiedFiles/'+correctedFileName)):
        inputFilePath="modifiedFiles/"
    else:
        print("No matching file specified for geometric distortion correction found: ",correctedFileName) 
    boxLenX=frameWidth/x_partition
    boxLenY=frameHeight/y_partition
    ioFile.read_mathcedStar(allHeadInfoList2,allFrameStarList2,correctedFileName,inputFilePath)
    'Use the obtained distortion correction coefficients to correct image distortion'
    allFrameStarList2=model_Correct_GD(allFrameStarList2,allMeanGDCCoeffX,allMeanGDCCoeffY)
    ioFile.write_processed_matchedStar(allHeadInfoList2,allFrameStarList2,
                                       "GDC_GDmodel_"+str(model)+"_"+correctedFileName,filePath='modifiedFiles/')  
  
    return allLists,allMeanGDCCoeffX,allMeanGDCCoeffY,allMeanGDCoeffX,allMeanGDCoeffY
def plotResiduals(x,y,title="Weighted Residuals",s =0.1,marker=".",fontScale=0.8):
    
    fig,ax=plt.subplots() 
    ax.scatter(x,y,s =s,marker=marker)
    plt.title(title, fontsize=24*fontScale)
    
    ax.set_xlabel('JD',size =26*fontScale)
    ax.set_ylabel('O-C (arcsec)',size =26*fontScale)
    ax.tick_params(labelsize=24*fontScale,length=14*fontScale,width=2.5*fontScale)
    # ax.set_ylim(-0.04,0.04)
    for sps in ax.spines.values():
        sps.set_linewidth(2.5*fontScale)
    plt.show()
def polynomialFunc3(xx,a,b,c,d,e,f,g,h,i,j):
    return a+b*xx[0]+c*xx[1]+d*xx[1]**2+e*xx[0]*xx[1]+f*xx[0]**2+g*xx[1]**3+h*xx[0]*xx[1]**2+i*xx[0]**2*xx[1]+j*xx[0]**3
    
def polynomialFunc2(xx,a,b,c,d,e,f):
    return a+b*xx[0]+c*xx[1]+d*xx[1]**2+e*xx[0]*xx[1]+f*xx[0]**2

 
def getGDCoeffs(xCoeff,yCoeff,X,Y,gridXCoeff,gridYCoeff,order,GDCModel,normalization=False,width=4096,height=4032):
    X,Y=np.meshgrid(X,Y)  
    X=np.array(X.flatten())
    Y=np.array(Y.flatten()) 
    xDistortedArray=np.array(np.mat(genPolynomialLeftMatrix(X,Y,order)).T *np.mat(xCoeff.T).T).flatten()
    yDistortedArray=np.array(np.mat(genPolynomialLeftMatrix(X,Y,order)).T *np.mat(yCoeff.T).T) .flatten()
    xCorrectedArray=np.array(np.mat(genPolynomialLeftMatrix(X,Y,order)).T *np.mat(gridXCoeff).T).flatten()
    yCorrectedArray=np.array(np.mat(genPolynomialLeftMatrix(X,Y,order)).T *np.mat(gridYCoeff).T ).flatten()
    if normalization:
       
        xDistortedArray=(xDistortedArray-width/2)/(width/2)
        yDistortedArray=(yDistortedArray-height/2)/(height/2)
        xCorrectedArray=(xCorrectedArray-width/2)/(width/2)
        yCorrectedArray=(yCorrectedArray-height/2)/(height/2)
        
    GDLeftMatrix=np.mat(genPolynomialLeftMatrix(xDistortedArray,yDistortedArray,GDCModel-1)) 
    xGDCCoeff=np.array(linalg.inv(GDLeftMatrix*GDLeftMatrix.T)*(GDLeftMatrix*scipy.matrix(xCorrectedArray).T)).flatten()
    yGDCCoeff =np.array(linalg.inv(GDLeftMatrix*GDLeftMatrix.T)*(GDLeftMatrix*scipy.matrix(yCorrectedArray).T)).flatten()
    
    from scipy.optimize import curve_fit
    p0=np.zeros(10)
    poptx, pcovx= curve_fit( polynomialFunc3, xdata=np.array([xCorrectedArray,yCorrectedArray]) , ydata=xDistortedArray,p0=p0,maxfev=1000000)
    popty, pcovy= curve_fit( polynomialFunc3, xdata=np.array([xCorrectedArray,yCorrectedArray]) , ydata=yDistortedArray,p0=p0,maxfev=1000000)
    
    global perrx,perry
    
    perrx = np.sqrt(np.diag(pcovx))
    perry = np.sqrt(np.diag(pcovy))
    
    GDLeftMatrix=np.mat(genPolynomialLeftMatrix(xCorrectedArray,yCorrectedArray,order)) 
    xGDCoeff=np.array(linalg.inv(GDLeftMatrix*GDLeftMatrix.T)*(GDLeftMatrix*scipy.matrix(xDistortedArray).T)).flatten()
    yGDCoeff=np.array(linalg.inv(GDLeftMatrix*GDLeftMatrix.T)*(GDLeftMatrix*scipy.matrix(yDistortedArray).T)).flatten()
    return xGDCCoeff,yGDCCoeff,xGDCoeff,yGDCoeff 

'''correctOrder=0 represents the four constant model, and the return value is the correctOrder order transformation parameter 
from the standard coordinate (C) lattice point to the image coordinate (affected by the image GD).'''
def alignModel(order,gridXArray,gridYArray,xCoeff,yCoeff,correctOrder ,isReverse):
    from decimal import localcontext,Decimal 
    if correctOrder>0:
        gridKsi       =np.mat(genPolynomialLeftMatrix(gridXArray,gridYArray,order),dtype=Decimal).T*np.mat(xCoeff,dtype=Decimal)
        gridEta       =np.mat(genPolynomialLeftMatrix(gridXArray,gridYArray,order),dtype=Decimal).T*np.mat(yCoeff,dtype=Decimal)
        gridLeftMatrix=np.mat(genPolynomialLeftMatrix(gridXArray,gridYArray, correctOrder),dtype=Decimal)
        
        with localcontext() as ctx:
            ctx.prec = 32  
            f_matrix = np.array((gridLeftMatrix*gridLeftMatrix.T).tolist(), dtype=float)
            f_inv = np.linalg.inv(f_matrix)
        
        
        gridXCoeff=np.append(np.array(f_inv*(gridLeftMatrix*gridKsi))\
                            ,np.zeros(((len(xCoeff)-int((correctOrder+1)*(correctOrder+2)/2))),dtype=Decimal)).astype(np.float64)
        gridYCoeff=np.append(np.array(f_inv*(gridLeftMatrix*gridEta))\
                            ,np.zeros(((len(yCoeff)-int((correctOrder+1)*(correctOrder+2)/2))),dtype=Decimal)).astype(np.float64)
    
    else: 
        gridCombineX=np.append(gridXArray,-gridYArray)
        gridCombineY=np.append(gridYArray,gridXArray)
        if isReverse:
            gridCombineX=np.append(gridXArray,gridYArray)
            gridCombineY=np.append(-gridYArray,gridXArray)
            
        gridLeftMatrix=np.mat([np.append(np.ones(len(gridXArray)),np.zeros(len(gridYArray))),gridCombineX,gridCombineY,
         np.append(np.zeros(len(gridXArray)),np.ones(len(gridYArray)))])
        gridStd=np.append(np.asarray(np.mat(genPolynomialLeftMatrix(gridXArray,gridYArray,order)).T*np.mat(xCoeff))
                        ,np.asarray(np.mat(genPolynomialLeftMatrix(gridXArray,gridYArray,order)).T*np.mat(yCoeff)))
        gridCoeff= np.array(linalg.inv(gridLeftMatrix*gridLeftMatrix.T)*(gridLeftMatrix*np.mat(gridStd).T)) 
        # Complete the number of fitting parameters to match the number of parameters in the plate constant model
        gridXCoeff=np.append(np.array([gridCoeff[0],xCoeff[1]/abs(xCoeff[1])*abs(gridCoeff[1]),xCoeff[2]/abs(xCoeff[2])*abs(gridCoeff[2])])\
                             ,np.zeros(((len(xCoeff)-3))))
        gridYCoeff=np.append(np.array([gridCoeff[3],yCoeff[1]/abs(yCoeff[1])*abs(gridCoeff[2]),yCoeff[2]/abs(yCoeff[2])*abs(gridCoeff[1])])\
                             ,np.zeros(((len(yCoeff)-3))))
        if isReverse:
            scaleSign=-1
        else:
            scaleSign= 1
            
        gridXCoeff=np.append(np.array([gridCoeff[0],gridCoeff[1],scaleSign*gridCoeff[2]]) ,np.zeros(((len(xCoeff)-3))))
        gridYCoeff=np.append(np.array([gridCoeff[3],gridCoeff[2],scaleSign*-gridCoeff[1]]) ,np.zeros(((len(yCoeff)-3))))
    return gridXCoeff,gridYCoeff


'Statistical accuracy, draw accuracy curves, and eliminate some stars with large measurement standard deviations'
def calcPrecisionCurve(bWriteStatResult,fittedOutputFileName='',excludeStarDict=dict(),brightStarMag=14,PAPER_SD_FIG_TITLE=""):
  
    outputList,meanRAOCDict,meanDEOCDict=MuApp.processFile(bWriteStatResult,fittedOutputFileName,leastNum=3)
    plsqPrecisionCurve,excludeStarDict=pltS.precisionCurve(outputList,excludeStarDict,"",PAPER_SD_FIG_TITLE,brightStarMag )
    return plsqPrecisionCurve,meanRAOCDict,meanDEOCDict,excludeStarDict,outputList


def ransacForTransform(xArray,yArray,ksiArray,etaArray,residual_threshold=0.2): 
    from skimage.measure import ransac 
    from GDSolver.mylib.ransacModels import MyModel
    data=np.c_[xArray.T,yArray.T]
    data=np.c_[data,ksiArray.T]
    data=np.c_[data,etaArray.T]  
    model_robust, inliers = ransac(data, MyModel , min_samples=15,
                               residual_threshold=residual_threshold, max_trials=1000)   
    return inliers,len(inliers)-inliers.sum() 
def updateProcessedFile(allLists,isMinusMeanOC,inputFileName,excludeStarDict,meanRAOCDict,meanDEOCDict):
    ioFile.write_processed_matchedStar(allLists[0],allLists[1],"intermediate_"+inputFileName,
                                       "tempFiles/",excludeStarDict, isMinusMeanOC,meanRAOCDict,meanDEOCDict) 
    inputFileNamePrefix="intermediate_"
    return inputFileNamePrefix

def copyFittedResultToStat(allLists,fittedFileName):
    import shutil 
    fittedOutputFileName=str(len(allLists[0]))+'_'+fittedFileName
    oldname='outputFiles/'+fittedOutputFileName
    newname='MutualApproximation/outputFiles/'+fittedOutputFileName
    shutil.copyfile(oldname,newname)
    return fittedOutputFileName
def updateMeanOCSum(meanRAOCDict,meanDEOCDict,meanRAOCDictSum,meanDEOCDictSum):
    if meanRAOCDictSum==dict():
        meanRAOCDictSum=meanRAOCDict
        meanDEOCDictSum=meanDEOCDict
    else:
        for key in meanRAOCDictSum:
            if meanRAOCDict.get(key)!=None and meanRAOCDictSum.get(key)!=None:
                meanRAOCDictSum[key][0]+=meanRAOCDict.get(key)[0]
                meanRAOCDictSum[key][1]=meanRAOCDict.get(key)[1]
                meanRAOCDictSum[key][2] =meanRAOCDict.get(key)[2]
                meanRAOCDictSum[key][3] =meanRAOCDict.get(key)[3]
            if meanDEOCDict.get(key)!=None and meanDEOCDictSum.get(key)!=None:
                meanDEOCDictSum[key][0]+=meanDEOCDict.get(key)[0]
                meanDEOCDictSum[key][1] =meanDEOCDict.get(key)[1] 
    return meanRAOCDictSum,meanDEOCDictSum
def reduceData(model,inputFileName,ITERATE_TIMES,bWriteFittedResult=True,bRefineResult=False ,bWriteStatResult=True):

    excludeStarDict=dict()
    plsqPrecisionCurve=plsq()
    inputFileNamePrefix=""
    meanRAOCDictSum=dict()
    meanDEOCDictSum=dict()
    for iterateFitTimes in range(ITERATE_TIMES):
    
        allLists,allMeanGDCCoeffX,allMeanGDCCoeffY,allMeanGDCoeffX,allMeanGDCoeffY\
        =solvePlateConstants(inputFileName, model , plsqPrecisionCurve, inputFileNamePrefix,isFitMagTerm=False
                                     ,isChangeTransform=True,minMag=MIN_PROCCESSING_MAG,maxMag=MAX_PROCCESSING_MAG)
        
        plt.rcParams['figure.figsize'] = (8.0, 4)
        plt.rcParams['figure.dpi'] = 80
        fontScale=0.8
        
        if bRefineResult==True:    
            refineFrameByCosts(allLists)
        
        ioFile.write_fitted_result( model,allLists,fittedFileName,bOutputHead=True)
        
        fittedOutputFileName=str(len(allLists[0]))+'_'+fittedFileName
        
        plsqPrecisionCurve,meanRAOCDict,meanDEOCDict,excludeStarDict,outputList=\
        calcPrecisionCurve(bWriteStatResult,fittedOutputFileName,excludeStarDict, 16,PAPER_FIG_TITLE)
        
        
        plt.rcParams['figure.figsize'] = (6.0, 6)
        plt.rcParams['figure.dpi'] = 80
        fontScale=0.8
        
        meanRAOCDictSum,meanDEOCDictSum=updateMeanOCSum(meanRAOCDict,meanDEOCDict,meanRAOCDictSum,meanDEOCDictSum)
        if iterateFitTimes !=ITERATE_TIMES-1: 
            isMinusMeanOC=False
            if iterateFitTimes>0: 
                isMinusMeanOC=True
            if NO_STAR_MODIFY:
                isMinusMeanOC=False 
            if NO_STAR_REJECT:
                excludeStarDict=dict() 
            inputFileNamePrefix=updateProcessedFile(allLists,isMinusMeanOC,inputFileName,excludeStarDict,meanRAOCDict,meanDEOCDict)
    with open("statisticRes/stat_meanOCSum_"+inputFileName,"w") as output:
        for dictKey in meanRAOCDictSum:
            if meanRAOCDictSum[dictKey][1]<1000:
                'ID mag meanRA meanDE stdRA stdDE num'
                outputStr="{}   {}  {:>10.4f}   {:>10.4f}  {:>10.4f}   {:>10.4f}    {}    {:>10.4f}".\
                format(dictKey,meanRAOCDictSum[dictKey][2],meanRAOCDictSum[dictKey][0],meanDEOCDictSum[dictKey][0],
                             meanRAOCDictSum[dictKey][1],meanDEOCDictSum[dictKey][1],meanRAOCDictSum[dictKey][3],\
                             ( meanRAOCDictSum[dictKey][1]**2+meanDEOCDictSum[dictKey][1]**2)**0.5)
                output.write(outputStr+"\n")

    return meanRAOCDictSum,meanDEOCDictSum,excludeStarDict,allMeanGDCCoeffX,allMeanGDCCoeffY,allMeanGDCoeffX,allMeanGDCoeffY
 
if __name__=='__main__':

    meanRAOCDictSum,meanDEOCDictSum,excludeStarDict,allMeanGDCCoeffX,allMeanGDCCoeffY,allMeanGDCoeffX,allMeanGDCoeffY\
    =reduceData(model,inputFileName,ITERATE_TIMES,  bWriteFittedResult=True,bRefineResult=False ,bWriteStatResult=True)
    print("allMeanCoeffX=",allMeanGDCoeffX)
    print("allMeanCoeffY=",allMeanGDCoeffY)
    print("allMeanGDCCoeffX=",allMeanGDCCoeffX)
    print("allMeanGDCCoeffY=",allMeanGDCCoeffY)
    print("perrx=",perrx)
    print("perry=",perry)
 