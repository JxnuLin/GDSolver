# -*- coding: utf-8 -*-
import  sys ,os
sys.path.append("..") 
import GDSolver.mylib.ioFile_ReadO_C as ioFile
import GDSolver.mylib.utils as cons
import statistics as stat
import numpy as np 

OUTPUT_PATH=os.path.dirname(__file__) +'/mutualAppximationResults/' 
inputFileName='38_38_Ref_very_17_c1_pickSameStar.txt'
 
MAG_PARTITION_RANGE=[14,15,17] 
MAX_MAGINDEX=30


 
         
def genMeanOCDict(allHeadInfoList,allFrameStarList,isStatRMS=False,leastNum=3):
    meanRAOCDict=dict()
    meanDEOCDict=dict() 
    for headInfo,curFrameStarList in zip(allHeadInfoList,allFrameStarList): 
        for num,starInfo in enumerate(curFrameStarList):
            if meanRAOCDict.get(starInfo.ID)==None:
                meanRAOCDict[starInfo.ID]=[[float(starInfo.dRA),float(starInfo.Mag),float(starInfo.pmRA),float(starInfo.pmDE)]]
                meanDEOCDict[starInfo.ID]=[float(starInfo.dDE)]
            else:
                meanRAOCDict[starInfo.ID].append([float(starInfo.dRA),float(starInfo.Mag),float(starInfo.pmRA),float(starInfo.pmDE)])
                meanDEOCDict[starInfo.ID].append(float(starInfo.dDE))
    for dictKey in meanRAOCDict:
        if len(meanRAOCDict[dictKey])>=leastNum and len(meanRAOCDict[dictKey])>1: 
            if stat.mean(np.array(meanRAOCDict[dictKey])[:,0]) >100 or stat.stdev(np.array(meanRAOCDict[dictKey])[:,0])>100\
                or stat.mean(meanDEOCDict[dictKey])>100 or stat.stdev(meanDEOCDict[dictKey])>100 \
                or stat.mean(np.array(meanRAOCDict[dictKey])[:,0])>100 or\
                        ((np.array(meanRAOCDict[dictKey])[:,0]**2).sum()/len(np.array(meanRAOCDict[dictKey])[:,0]))**0.5>100\
                or stat.mean(meanDEOCDict[dictKey]) >100 or\
                 ((np.array(meanDEOCDict[dictKey])**2).sum()/len(meanDEOCDict[dictKey]))**0.5>100:
                     print(dictKey,'Statistics greater than 100 will not be returned',end="")
            inliersIndexs,deleteCount=cons.DeleteOutliers(np.array(meanRAOCDict[dictKey])[:,0],3)
            if cons.isTarget(dictKey) and deleteCount>0:
                print("The measurement of RA at position ",dictKey," of the target with ",deleteCount," frames may be incorrect. Consider removing it")
                print(inliersIndexs)
            if not cons.isTarget(dictKey):
                
                meanRAOCDict[dictKey]=list(np.array(meanRAOCDict[dictKey])[inliersIndexs])
                
            if not isStatRMS:
                meanRAOCDict[dictKey]=[stat.mean(np.array(meanRAOCDict[dictKey])[:,0]),
                        stat.stdev(np.array(meanRAOCDict[dictKey])[:,0]),meanRAOCDict[dictKey][0][1],len(meanRAOCDict[dictKey])
                        ,meanRAOCDict[dictKey][0][2],meanRAOCDict[dictKey][0][3]] 
            else:
                meanRAOCDict[dictKey]=[stat.mean(np.array(meanRAOCDict[dictKey])[:,0]),
                        ((np.array(meanRAOCDict[dictKey])[:,0]**2).sum()/len(np.array(meanRAOCDict[dictKey])[:,0]))**0.5
                         ,meanRAOCDict[dictKey][0][1],len(meanRAOCDict[dictKey],meanRAOCDict[dictKey][0][2],meanRAOCDict[dictKey][0][3])]
            inliersIndexs,deleteCount=cons.DeleteOutliers(np.array(meanDEOCDict[dictKey]),3)
            if cons.isTarget(dictKey) and deleteCount>0:
                print("The measurement of DE at position ",dictKey," of the target with ",deleteCount," frames may be incorrect. Consider removing it")
                print(inliersIndexs)
            
            if not  cons.isTarget(dictKey):           
                meanDEOCDict[dictKey]=list(np.array(meanDEOCDict[dictKey])[inliersIndexs])
            
            if not isStatRMS:
                meanDEOCDict[dictKey]=[stat.mean(meanDEOCDict[dictKey]),stat.stdev(meanDEOCDict[dictKey])] 
            else:
                meanDEOCDict[dictKey]=[stat.mean(meanDEOCDict[dictKey])
                , ((np.array(meanDEOCDict[dictKey])**2).sum()/len(meanDEOCDict[dictKey]))**0.5]
            
        else:
            'Mean, standard deviation, magnitude, number of statistics, proper motion in RA, proper motion in DE'
            meanRAOCDict[dictKey]=[100,100,meanRAOCDict[dictKey][0][1],len(meanRAOCDict[dictKey]),1000,1000]
            meanDEOCDict[dictKey]=[100,100]
            
    return meanRAOCDict,meanDEOCDict   
def processFile(bWriteResult=True,inputFileName='' ,isStatRMS=False,leastNum=3):
    import os
    allFrameStarList=[] 
    allHeadInfoList=[]
    outputList=[] 
    
    ioFile.read_mathcedStar(allHeadInfoList,allFrameStarList,inputFileName)
    meanRAOCDict,meanDEOCDict=genMeanOCDict(allHeadInfoList,allFrameStarList,isStatRMS,leastNum)
    if bWriteResult:
        with open(os.path.dirname(os.path.dirname(__file__))+"/"+"statisticRes/stat_"+inputFileName,"w") as output:
            outputStr='ID   mag     meanRA  meanDE  stdRA   stdDE   num     std_all     pmRA    pmDE'
            output.write(outputStr+"\n")
            for dictKey in meanRAOCDict:
                'ID mag meanRA meanDE stdRA stdDE num'
                outputStr="{}   {}  {:>10.4f}   {:>10.4f}  {:>10.4f}   {:>10.4f}    {}    {:>10.4f}  {:>10.3f}  {:>10.3f}".format(dictKey
                               ,meanRAOCDict[dictKey][2],meanRAOCDict[dictKey][0],meanDEOCDict[dictKey][0],
                                 meanRAOCDict[dictKey][1],meanDEOCDict[dictKey][1],meanRAOCDict[dictKey][3]
                                 ,( meanRAOCDict[dictKey][1]**2+meanDEOCDict[dictKey][1]**2)**0.5,meanRAOCDict[dictKey][4],meanRAOCDict[dictKey][5])
                 
                if meanRAOCDict[dictKey][1]<100:
                    output.write(outputStr+"\n")
                    outputList.append(outputStr)
    else: 
        for dictKey in meanRAOCDict:
            if meanRAOCDict[dictKey][1]<100: 
                'ID mag meanRA meanDE stdRA stdDE num'
                outputStr="{}   {}  {:>10.4f}   {:>10.4f}  {:>10.4f}   {:>10.4f}    {}    {:>10.4f}  {:>10.3f}  {:>10.3f}".format(dictKey
                       ,meanRAOCDict[dictKey][2],meanRAOCDict[dictKey][0],meanDEOCDict[dictKey][0],
                         meanRAOCDict[dictKey][1],meanDEOCDict[dictKey][1],meanRAOCDict[dictKey][3]
                         ,( meanRAOCDict[dictKey][1]**2+meanDEOCDict[dictKey][1]**2)**0.5,meanRAOCDict[dictKey][4],meanRAOCDict[dictKey][5])
                outputList.append(outputStr)
    print("mean diff")
    return outputList,meanRAOCDict,meanDEOCDict
if __name__=='__main__':
    processFile(True,inputFileName)
    
sys.path.pop()