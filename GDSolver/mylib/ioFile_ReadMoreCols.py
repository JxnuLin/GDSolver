# -*- coding: utf-8 -*-
"""
Created on Wed Nov 29 15:28:46 2017

@author: Administrator
"""

import re 
import os
import sys 
from collections import namedtuple
import statistics as stat
import GDSolver.mylib.utils as cons
sys.path.append("..")
HeadInfo=namedtuple('HeadInfo',['FILENAME','STARCOUNT','YEAR','MONTH','DAY','HOUR','MINUTE','SECOND','EXPOSURE'])
StarInfo=namedtuple('StarInfo',['ID','x','y','RA','DE','Mag','pmRA','pmDE','plx','rv']) 
def renameFile(fileName,addStr): 
    fileName=fileName[:fileName.rfind('.')]
    if os.access(os.path.dirname(os.path.dirname(__file__))+'/outputFiles/'+addStr+'_' +fileName+'.txt',os.F_OK):
        os.remove(os.path.dirname(os.path.dirname(__file__))+'/outputFiles/' +addStr+'_' +fileName+'.txt')
    os.rename(os.path.dirname(os.path.dirname(__file__))+'/outputFiles/'+fileName+'.txt'
              ,os.path.dirname(os.path.dirname(__file__))+'/outputFiles/'+addStr+'_' +fileName+'.txt')
    return addStr+'_'+fileName+'.txt'

def write_processed_matchedStar(allHeadInfoList,allFrameStarList,fileName,filePath=os.path.dirname(os.path.dirname(__file__))+"/outputFiles/"
                                ,excludeStarDict=dict(),isMinusMeanOC=False,meanRAOCDict=dict(),meanDEOCDict=dict(),allMoonsList=[]): 
    frameCount=0 
    with open(filePath+fileName, 'wt') as file:
        if allMoonsList==[]:
            zipList=zip(allHeadInfoList,allFrameStarList,[[] for i in range(len(allHeadInfoList))])
        
        elif allMoonsList!=[] and len(allMoonsList)!=len(allHeadInfoList):
            print("Satellite list error, a frame is missing satellites") 
        else:
            zipList=zip(allHeadInfoList,allFrameStarList,allMoonsList)
        for headInfo,curFrameStarList,curFrameMoonList in zipList: 
            frameCount+=1
            curFrameNewStarCount=0 
            for starInfo in curFrameStarList:
                if excludeStarDict.get(starInfo.ID)==None:
                    curFrameNewStarCount+=1
            file.write('{FILENAME}    {newSTARCOUNT}    {YEAR}    {MONTH:<10s}{DAY:<10s}{HOUR:<10s}{MINUTE:<10s}{SECOND:<10s}{EXPOSURE:<10s}\n'.format(**headInfo._asdict(),newSTARCOUNT=curFrameNewStarCount))
 
            for starInfo in curFrameStarList:
                if excludeStarDict.get(starInfo.ID)==None:
                    if isMinusMeanOC:
                        'Set the average O-C of stars with less than 5 statistics in the statistical function to 100'
                        if meanRAOCDict[starInfo.ID][0]<100 and int(starInfo.ID)>1e10:
                            fRA=float(starInfo.RA)+meanRAOCDict[starInfo.ID][0]/3600
                            fDE=float(starInfo.DE)+meanDEOCDict[starInfo.ID][0]/3600
                        else:
                            fRA=float(starInfo.RA)
                            fDE=float(starInfo.DE)
                    else:
                        fRA=float(starInfo.RA)
                        fDE=float(starInfo.DE)
                    file.write('{ID:>20} {x:>12.9s} {y:>12.9s}    {fRA:>22.11f}    {fDE:>+22.11f}    {Mag}     {fPmRA:>6.3f}     {fPmDE:>6.3f}     {plx:>10}     {rv:>10}\n'
                               .format(**starInfo._asdict(),fRA=fRA,fDE=fDE,fPmRA=float(starInfo.pmRA), fPmDE=float(starInfo.pmDE)))
            for moons in curFrameMoonList:
                file.write(moons)
           
            
#    return renameFile(fileName,str(frameCount))

def write_fitted_result(isFourParam ,allLists,fileName,bOutputHead=True,bSeparateStars=True):
    '''Output the fitting results, where the parameters fittingResultList, isFourParam, allHeadInfoList, 
    and allFrameStarList contain the model fitting results. bOutputHead represents whether to output header file information for each frame.
    If not, sort the output stars by their IDs, and use bSeparateStars to indicate whether to use split lines to segment different stars'''
    [allHeadInfoList,allFrameStarList,fittingResultList,allKsiList,allEtaList]=allLists
    frameCount=0  
    if bOutputHead:
        with open(os.path.dirname(os.path.dirname(__file__))+'/outputFiles/'+fileName, 'wt') as file:
            for headInfo,curFrameStarList,curFrameFittedResult,ksiArray,etaArray in zip(allHeadInfoList,allFrameStarList,fittingResultList,allKsiList,allEtaList):
                frameCount+=1
#            
                if isFourParam==1:  
                    file.write('{FILENAME}    {STARCOUNT}    {YEAR}  {MONTH:<5s} {DAY:<5s} {HOUR:<5s} {MINUTE:<5s} {SECOND:<6s} {EXPOSURE:<6s} {JD:>10.8f} {a:>10.8f} {b:>10.8f} {c:>10.8f} {d:>10.8f}\n'
                               .format(**headInfo._asdict(),a=float(curFrameFittedResult.x[0]),b=float(curFrameFittedResult.x[1]),
                               c=float(curFrameFittedResult.x[2]),d=float(curFrameFittedResult.x[3]),
                               JD=cons.julian_date(headInfo.YEAR,headInfo.MONTH,headInfo.DAY,float(headInfo.HOUR)+float(headInfo.MINUTE)/60+(float(headInfo.SECOND)+float(headInfo.EXPOSURE)/2)/3600)))
                    for starInfo,residualsRA,residualsDE,ksi,eta in \
                    zip(curFrameStarList,curFrameFittedResult.fun[:len(curFrameFittedResult.fun)//2],curFrameFittedResult.fun[len(curFrameFittedResult.fun)//2:len(curFrameFittedResult.fun)],ksiArray,etaArray):
 
                        file.write('{ID:>20}{x:>13.8s}{y:>13.8s}    {RA:0<14s}    {DE:0<14s} {fMag:>8.4f} {fPmRA:>10.3f} {fPmDE:>10.3f} {ksi:>15.4f} {eta:>15.4f} {residualsRA:>11.4f} {residualsDE:>11.4f}\n'
                                   .format(**starInfo._asdict(),fPmRA=float(starInfo.pmRA),fMag=float(starInfo.Mag), fPmDE=float(starInfo.pmDE),ksi=ksi,eta=eta,
                                   residualsRA=residualsRA,residualsDE=float(residualsDE)))
                else:
                    file.write('{FILENAME}    {STARCOUNT}    {YEAR}  {MONTH:<5s} {DAY:<5s} {HOUR:<5s} {MINUTE:<5s} {SECOND:<6s} {EXPOSURE:<6s} {JD:>10.8f} {a:>10.8f} {b:>10.8f} {c:>10.8f} {d:>10.8f}\n'
                               .format(**headInfo._asdict(),a=float(curFrameFittedResult[0].x[1]),b=float(curFrameFittedResult[0].x[2]),
                               c=float(curFrameFittedResult[1].x[1]),d=float(curFrameFittedResult[1].x[2]),
                               JD=cons.julian_date(headInfo.YEAR,headInfo.MONTH,headInfo.DAY,float(headInfo.HOUR)+float(headInfo.MINUTE)/60+(float(headInfo.SECOND)+float(headInfo.EXPOSURE)/2)/3600)))
                    for starInfo,residualsRA,residualsDE,ksi,eta in zip(curFrameStarList,curFrameFittedResult[0].fun,curFrameFittedResult[1].fun,ksiArray,etaArray):
                        file.write('{ID:>20}{x:>13.8s}{y:>13.8s}    {RA:0<14s}    {DE:0<14s} {fMag:>8.4f} {fPmRA:>10.3f} {fPmDE:>10.3f} {ksi:>15.4f} {eta:>15.4f} {residualsRA:>11.4f} {residualsDE:>11.4f}\n'
                                   .format(**starInfo._asdict(),fPmRA=float(starInfo.pmRA), fMag=float(starInfo.Mag),fPmDE=float(starInfo.pmDE),ksi=ksi,eta=eta,
                                   residualsRA=residualsRA,residualsDE=float(residualsDE)))
    else: 
        starAllResidualsList=[]
        oriFileName=fileName[:fileName.rfind('.')]
        for headInfo,curFrameStarList,curFrameFittedResult in zip(allHeadInfoList,allFrameStarList,fittingResultList):
            frameCount+=1
            
            if isFourParam==1: 
                for starInfo,residualsRA,residualsDE in zip(curFrameStarList,curFrameFittedResult.fun[:len(curFrameFittedResult.fun)//2],curFrameFittedResult.fun[len(curFrameFittedResult.fun)//2:len(curFrameFittedResult.fun)]):
                    starAllResidualsList.append([starInfo.ID,format(residualsRA,'>6.3f'),format(residualsDE, '>6.3f'),format(starInfo.x,'<10.8s'),format(starInfo.y,'<10.8s'),format(starInfo.RA,'<13.13s'),format(starInfo.DE,'<13.13s'),
                                                 format(starInfo.Mag,'<10.8s'),format(starInfo.pmRA,'<10.8s'),format(starInfo.pmDE,'<10.8s')])
            else: 
                for starInfo,residualsRA,residualsDE in zip(curFrameStarList,curFrameFittedResult[0].fun,curFrameFittedResult[1].fun):
                    starAllResidualsList.append([starInfo.ID,format(residualsRA,'>6.3f'),format(residualsDE, '>6.3f'),format(starInfo.x,'<10.8s'),format(starInfo.y,'<10.8s'),format(starInfo.RA,'<13.13s'),format(starInfo.DE,'<13.13s'),
                                                 format(starInfo.Mag,'<10.8s'),format(starInfo.pmRA,'<10.8s'),format(starInfo.pmDE,'<10.8s')])
        lastStarID=0
        starAllResidualsList.sort(key=lambda starsResidualList:starsResidualList[0]) 
        fileName=oriFileName+'_StarsAllResiduals.txt'
        with open(os.path.dirname(os.path.dirname(__file__))+"/outputFiles/"+fileName, 'wt') as file:
            for starResidual in starAllResidualsList:
                if(lastStarID!=starResidual[0]) and bSeparateStars:
                    file.write(format("","->200s")+"\n")
                file.write(str(starResidual).replace("'",'  ').replace(',','  ').replace('[','').replace(']','')+'\n')
                lastStarID=starResidual[0]
                
        starStatResidualsList=[] 
        fileName=oriFileName+'_StarsStatResiduals.txt'
        curStarRecidualsRA=[]
        curStarRecidualsDE=[]
        lastStarResidual=starAllResidualsList[0] 
        for starResidual in starAllResidualsList:
#            print(lastStarID)
            if lastStarResidual[0]==starResidual[0]:
                curStarRecidualsRA.append(float(starResidual[1]))
                curStarRecidualsDE.append(float(starResidual[2]))
            else: 
                curStarStats=len(curStarRecidualsRA)
                if curStarStats>1:
                    starStatResidualsList.append([lastStarResidual[0],stat.stdev(curStarRecidualsRA),stat.stdev(curStarRecidualsDE),curStarStats,
                                                  stat.mean(curStarRecidualsRA),stat.mean(curStarRecidualsDE),float(lastStarResidual[3]),float(lastStarResidual[4]),
                                                  lastStarResidual[5],lastStarResidual[6] ,lastStarResidual[7], lastStarResidual[8], lastStarResidual[9]])
                else:
                    starStatResidualsList.append([lastStarResidual[0],0,0,curStarStats,0,0,float(lastStarResidual[3]),float(lastStarResidual[4]),
                                                  lastStarResidual[5],lastStarResidual[6] ,lastStarResidual[7], lastStarResidual[8], lastStarResidual[9]])
                curStarRecidualsRA=[]
                curStarRecidualsDE=[]
                curStarRecidualsRA.append(float(starResidual[1]))
                curStarRecidualsDE.append(float(starResidual[2]))
            lastStarResidual=starResidual   
        curStarStats=len(curStarRecidualsRA)
        if curStarStats>1:
            starStatResidualsList.append([lastStarResidual[0],stat.stdev(curStarRecidualsRA),stat.stdev(curStarRecidualsDE),curStarStats,
                                          stat.mean(curStarRecidualsRA),stat.mean(curStarRecidualsDE),float(lastStarResidual[3]),float(lastStarResidual[4]),
                                          lastStarResidual[5],lastStarResidual[6] ,lastStarResidual[7], lastStarResidual[8], lastStarResidual[9]])
        else:
            starStatResidualsList.append([lastStarResidual[0],0,0,curStarStats,0,0,float(lastStarResidual[3]),float(lastStarResidual[4]),
                                          lastStarResidual[5],lastStarResidual[6] ,lastStarResidual[7], lastStarResidual[8], lastStarResidual[9]])
        with open('outputFiles/'+fileName, 'wt') as file:
            for starStatResidual in starStatResidualsList:
                file.write("{} {:>6.3f} {:>6.3f} {:>4} {:>8.3f} {:>8.3f} {:>12.4f} {:>12.4f}  {}  {}  {} {} {}".format(*starStatResidual)+'\n')
            
    return renameFile(fileName,str(frameCount))

''' Read the matching result file and save it in allHeadInfoList, allFrameStarList,
 where each element in allHeadInfoList is a header file information.
 Each element in allFrameStarList is a list of all stars matched in a frame, 
 where each element represents a matching star information (refer to the format of the matching result text) '''
def read_mathcedStar(allHeadInfoList,allFrameStarList,fileName,filePath=os.path.dirname(os.path.dirname(__file__))+'/outputFiles/'
                     ,allMoonsList=[],excludeReadStarDict=dict(),minMag=-9999,maxMag=9999):
    
    curFrameStarList=[] 
    moonsList=[]
    curFrameStarCount=0
    curHeadInfo= HeadInfo('FILENAME', 'STARCOUNT', 'YEAR', 'MONTH', 'DAY', 'HOUR', 'MINUTE', 'SECOND',   'EXPOSURE') 
    with open(filePath+fileName, 'rt') as file:
        for line in file:
            line=line.strip() 
            fields= re.split(r'\s+',line) 
            if len(fields) ==9: 
                if curFrameStarCount!=0: 
                    newHeadInfo =HeadInfo(*curHeadInfo[:1],curFrameStarCount,*curHeadInfo[2:9])
                    if curFrameStarCount!=int(curHeadInfo.STARCOUNT):
                        print("The number of rows and stars in image"+curHeadInfo.FILENAME+" does not match! The number of stars has been corrected to the number of rows in the read in result!" )
                    allHeadInfoList.append(newHeadInfo)
                    allFrameStarList.append(curFrameStarList)
                    allMoonsList.append(moonsList)
                    curFrameStarList=[]
                    moonsList=[]
                curFrameStarCount=0
                curHeadInfo= HeadInfo(*fields)
                #print(curHeadInfo)
            elif len(fields)==10:
                'Exclude stars specified in excludeReadStarDict and do not read them'
                if excludeReadStarDict.get(fields[0])!=None:
                    continue
                if float(fields[5])>=maxMag or float(fields[5])<=minMag:
                    continue
                
                curFrameStarList.append(StarInfo(*fields))
                curFrameStarCount+=1
            elif len(fields)==2:
                moonsList.append(line+"\n")
                
        newHeadInfo =HeadInfo(*curHeadInfo[:1],curFrameStarCount,*curHeadInfo[2:9])
        allHeadInfoList.append(newHeadInfo)
        allFrameStarList.append(curFrameStarList)
        allMoonsList.append(moonsList)
        if len(curFrameStarList)==0:
            print("No valid lines were read, possibly due to incorrect input file format!")

 