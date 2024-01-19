# -*- coding: utf-8 -*-
from matplotlib import pyplot as plt
import re, sys,os
sys.path.append("..")
from scipy.optimize import least_squares,curve_fit
import math 
import numpy as np 
import GDSolver.mylib.ioFile_ReadMoreCols  as ioFile 
import statistics as stat 
from GDSolver.mylib.utils import isTarget,genExcludeStarDict,plsq 
from matplotlib.ticker import AutoMinorLocator ,MultipleLocator  
plt.rcParams['figure.figsize'] = (8.0, 4)
plt.rcParams['figure.dpi'] = 80
fontScale=0.8
pcRangeFactor=1 
 

'y = A2 + (A1-A2)/(1 + exp((x-x0)/dx))'
def fittingFunc(p0,x):
    return   p0[1] + (p0[0]-p0[1])/(1 + math.e**((x-p0[2])/p0[3])) 
def fittingFuncResidual(p,x,obj):
    return   fittingFunc(p,x)- obj     

def plotMean(ax,statListSorted,lowerLimit=0,upperLimit=0.1,picTitle="all stars mean o-c in RA or DE",\
             vLinePos=100,medianText="real",moonColors= ['b','g','r','c','y','m'],scatterColor="black",dataLabel="",\
             hInfoColor="green",plotStatMore=True,pSize=8,hideTargetText=False,yLabel='$\\langle O-C \\rangle$ (mas)'):
    unitScale=1000
    lowerLimit*=unitScale
    upperLimit*=unitScale
    allStarMeanMedian=unitScale*stat.median(np.array(statListSorted,dtype=float)[ :,3]) 
    ax.scatter(np.array(statListSorted,dtype=float)[:,0],unitScale*np.array(statListSorted,dtype=float)[:,3],color=scatterColor,s =pSize*fontScale,label=dataLabel )
    if plotStatMore:
        'Highlighting the average O-C of natural satellites, for gaia data, targets with IDs less than 1e10 are considered as targets'
        moonNo=0
        if not hideTargetText:
            for moonIndex in range(len(statListSorted)):
                if isTarget(statListSorted[moonIndex][2]):
                    ax.plot(float(statListSorted[moonIndex][0]),unitScale*float(statListSorted[moonIndex][3]), moonColors[moonNo%len(moonColors)]+'.')#,markersize=pSize*fontScale/3.6)
                    ax.text(float(statListSorted[moonIndex][0]) ,unitScale*(float(statListSorted[moonIndex][3])+0.02*moonNo),\
                            str(round(unitScale*float(statListSorted[moonIndex][3]),int(4-math.log10(unitScale)))).ljust(int(6-math.log10(unitScale)),"0"), \
                            fontsize=25*fontScale,color=moonColors[moonNo%len(moonColors)])
                    moonNo+=1
        
        if vLinePos<99:
            ax.axvline(vLinePos, color='black',linestyle='-', lw=2.5*fontScale )
        if medianText=="real":
            ax.axhline(allStarMeanMedian, color=hInfoColor,linestyle='--', lw=2.5*fontScale )
            ax.text(ax.get_xlim()[1]+0.1, allStarMeanMedian, str(round(allStarMeanMedian,4)).ljust(int(6-math.log10(unitScale)),"0"), fontsize=24*fontScale,color=hInfoColor) 
        elif medianText=="mid":
            ax.text(ax.get_xlim()[1]+0.1, (lowerLimit+upperLimit)/2, str(round(allStarMeanMedian,4)).ljust(int(6-math.log10(unitScale)),"0"), fontsize=24*fontScale,color="black") 
        
        
    plt.title(picTitle,  fontsize=24*fontScale)
    ax.set_xlabel('$G$ mag',size =30*fontScale)
    ax.set_ylabel(yLabel,size =26*fontScale)
    ax.tick_params(labelsize=24*fontScale,length=14*fontScale,width=2.5*fontScale)
    x1, x2 = ax.get_xlim()
    ax.set_xlim(math.floor(x1*2)/2, math.ceil(x2*2)/2)
    
    ax.set_ylim(lowerLimit,upperLimit)
    y1, y2 = ax.get_ylim()
    ax.set_yticks(np.arange(y1,y2+0.00001,unitScale*0.02*pcRangeFactor)) 
    ax.xaxis.set_minor_locator(MultipleLocator(1))  
    ax.yaxis.set_minor_locator(MultipleLocator(unitScale*0.01*pcRangeFactor))
    ax.tick_params(which= 'minor',length=7*fontScale,width=2.5*fontScale)   
    if upperLimit-lowerLimit>0.2*unitScale:
        ax.set_yticks(np.arange(lowerLimit,upperLimit+0.0001,0.05*unitScale))
    
    
    ax_c = ax.twinx() 
    y1, y2 = ax.get_ylim()
    ax_c.set_ylim(y1, y2) 
    ax_c.set_yticks(np.arange(y1,y2+0.00001,unitScale*0.05*pcRangeFactor)) 
    ax_c.tick_params(labelsize=0,length=14*fontScale,width=2.5*fontScale,direction ='in')
    ax_c.yaxis.set_minor_locator(MultipleLocator(unitScale*0.01*pcRangeFactor))
    ax_c.tick_params(which= 'minor',length=7*fontScale,width=2.5*fontScale,direction ='in')
    ax_c.figure.canvas.draw()   
    
    ax_d = ax.twiny() 
    x1, x2 = ax.get_xlim()
    ax_d.set_xlim(x1, x2)
    ax_d.tick_params(labelsize=0*fontScale,length=14*fontScale,width=2.5*fontScale,direction ='in')

    ax_d.xaxis.set_minor_locator(MultipleLocator(1))
    ax_d.tick_params(which= 'minor',length=7*fontScale,width=2.5*fontScale,direction ='in')
    ax_d.figure.canvas.draw()   
    for sps in ax.spines.values():
        sps.set_linewidth(2.5*fontScale)

def plotSD(ax,statListSorted,plsqPrecisionCurve,IS_PAPER_PLOT=False,PAPER_FIG_TITLE=""\
           ,curIterExcludeStarDict=dict(),pSize=8, upperLimit=0.1,moonColors= ['b','g','r','c','y','m']\
           ,dataLabel="",scatterColor="black",textColorOff=0,textYOff=0):
    
    plt.title(PAPER_FIG_TITLE,  fontsize=24*fontScale)
    ax.tick_params(labelsize=24*fontScale,length=14*fontScale,width=2.5*fontScale)  
    unitScale=1000
    
    monorTicksInterval=math.ceil(math.ceil(unitScale*upperLimit/100)/200*unitScale/10)/100
    ax.xaxis.set_minor_locator( MultipleLocator(1)) 
    ax.yaxis.set_minor_locator( MultipleLocator(unitScale*monorTicksInterval*pcRangeFactor))
    ax.scatter(np.array(statListSorted,dtype=float)[:,0],unitScale*np.array(statListSorted,dtype=float)[:,1],color=scatterColor,s =pSize*fontScale,label=dataLabel)
    if not IS_PAPER_PLOT:
        moonNo=0
        for moonIndex in range(len(statListSorted)):
            'Highlighting the average O-C of natural satellites, for gaia data, targets with IDs less than 1e10 are considered as targets'
            if isTarget(statListSorted[moonIndex][2]):
                ax.plot(float(statListSorted[moonIndex][0]),unitScale*float(statListSorted[moonIndex][1]), moonColors[moonNo%len(moonColors)]+'.',markersize =pSize*fontScale/3.6)
                ax.text(float(statListSorted[moonIndex][0]), unitScale*(textYOff+float(statListSorted[moonIndex][1])+0.02*moonNo),\
                        str(round(unitScale*float(statListSorted[moonIndex][1]),int(4-math.log10(unitScale)))).ljust(int(6-math.log10(unitScale)),"0"), \
                        fontsize=25*fontScale,color=moonColors[moonNo%len(moonColors)])#color=scatterColor)#
                
                moonNo+=1

        rejectedStarInfoList=[]
        for curIterStar in statListSorted:
            if curIterExcludeStarDict.get(curIterStar[2])!=None:
                rejectedStarInfoList.append([float(curIterStar[0]),float(curIterStar[1])]) 
        if rejectedStarInfoList!=[] :
            ax.plot(np.array(rejectedStarInfoList)[:,0],unitScale*np.array(rejectedStarInfoList)[:,1], 'xr', markersize=4)

    ax.set_xlabel('$G$ mag',size =30*fontScale)
    ax.set_ylabel('SD (mas)',size =30*fontScale) 
    ax.set_ylim(0, unitScale*upperLimit*pcRangeFactor)  
#    ax.set_xlim(xLeft, xRight)  
    x1, x2 = ax.get_xlim()
    ax.set_xlim(math.floor(x1*2)/2, math.ceil(x2*2)/2)
    ax.set_yticks(np.arange(0,unitScale*upperLimit*pcRangeFactor+0.00001,unitScale*monorTicksInterval*3*pcRangeFactor))
 
    # ax.set_ylim( 1, 400)
    # ax.set_yscale('log') 
    
    ax.tick_params(which= 'minor',length=7*fontScale,width=2.5*fontScale)
     
    ax_c = ax.twinx() 
    y1, y2 = ax.get_ylim()
    ax_c.set_ylim(y1, y2) 
    
    ax_c.set_yticks(np.arange(y1,y2+0.00001,unitScale*monorTicksInterval*3*pcRangeFactor)) 
    ax_c.tick_params(labelsize=0,length=14*fontScale,width=2.5*fontScale,direction ='in')
    ax_c.yaxis.set_minor_locator(MultipleLocator(unitScale*monorTicksInterval*pcRangeFactor))
    ax_c.tick_params(which= 'minor',length=7*fontScale,width=2.5*fontScale,direction ='in')
    ax_c.figure.canvas.draw()   
    
    ax_d = ax.twiny() 
    x1, x2 = ax.get_xlim()
    ax_d.set_xlim(x1, x2)
    ax_d.tick_params(labelsize=0*fontScale,length=14*fontScale,width=2.5*fontScale,direction ='in')
    ax_d.xaxis.set_minor_locator(MultipleLocator(1))
    ax_d.tick_params(which= 'minor',length=7*fontScale,width=2.5*fontScale,direction ='in')
    ax_d.figure.canvas.draw()   
    for sps in ax.spines.values():
        sps.set_linewidth(2.5*fontScale)

 
#    ax.legend(loc=(0.05,0.76), shadow=True, fontsize=24*fontScale)  
    if len(statListSorted)!=0 and len(plsqPrecisionCurve.x)!=0:
        fittedCurveX=np.arange(np.array(statListSorted,dtype=float)[:,0].min()  ,np.array(statListSorted,dtype=float)[:,0].max() ,0.05)
        ax.plot(fittedCurveX,unitScale*fittingFunc(plsqPrecisionCurve.x,fittedCurveX),color='tab:orange',linewidth=5.5 *fontScale ,label="Fitted Curve")
    ax.legend(loc='best', shadow=True, fontsize=24*fontScale) 
    return math.ceil(x2*2)/2
def plotPrecisionResiduals(ax,statListSorted,plsqPrecisionCurve,PAPER_FIG_TITLE=""):    
    plt.title(PAPER_FIG_TITLE,  fontsize=24*fontScale)
    ax.tick_params(labelsize=24*fontScale,length=14*fontScale,width=2.5*fontScale) 

    ax.xaxis.set_minor_locator( MultipleLocator(1)) 
    ax.yaxis.set_minor_locator( MultipleLocator(0.01)) 
    ax.scatter(np.array(statListSorted,dtype=float)[:,0],
               np.array(statListSorted,dtype=float)[:,1]-fittingFunc(plsqPrecisionCurve.x,np.array(statListSorted,dtype=float)[:,0]),
               color="red",s =25*fontScale ,label="Residual",marker="x") 
    ax.axhline(0, color='black',linestyle='-', lw=2.5*fontScale )
    ax.set_xlabel('$G$ mag',size =30*fontScale)
    ax.set_ylabel('Residual (arcsec)',size =30*fontScale) 
    ax.set_ylim(-0.06, 0.06 )  
    ax.set_xlim(7.5, 20)  
    ax.set_yticks(np.arange(-0.06,0.061,0.02))
    x1, x2 = ax.get_xlim()
    ax.set_xlim(math.floor(x1*2)/2, math.ceil(x2*2)/2)
    
    ax.tick_params(which= 'minor',length=7*fontScale,width=2.5*fontScale)
    
    ax_c = ax.twinx() 
    y1, y2 = ax.get_ylim()
    ax_c.set_ylim(y1, y2) 
    ax_c.set_yticks(np.arange(y1,y2+0.00001,0.02)) 
    ax_c.tick_params(labelsize=0,length=14*fontScale,width=2.5*fontScale,direction ='in')
    ax_c.yaxis.set_minor_locator(MultipleLocator(0.01))
    ax_c.tick_params(which= 'minor',length=7*fontScale,width=2.5*fontScale,direction ='in')
    ax_c.figure.canvas.draw()   
    
    ax_d = ax.twiny() 
    x1, x2 = ax.get_xlim()
    ax_d.set_xlim(x1, x2)
    ax_d.tick_params(labelsize=0*fontScale,length=14*fontScale,width=2.5*fontScale,direction ='in')
    ax_d.xaxis.set_minor_locator(MultipleLocator(1))
    ax_d.tick_params(which= 'minor',length=7*fontScale,width=2.5*fontScale,direction ='in')
    ax_d.figure.canvas.draw()   
    for sps in ax.spines.values():
        sps.set_linewidth(2.5*fontScale)


def curveFittingFuncResidual(x,a,b,c,d):
    return   b + (a-b)/(1 + math.e**((x-c)/d))

def precisionCurve(outputList,excludeStarDict,IS_PAPER_PLOT=False,PAPER_FIG_TITLE=""\
                   ,brightStarMag=15,No_Weight_In_reduction=False,unitScale=1000): 
    statList=[]
    brightStarStatList=[]
    statListRA=[]
    statListDE=[]
    plsqPrecisionCurve=plsq()
    
    upperLimitRA=0.1
    upperLimitDE=0.1
    for line in outputList:
        line=line.strip() 
        fields= re.split(r'\s+',line)
        'mag, sd, id, mean(o-c)'
        statList.append([fields[1],(float(fields[4])**2+float(fields[5])**2)**0.5,fields[0],(float(fields[2])**2+float(fields[3])**2)**0.5])
        statListRA.append([fields[1], fields[4],fields[0], fields[2]])
        statListDE.append([fields[1], fields[5],fields[0], fields[3]])
        if float(fields[4])>upperLimitRA:
            upperLimitRA=float(fields[4])
        if float(fields[5])>upperLimitDE:
            upperLimitDE=float(fields[5])
        if float(fields[1])<brightStarMag:
            brightStarStatList.append([fields[1],(float(fields[4])**2+float(fields[5])**2)**0.5,\
                                       fields[0],(float(fields[2])**2+float(fields[3])**2)**0.5])

    
    upperLimitRA+=0.005
    upperLimitDE+=0.005
    if brightStarStatList==[]:
        print("The parameter setting for brightStarMag is too small, and there are no stars brighter than that magnitude!") 
        sys.exit()  
         
        
    '平均O-C的显示范围'
    OMCLowerRange=-0.08
    OMCUpperRange=0.08
    statListSortedRA=sorted(statListRA,key=lambda a:float(a[0])) 
    statListSortedDE=sorted(statListDE,key=lambda a:float(a[0])) 
    
    fig,ax=plt.subplots() 
    plotMean(ax,statListSortedRA,OMCLowerRange,OMCUpperRange,\
             picTitle="all stars mean o-c in RA ("+str(OMCLowerRange)+","+str(OMCUpperRange)+")")
    
    fig,ax=plt.subplots() 
    plotMean(ax,statListSortedDE,OMCLowerRange,OMCUpperRange,\
             picTitle="all stars mean o-c in DE ("+str(OMCLowerRange)+","+str(OMCUpperRange)+")")
     
    plt.show()
     
    statListSorted=sorted(statList,key=lambda a:float(a[0]))  
    brightStarMeanMedian=unitScale*stat.median(np.array(brightStarStatList,dtype=float)[ :,3]) 
    
    fig,ax=plt.subplots() 
    plotMean(ax,statListSorted,picTitle="all stars mean O$-$C",vLinePos=brightStarMag,medianText=False)  
    ax.axhline(brightStarMeanMedian, color='red',linestyle='--', lw=2.5*fontScale )
    ax.text(ax.get_xlim()[1]+0.1, brightStarMeanMedian, str(round(brightStarMeanMedian,int(4-math.log10(unitScale)))).ljust(int(6-math.log10(unitScale)),"0"), fontsize=24*fontScale,color='red')    
    plt.show()
    
     
     
    curIterExcludeStarDict=dict()
    
    meanStdList=genExcludeStarDict(statListSorted,curIterExcludeStarDict)
    genExcludeStarDict(sorted(statListRA,key=lambda a:float(a[0]))[2:len(statListRA)],curIterExcludeStarDict)
    genExcludeStarDict(sorted(statListDE,key=lambda a:float(a[0]))[2:len(statListDE)],curIterExcludeStarDict) 
    excludeStarDict=dict(curIterExcludeStarDict, **excludeStarDict) 
    
    
    
    if statListSorted!=[] and not No_Weight_In_reduction:
        
        if len(meanStdList)>5:
            
            if len(statListSorted)>20 and float(statListSorted[-10][0])<int(np.array(meanStdList,dtype=float)[:,0][-1]):
                endIndex=-1
            else:
                endIndex=len(np.array(meanStdList,dtype=float)[:,0])
            boltzmannX=np.array(meanStdList,dtype=float)[1:endIndex,0]
            boltzmannY=np.array(meanStdList,dtype=float)[1:endIndex,1]
        else:
            boltzmannX=np.array(meanStdList,dtype=float)[:,0] 
            boltzmannY=np.array(meanStdList,dtype=float)[:,1] 
            
            
        
        p=[boltzmannY[1], boltzmannY[-1],float(statListSorted[-1][0]), float(statListSorted[-1][0])-4]
        plsqPrecisionCurve= least_squares(fittingFuncResidual, p, args=(boltzmannX, boltzmannY)) 

        
        popt, pcov = curve_fit(curveFittingFuncResidual,xdata=boltzmannX,ydata= boltzmannY,p0=plsqPrecisionCurve.x,maxfev=100000 )
        perr = np.sqrt(np.diag(pcov))
        print("perr:",perr,"      popt:",popt)
        print("meanStd plsqPrecisionCurve.x:", plsqPrecisionCurve.x )
        plsqPrecisionCurve.x=popt
    
    
        if  plsqPrecisionCurve.x[0]<0.001 or plsqPrecisionCurve.x[2]<14\
            or plsqPrecisionCurve.x[1]*plsqPrecisionCurve.x[-1]<plsqPrecisionCurve.x[0]*plsqPrecisionCurve.x[-1]\
        or any(plsqPrecisionCurve.x[:3])<0:
            print("The precision curve may be incorrect, and the weighted solution for the film constant will not be used!")
            plsqPrecisionCurve=plsq()
    else: 
        plsqPrecisionCurve=plsq()
#        plsqPrecisionCurve=[]
     
    brightStarSDMedian=unitScale*stat.median(np.array(brightStarStatList,dtype=float)[ :,1]) 
    fig,ax=plt.subplots()  
    plotSD(ax,statListSortedRA,plsq(),IS_PAPER_PLOT,PAPER_FIG_TITLE,curIterExcludeStarDict,upperLimit=upperLimitRA,scatterColor="black",dataLabel="RA")
    plotSD(ax,statListSortedDE,plsq(),IS_PAPER_PLOT,PAPER_FIG_TITLE,curIterExcludeStarDict,upperLimit=upperLimitDE,scatterColor="red",dataLabel="DE",textYOff=0.03)
    plt.show()
    
    
    fig,ax=plt.subplots() 
    xRight=plotSD(ax,statListSorted,plsqPrecisionCurve,IS_PAPER_PLOT,PAPER_FIG_TITLE,curIterExcludeStarDict, upperLimit=(upperLimitRA**2+upperLimitDE**2)**0.5)
    ax.axhline(brightStarSDMedian, color='red',linestyle='--', lw=2.5*fontScale )
    ax.axvline(brightStarMag, color='black',linestyle='-', lw=2.5*fontScale )
    ax.text(xRight+0.1, brightStarSDMedian, str(round(brightStarSDMedian,int(4-math.log10(unitScale)))).ljust(int(6-math.log10(unitScale)),"0"), fontsize=24*fontScale,color='red')
    
    plt.show()
    
    for moonIndex in range(len(statListSorted)): 
        if isTarget(statListSorted[moonIndex][2]):
            print("target",statListSorted[moonIndex][2]," precision：", str(round(float(statListSortedRA[moonIndex][1]),3)), str(round(float(statListSortedDE[moonIndex][1]),3)))
         
    curIterExcludeStarDict=dict()
       
    
    ax.legend(loc=(0.05,0.76), shadow=True, fontsize=24*fontScale)
    
    plt.show() 
    if statListSorted!=[] and plsqPrecisionCurve!=[] and IS_PAPER_PLOT:
        fig,ax=plt.subplots() 
        plotPrecisionResiduals(ax,statListSorted,plsqPrecisionCurve,PAPER_FIG_TITLE )
        plt.show()
    return plsqPrecisionCurve,excludeStarDict
    
 
    
if __name__=="__main__": 
    sys.path.pop()