# -*- coding: utf-8 -*-
"""
Created on Mon Nov 27 16:25:19 2017

@author: Administrator
"""
from scipy import optimize 
import numpy as np
import cmath
def func(param,x,y,model,mag=[],fitMagTermOrder=0): 
    
    if model==1:
        combineLen=x.size
        coeffZerosOnesArray=np.zeros(combineLen)
        coeffZerosOnesArray[int(combineLen/2):combineLen]=np.ones(int(combineLen/2))
        coeffOnesZerosArray=np.ones(combineLen)
        coeffOnesZerosArray[int(combineLen/2):combineLen]=np.zeros(int(combineLen/2))
        return param[0]*coeffOnesZerosArray+param[1]*x+param[2]*y  +param[3]*coeffZerosOnesArray

          
    else:      
        result=genPolynomial(param,x,y,model-1)
        if mag!=[]:
            for magTermNum in range(fitMagTermOrder):
                result+=param[int(model*(model+1)/2)+magTermNum]*np.array(mag)**(magTermNum+1)
        return  result
def genPolynomial(p,x,y,order):
    if order>1:
        lowOrderFunc=genPolynomial(p,x,y,order-1)
        for itemI in range(order+1):
            lowOrderFunc+=p[int(order*(order+1)/2+itemI)]*x**itemI*y**(order-itemI)
        return lowOrderFunc
    else:
        return p[0]+p[1]*x+p[2]*y 
    
    
    
def residuals(p,x,y,obj,model,weights,mag=[],fitMagTermOrder=0):  
    return  weights*(func( p,x,y,model,mag,fitMagTermOrder)- obj)


def fit_model(x,y,obj,obj2,model,weights=[]):
    
    if model==1:
        p0=[0,0,0 ,0]
        newX=np.resize(x,x.size*2)
        newX[x.size:x.size*2]= -y
        newY=np.resize( y,y.size*2)
        newY[y.size:y.size*2]=x
        newObj=np.resize(obj,obj.size*2)
        newObj[obj.size:obj.size*2]=obj2
        if len(weights)==0:
            weights=np.ones(len(newX))
        plsq = optimize.least_squares(residuals,  p0, args=(newX,newY,newObj,model,weights))
        "In the result of four constant fitting, the four parameters are c, a, b, and d in ksi=a * x+b * y+c, eta=b * x-a * y+d"
        return plsq
    else:
        
        if len(weights)==0:
            weights=np.ones(len(x))
        if model<6:
            if model==2:
                #6常数
                p0=[0,0,0]
            elif model==3:
                p0=[0,0,0 ,0,0,0]
            elif model==4:
                p0=[0,0,0 ,0,0,0,0,0,0,0]
            elif model==5:
                p0=[0,0,0 ,0,0,0,0,0,0,0,0,0,0,0,0]
            plsq = optimize.least_squares(residuals, p0, args=(x,y,obj,model,weights)) 
            plsq2 = optimize.least_squares(residuals, p0, args=(x,y,obj2,model,weights))  
        else: 
            plsqtemp = optimize.least_squares(residuals,[0,0,0 ,0,0,0,0,0,0,0,0,0,0,0,0], args=(x,y,obj,5,weights)) 
            plsqtemp2 = optimize.least_squares(residuals,[0,0,0 ,0,0,0,0,0,0,0,0,0,0,0,0], args=(x,y,obj2,5,weights)) 
            
            p0=np.append(plsqtemp.x[:3],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])
            p02=np.append(plsqtemp2.x[:3],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])
            if model==6:
                pass
            elif model==7:
                p0=np.append(p0,[0,0,0,0,0,0,0])
                p02=np.append(p02,[0,0,0,0,0,0,0])
            plsq = optimize.least_squares(residuals, p0, args=(x,y,obj,model,weights)) 
            plsq2 = optimize.least_squares(residuals, p02, args=(x,y,obj2,model,weights))  
        return plsq,plsq2
     
'y = A2 + (A1-A2)/(1 + exp((x-x0)/dx))'
def precisionPremiumCurve(p,x): 
    'Force parameter x0 to 0 to intercept the upper half of the Boltzmann curve'
    return   p[1] + (p[0]-p[1])/(1 + cmath.e**(x/p[2]))
#    return p[0]+p[1]*x+p[2]*x**2+p[3]*x**3+p[4]*x**4
    

def sigmoidCurve(p,x):
    import math
    return   p[1] + (p[0]-p[1])/(1 + math.e**((x-p[2])/p[3])) 