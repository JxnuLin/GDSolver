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
HeadInfo=namedtuple('HeadInfo',['FILENAME','STARCOUNT','YEAR','MONTH','DAY','HOUR','MINUTE','SECOND','EXPOSURE',   'JD',   'a',   'b',   'c',   'd'])
StarInfo=namedtuple('StarInfo',['ID','x','y','RA','DE','Mag','pmRA','pmDE','ksi','eta',     'dRA','dDE'])
PpStarInfo=namedtuple('PpStarInfo',['ID','x','y','RA','DE','Mag','pmRA','pmDE','ksi','eta', 'dRA','dDE' ,'refStarID','sepa'])
 

def read_mathcedStar(allHeadInfoList,allFrameStarList,fileName,region=[0,0,0,0],filePath=os.path.dirname(os.path.dirname(__file__))+'/outputFiles/'):
    "读取匹配结果的文件， 保存在allHeadInfoList,allFrameStarList中，其中allHeadInfoList中每个元素是一行头文件信息，allFrameStarList中每个元素是一帧中匹配到的所有星的列表，其中的每个元素表示的是一颗匹配星信息（参考匹配结果文本的格式）"
    curFrameStarList=[] 
    curFrameStarCount=0
    curHeadInfo= HeadInfo('FILENAME', 'STARCOUNT', 'YEAR', 'MONTH', 'DAY', 'HOUR', 'MINUTE', 'SECOND',   'EXPOSURE',   'JD',   'a',   'b',   'c',   'd')
    with open(filePath+fileName, 'rt') as file:
        for line in file:
            line=line.strip() 
            fields= re.split(r'\s+',line) 

            
            if len(fields) ==14:
                
                """当前行读取的是信息头"""
                if curFrameStarCount!=0:
                    """处理上一帧的信息"""
                    newHeadInfo =HeadInfo(*curHeadInfo[:1],curFrameStarCount,*curHeadInfo[2:14])
                    if curFrameStarCount!=int(curHeadInfo.STARCOUNT):
                        print("图像："+curHeadInfo.FILENAME+"行数和星数不匹配！读入的结果中已将星数修正为行数！")
                    allHeadInfoList.append(newHeadInfo)
                    allFrameStarList.append(curFrameStarList)
                    curFrameStarList=[]
                curFrameStarCount=0
                curHeadInfo= HeadInfo(*fields)
                #print(curHeadInfo)
            elif len(fields)==12:
                if region[2]!=0 and (float(fields[1])<region[0] or float(fields[1])>region[2] or float(fields[2])<region[1] or float(fields[2])>region[3]):
                    continue
                "读取的匹配文件中星信息的每一行数据列数，根据不同的匹配文件，该参数可能不同"
                curFrameStarList.append(StarInfo(*fields))
                curFrameStarCount+=1
        """添加最后一幅图的信息"""
        newHeadInfo =HeadInfo(*curHeadInfo[:1],curFrameStarCount,*curHeadInfo[2:14])
        allHeadInfoList.append(newHeadInfo)
        allFrameStarList.append(curFrameStarList)
        if len(curFrameStarList)==0:
            print("未读取到有效行，可能是输入文件格式有误！")