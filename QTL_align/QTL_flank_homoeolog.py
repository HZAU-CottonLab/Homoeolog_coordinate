'''
Descripttion: 
version: 
Author: zpliu
Date: 2021-07-13 15:03:26
LastEditors: zpliu
LastEditTime: 2021-07-13 21:12:35
@param: 
'''
import pandas as pd 
import numpy as np 
import sys


homoeologCoordinate=pd.read_csv("homoeolog_coordinate.txt",header=None,index_col=None,sep="\t")
'''homoeolog_coordinate
    Ghir_A01        70889   74180   Ghir_A01G000040-Ghir_D01G000060
    Ghir_D01        40781   44323   Ghir_D01G000060-Ghir_A01G000040
    Ghir_A01        87326   100744  Ghir_A01G000070-Ghir_D01G000110
    Ghir_D01        78908   83636   Ghir_D01G000110-Ghir_A01G000070
'''
QTLcoordinate=pd.read_csv(sys.argv[1],header=None,index_col=None,sep="\t")
'''QTL coordinate File
    Ghir_A01        100011098       100011098
    Ghir_A01        10005837        10005837
    Ghir_A01        100086763       100086763
'''


out=[]
for QTL in QTLcoordinate.values:
    chrom=QTL[0]
    start=QTL[1]
    tmpData=homoeologCoordinate.loc[homoeologCoordinate[0]==chrom].sort_values(by=[1])
    tmpArry=np.sort(np.append(tmpData[1],start))
    #! QTL location in the chromsome
    QTLIndex=np.argwhere(tmpArry==start)[0][0]
    tmpDatalength=tmpData.shape[0]
    if QTLIndex==0:
        #! the first gene ID after QTL
        geneId1=tmpData.iloc[0,3]
        #* Ghir_A01G000040-Ghir_D01G000060
        chromsome1=geneId1.split("-")[0][0:8]+","
        chromsome2=geneId1.split("-")[1][0:8]+","
        gene1=geneId1.split("-")[0]
        gene2=geneId1.split("-")[1]
    elif QTLIndex==tmpDatalength:
        #! the last geneID  before QTL
        geneId1=tmpData.iloc[tmpDatalength-1,3]
        chromsome1=geneId1.split("-")[0][0:8]+","
        chromsome2=geneId1.split("-")[1][0:8]+","
        gene1=geneId1.split("-")[0]
        gene2=geneId1.split("-")[1]
    else:
        geneId1=tmpData.iloc[QTLIndex,3]
        geneId2=tmpData.iloc[QTLIndex-1,3]
        chromsome1=",".join(list(set([geneId1.split("-")[0][0:8],geneId2.split("-")[0][0:8],])))
        chromsome2=",".join(list(set([geneId1.split("-")[1][0:8],geneId2.split("-")[1][0:8],])))
        gene1=geneId1.split("-")[0]+","+geneId2.split("-")[0]
        gene2=geneId1.split("-")[1]+","+geneId2.split("-")[1]
    out.append((gene1,gene2,chromsome1,chromsome2))
geneId1,geneId2,chr1,chr2=zip(*out)    
QTLcoordinate[3]=geneId1
QTLcoordinate[4]=geneId2
QTLcoordinate[5]=chr1
QTLcoordinate[6]=chr2
QTLcoordinate.to_csv(sys.argv[2],header=False,index=False,sep="\t")

'''flank homoeolog geneId of QTL
Ghir_A01        164527  164527  Ghir_A01G000250,Ghir_A01G000220 Ghir_D01G000260,Ghir_D01G000230 Ghir_A01        Ghir_D01
Ghir_A01        255346  255346  Ghir_A01G000410,Ghir_A01G000380 Ghir_D01G000410,Ghir_D01G000380 Ghir_A01        Ghir_D01
Ghir_A01        269159  269159  Ghir_A01G000410,Ghir_A01G000380 Ghir_D01G000410,Ghir_D01G000380 Ghir_A01        Ghir_D01
Ghir_A01        277130  277130  Ghir_A01G000410,Ghir_A01G000380 Ghir_D01G000410,Ghir_D01G000380 Ghir_A01        Ghir_D01
'''
    


        
        
        
        
    

    


