'''
Descripttion: 
version: 
Author: zpliu
Date: 2021-07-05 23:33:26
LastEditors: zpliu
LastEditTime: 2021-07-06 11:23:48
@param: 
'''
'''
inputData:
    align region
    Ghir_A01        72995   76680   Ghir_D01        38281   42200   Ghir_A01G000040*-       Ghir_D01G000060*+       PromoterRegion  262472
    Ghir_A01        84826   96240   Ghir_D01        76408   79011   Ghir_A01G000070*+       Ghir_D01G000110*+       PromoterRegion  97892
    Ghir_A01        86191   88883   Ghir_D01        65254   67887   Ghir_A01G000080*+       Ghir_D01G000100*+       PromoterRegion  0
    match region:
    Ghir_A01        72994   73405   Ghir_D01        41741   42152   Ghir_A01G000040*-       Ghir_D01G000060*+       sameDirection   Match
    Ghir_A01        73405   73405   Ghir_D01        42562   42563   Ghir_A01G000040*-       Ghir_D01G000060*+       sameDirection   Insert
    Ghir_A01        73405   73662   Ghir_D01        42307   42564   Ghir_A01G000040*-       Ghir_D01G000060*+       sameDirection   Match
return: 
    diversity region:
    Ghir_A01        76151   76680   Ghir_A01G000040*-
    Ghir_D01        38281   41740   Ghir_D01G000060*+
    Ghir_D01        42153   42200   Ghir_D01G000060*+
'''
import sys 
import pybedtools
import pandas as pd 
from tempfile import NamedTemporaryFile
#############################################
#根据promoter、code区域和downstream区域的比对情况
#使用bedtools将总的区域substract比对好的区域
#############################################

def getDiversityRegion(alignInterVal,matchInterval):
    out=""
    tmpFile=NamedTemporaryFile(mode='w+t',encoding='utf-8')
    border1=[int(i[1]) for i in matchInterval ]
    border2=[int(i[2]) for i in matchInterval ]
    border=border1+border2
    DiversityInterval=alignInterVal.subtract(matchInterval.sort().merge())
    if DiversityInterval:
        DiversityInterval.saveas(tmpFile.name)
        DiversityDataFrame=pd.read_csv(tmpFile.name,header=None,index_col=None,sep="\t")
        for value  in DiversityDataFrame.values:
            if value[1]==value[2] and value[1] in border:
                pass
            elif value[1] in border and value[2] in border:
                out+="\t".join([value[0],str(value[1]+1),str(value[2]-1),value[3]])+"\n"
            elif value[1] in border and value[2] not in border:
                out+="\t".join([value[0],str(value[1]+1),str(value[2]),value[3]])+"\n"
            else:
                out+="\t".join([value[0],str(value[1]),str(value[2]-1),value[3]])+"\n"
        return out
    else:
        #! no diversity Interval
        return None



Alignment_region=pd.read_csv(sys.argv[1],header=None,index_col=None,sep="\t")
matchRegion=pd.read_csv(sys.argv[2],header=None,index_col=None,sep="\t")

Diverregion=[]
for index in range(0,Alignment_region.shape[0],1):
    #!获取比对的区间
    AtInterVal=pybedtools.BedTool("\t".join([
        Alignment_region.iloc[index,0],
        str(Alignment_region.iloc[index,1]),
        str(Alignment_region.iloc[index,2]),
        Alignment_region.iloc[index,6]
    ]),from_string=True)
    DtInterVal=pybedtools.BedTool("\t".join([
        Alignment_region.iloc[index,3],
        str(Alignment_region.iloc[index,4]),
        str(Alignment_region.iloc[index,5]),
        Alignment_region.iloc[index,7]
    ]),from_string=True)
    AtgeneId=Alignment_region.iloc[index,6]
    DtgeneId=Alignment_region.iloc[index,7]
    #! 获取At和Dt成对的区间
    AtmatchRegionInterval=pybedtools.BedTool("\n".join(
        matchRegion.loc[matchRegion[6]==AtgeneId][[0,1,2,6]].apply(
            lambda x:"\t".join([str(i) for i in x]),axis=1)),
            from_string=True
            )
    DtmatchRegionInterval=pybedtools.BedTool(
        "\n".join(matchRegion.loc[matchRegion[7]==DtgeneId][[3,4,5,7]].apply(
            lambda x:"\t".join([str(i) for i in x]),axis=1)),
            from_string=True
    )
    #获取特有的区间 
    #! 有的基因完全匹配上了
    tmpA=getDiversityRegion(AtInterVal,AtmatchRegionInterval)
    tmpD=getDiversityRegion(DtInterVal,DtmatchRegionInterval)
    if tmpA:
        Diverregion.append(tmpA)
    if tmpD:
        Diverregion.append(tmpD)
with open(sys.argv[3],'w') as File:
    for line in Diverregion:
        File.write(line)
    