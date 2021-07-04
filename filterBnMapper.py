'''
Descripttion: 
version: 
Author: zpliu
Date: 2021-07-02 21:57:39
LastEditors: zpliu
LastEditTime: 2021-07-02 23:27:38
@param: 
'''
import logging
import pandas as pd 
import numpy as np
import pybedtools
import re
import pysam 
import sys 
logger=logging.getLogger()
logger.setLevel(logging.INFO)
logger.info("loading module...")
##########################################
#get homolog paired region
#! 根据bnMapper的结果将lastz比对准确的区域对应起来，过滤不在同源基因区间的映射
#! 对于中间的间隔区域，同样使用muscle对应起来.
#! 由于muscle在比对的时候无法区分正负链
##########################################
def merge_filterInterval(IntervalList):
    '''
     #! get all interval and sorted 
     args:
        - IntervalList: BedTools list 
     return:
        - list region
    '''
    tmp=np.array([i[0][:] for i in IntervalList])
    return [tmp[0,0],tmp[0,1],tmp[-1,2]]
homologGeneID_stand=pd.read_csv("./gene_region/homoeolog_promoter_geneBody.bed",header=None,index_col=None,sep="\t")
##################################################
#filter the mappter Interval
##################################################
out=[]
with open(sys.argv[1],'r') as File:
    count=1
    for line in File:
        if count%10000==0:
            logger.info(count)
        count+=1
        line=line.strip("\n").split("\t")
        if line[4]=="None":
            #! without mapper region
            continue
        else:
            homoeologRegionData=homologGeneID_stand.loc[(homologGeneID_stand[3]==line[3].strip("*[+-]"))|(homologGeneID_stand[8]==line[3].strip("*[+-]"))]
            At_geneRegion="\t".join(homoeologRegionData.iloc[0,0:4].astype(str))+"*"+homoeologRegionData.iloc[0,4]
            Dt_geneRegion="\t".join(homoeologRegionData.iloc[0,5:-1].astype(str))+"*"+homoeologRegionData.iloc[0,-1]
            At_geneInterval=pybedtools.BedTool(At_geneRegion,from_string=True)
            Dt_geneInterval=pybedtools.BedTool(Dt_geneRegion,from_string=True)
            #! retain the interval intersect with homoeolog region 
            #! get Mapper interval from string 
            #! example 'Ghir_D10:38086-38275'
            getInterval=lambda x: pybedtools.BedTool(re.sub(r'[:-]',"\t",x),from_string=True)
            if re.match(r'^Ghir_A',line[0]):
                #! At request and mapper to Dt, filter the mapper
                #! a request may have more than one mapping
                #! todo(1): the mapprer may have a large interval between them
                mapper=Dt_geneRegion.split("\t")[-1]
                filterInterval=[getInterval(i) for i in line[4:] if getInterval(i).intersect(Dt_geneInterval)] 
            if re.match(r'^Ghir_D',line[0]):
                mapper=At_geneRegion.split("\t")[-1]
                #! Dt request and mapper to At, filter the mapper
                filterInterval=[getInterval(i) for i in line[4:] if getInterval(i).intersect(Dt_geneInterval)] 
            ###########################################
            # filter interval in the mapper gene region
            ###########################################
            if filterInterval and re.match("Ghir_D",mapper):
                #change order to At-Dt
                out.append("\t".join(line[0:4])+"\t"+"\t".join(merge_filterInterval(filterInterval))+"\t"+mapper+"\n")
            elif filterInterval and re.match("Ghir_A",mapper):
                #! change order to At-Dt
                out.append("\t".join(merge_filterInterval(filterInterval))+"\t"+mapper+"\t".join(line[0:4])+"\n")
            else:
                #! the mapper out of 2k+gene body
                continue
            # break
with open(sys.argv[2],'w') as File:
    for line in out:
        File.write(line)
            
            


