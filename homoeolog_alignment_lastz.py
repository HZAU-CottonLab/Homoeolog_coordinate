'''
Descripttion: 
version: 
Author: zpliu
Date: 2021-07-04 09:52:24
LastEditors: zpliu
LastEditTime: 2021-07-04 11:31:26
@param: 
'''
import logging
import pandas as pd 
import numpy as np
import pybedtools
import re
import pysam
from tempfile import NamedTemporaryFile
import os 
from lastz_script.lastz import run_lastz_result 
import sys 
logger=logging.getLogger()
logger.setLevel(logging.INFO)
logger.info("loading module...")



####################################################
#! 根据bnMapper过滤后的结果，将整个同源基因区域的对应关系补充完整
###################################################
def get_flankbnMApper_region(bnMapperregion,geneRegion):
    '''
    according the bnMapper region to get flank interval
    #! don't consider the stand; gene may transfor 
    args:
        -bnMapperregion: bnMapper gene between homoeolog
        -geneRegion: 2k+gene body region
    return:
        -lastzRegion: @list
    '''
    out=[bnMapperregion]
    #! the remains bp compared to gene boundary
    #! withou left  get 0,else 延伸两者中最长的长度
    At=[geneRegion[1]-bnMapperregion[1],geneRegion[2]-bnMapperregion[2]]
    Dt=[geneRegion[4]-bnMapperregion[4],geneRegion[5]-bnMapperregion[5]]
    #########################
    # left
    #########################
    if not At[0] and not  Dt[0]:
        pass
        # out.append([geneRegion[0],geneRegion[1],bnMapperregion[1]-1,geneRegion[3],geneRegion[4],bnMapperregion[4]-1])
    else:
        interValLength=min([geneRegion[1]-bnMapperregion[1],geneRegion[4]-bnMapperregion[4]])
        out.append([geneRegion[0],bnMapperregion[1]+interValLength,bnMapperregion[1]-1,geneRegion[3],bnMapperregion[4]+interValLength,bnMapperregion[4]-1])

    ############################
    # right
    ############################
    if not At[1] and not Dt[1]:
        pass
        # out.append([geneRegion[0],bnMapperregion[2]+1,geneRegion[2],geneRegion[3],bnMapperregion[5]+1,geneRegion[5]])
    else :
        interValLength=max([geneRegion[2]-bnMapperregion[2],geneRegion[5]-bnMapperregion[5]])
        out.append([geneRegion[0],bnMapperregion[2]+1,bnMapperregion[2]+interValLength,geneRegion[3],bnMapperregion[5]+1,bnMapperregion[5]+interValLength])
    return out

#########################################
# 获取lastz需要比对的区域    
#########################################    
filter_bnMapper=pd.read_csv("/data/cotton/zhenpingliu/LZP_fiberFullPopulationRNAseq/04homologExpressionPattern/sub-genome-chain/homoeologVariant/gene_region/homoeolog_bnMapper_filter.txt",header=None,index_col=None,sep="\t")
##############################
# extract genome sequence
##############################
genomeFile='/data/cotton/MaojunWang/WMJ_fiberFullPopulationRNAseq/MappingFPKM/Ghir_Genome_Index/Ghirsutum_genome.fasta'
genomeObject=pysam.FastaFile(genomeFile)
###############################
# homoeolog gene ID 
###############################
homologGeneID_stand=pd.read_csv(sys.argv[1],header=None,index_col=None,sep="\t")
out=[]
for homoeologGeneId,stand in homologGeneID_stand[[3,4]].values:
    # print(homoeologGeneId)
    homoeologGeneInfo=homologGeneID_stand.loc[homologGeneID_stand[3]==homoeologGeneId]
    homoeolog_bnMapperInfo=filter_bnMapper.loc[filter_bnMapper[3]==homoeologGeneId+"*"+stand]
    if homoeolog_bnMapperInfo.empty:
        print("this homoeolog gene do'nt have chain flagment {}".format(homoeologGeneId))
        continue
    #! get bnMApper pairs region
    Atchrom,Atstart,Atend=homoeolog_bnMapperInfo.iloc[0,0],homoeolog_bnMapperInfo.sort_values(by=[1]).iloc[0,1],homoeolog_bnMapperInfo.sort_values(by=[1]).iloc[-1,2]
    Dtchrom,Dtstart,Dtend=homoeolog_bnMapperInfo.iloc[0,4],homoeolog_bnMapperInfo.sort_values(by=[5]).iloc[0,5],homoeolog_bnMapperInfo.sort_values(by=[5]).iloc[-1,6]
    bnMapperRegion=[Atchrom,Atstart,Atend,Dtchrom,Dtstart,Dtend]
    geneRegion=[homoeologGeneInfo.iloc[0,0],homoeologGeneInfo.iloc[0,1],homoeologGeneInfo.iloc[0,2],homoeologGeneInfo.iloc[0,5],homoeologGeneInfo.iloc[0,6],homoeologGeneInfo.iloc[0,7]]
    #######################################
    ##get paired region
    #######################################
    AliginPairedRegion=get_flankbnMApper_region(bnMapperRegion,geneRegion)
    #! get mapper region
    out.append(run_lastz_result(AliginPairedRegion,homoeologGeneInfo,genomeObject))
    
with open(sys.argv[2],'w') as File :
    for line in out:
        File.write(line+"\n")
