'''
Descripttion: 
version: 
Author: zpliu
Date: 2021-07-13 19:30:22
LastEditors: zpliu
LastEditTime: 2021-07-13 21:16:24
@param: 
'''
import logging
import pandas as pd
import numpy as np
import pybedtools
import re
import pysam
import sys 
import os 
from tempfile import NamedTemporaryFile
while True:
    #! load module 
    filePath=os.path.abspath(__file__)
    dirPath,fileName=os.path.split(filePath)
    modulePath=os.path.join(dirPath,"../../")
    sys.path.append(modulePath)
    from Homoeolog_coordinate.lastz_script.lastz import parse_lastz_cigra 
    break

def run_QTL_lastz(QTL_pairedRegion,genomeObject):
    '''run lastz to align the QTL region
    args:
        -QTL_pairedRegion: @pd.DataFrame
        -genomeObject: genome sequence file(@pysam)
    return:
        
    '''
    out = []
    Matchscore=[]
    for  pairedRegion in QTL_pairedRegion:
        #! align each region 
        AtFile = NamedTemporaryFile(mode='w+t', encoding='utf-8')
        DtFile = NamedTemporaryFile(mode='w+t', encoding='utf-8')
        if pairedRegion[1] >= pairedRegion[2] or pairedRegion[5] >= pairedRegion[6]:
            #! the Interval only 1bp
            print("Error with:\n",pairedRegion)
            continue
        try:
            #*index start with 0
            AtFile.write(
                ">At\n"+genomeObject.fetch(start=pairedRegion[1]-1, end=pairedRegion[2], region=pairedRegion[0]))
            DtFile.write(
                ">Dt\n"+genomeObject.fetch(start=pairedRegion[5]-1, end=pairedRegion[6], region=pairedRegion[4]))
        except ValueError:
            #! the paired region low than 0 when extend out of chromosome
            print('the region out of chromsome', end="\t")
            print(pairedRegion)
            continue 
        #! read and write model so back the point to first line
        AtFile.seek(0)
        DtFile.seek(0)
        # print(DtFile.read())
        # print(AtFile.read())
        ############################################
        # begin to lastz
        #! hspthresh=3000 filter diversity
        #! 保留片段比对之间的得分
        ############################################
        try:
            cigraFlag = os.popen(
                "/data/cotton/zhenpingliu/LZP_fiberFullPopulationRNAseq/software/lastz-1.04.03/lastz-distrib/bin/lastz --strand=both --step=10  --hspthresh=3000 --format=cigar {} {}".format(AtFile.name, DtFile.name)).read()
        except:
            #! if region out of chromsome the cigra will be ''
            print('the region out of chromsome', end="\t")
            print(pairedRegion)
            continue
        
        QTLchrom,QTLstart,QTLend,QTLId,targetChr,targetstart,targetend,targetId=pairedRegion
        homoeologGeneCoordinate=pd.DataFrame(
            [
                [QTLchrom,QTLstart,QTLend,QTLId.strip("*+"),'+',targetChr,targetstart,targetend,targetId.strip("*+"),'+'],
            ]
        )
        #! parse the cigra
        tmpPairRegion=[QTLchrom,QTLstart,QTLend,targetChr,targetstart,targetend]
        tmpout,tmpscore= parse_lastz_cigra(tmpPairRegion, cigraFlag,
                                 homoeologGeneCoordinate)
        out.append(tmpout+"\n")
        Matchscore.append("\t".join([str(i) for i in pairedRegion])+"\t"+str(tmpscore)+"\n")   
        
    return Matchscore,out

def test_run_QTL_lastz():
    genomeFile = '/data/cotton/MaojunWang/WMJ_fiberFullPopulationRNAseq/MappingFPKM/Ghir_Genome_Index/Ghirsutum_genome.fasta'
    genomeObject = pysam.FastaFile(genomeFile)
    #! QTL align coordinate
    QTLregion=pd.DataFrame(
        [
            ['Ghir_A01',164027,165027,'QTL1*+','Ghir_D01',145367,162194,'Ghir_D01G000260-Ghir_D01G000230*+'],
            ['Ghir_A01',254846,255846,'QTL2*+','Ghir_D01',243712,268247,'Ghir_D01G000410-Ghir_D01G000380*+']
        ]
    ) 
    score,out=run_QTL_lastz(QTLregion.values,genomeObject) 
    print(out,score)  



if __name__ == "__main__":
    genomeFile = '/data/cotton/MaojunWang/WMJ_fiberFullPopulationRNAseq/MappingFPKM/Ghir_Genome_Index/Ghirsutum_genome.fasta'
    genomeObject = pysam.FastaFile(genomeFile)
    QTLregion=pd.read_csv(sys.argv[1],header=None,index_col=None,sep="\t") 
    Matchscore,out=run_QTL_lastz(QTLregion.values,genomeObject)   
    #! QTL match coordinate 
    with open(sys.argv[2],'w') as File:
        for line in out:
            File.write(line)
    #! QTL region match score
    with open(sys.argv[3],'w') as File:
        for line in Matchscore:
            File.write(line)


