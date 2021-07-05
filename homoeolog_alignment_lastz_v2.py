'''
Descripttion: 
version: 
Author: zpliu
Date: 2021-07-05 17:04:33
LastEditors: zpliu
LastEditTime: 2021-07-05 20:12:16
@param: 
'''
import logging
import pandas as pd
import numpy as np
import pybedtools
import re
import pysam
from lastz_script.lastz import run_lastz_result
import sys
'''
#! test data form 
Ghir_A01	71162	72994	Ghir_A01G000040*-	Ghir_D01	42201	43973	Ghir_D01G000060*+	codeRegion
Ghir_A01	96241	99719	Ghir_A01G000070*+	Ghir_D01	79012	82604	Ghir_D01G000110*+	codeRegion
Ghir_A01	68889	71161	Ghir_A01G000040*-	Ghir_D01	43974	46323	Ghir_D01G000060*+	downRegion
Ghir_A01	99720	102744	Ghir_A01G000070*+	Ghir_D01	82605	85636	Ghir_D01G000110*+	downRegion
Ghir_A01	72995	76180	Ghir_A01G000040*-	Ghir_D01	38781	42200	Ghir_D01G000060*+	PromoterRegion
'''
##############################
# extract genome sequence
##############################
genomeFile = '/data/cotton/MaojunWang/WMJ_fiberFullPopulationRNAseq/MappingFPKM/Ghir_Genome_Index/Ghirsutum_genome.fasta'
genomeObject = pysam.FastaFile(genomeFile)
alignRegion = pd.read_csv(sys.argv[1], header=None, index_col=None, sep="\t")

matchResult = []  # 碱基水平的比对结果
matchRegion = []  # 两端是否会发生延伸 
##延伸后scode相差500bp以上
def f(x, y): return True if int(x)-int(y) >= 10000 else False
for pairdRegion in alignRegion.values:
    # 对每对同源基因进行比较
    AtgeneId, Atstand = pairdRegion[3].split("*")
    DtgeneId, Dtstand = pairdRegion[7].split("*")
    Regiontype = pairdRegion[-1]
    homoeologGeneInfo = pd.DataFrame(
        [list(pairdRegion[0:3])+pairdRegion[3].split("*")+list(pairdRegion[4:7])+pairdRegion[7].split("*")])
    alignmentRegion = np.append(
        pairdRegion[0:3],pairdRegion[4:7]
    )
    if Regiontype == 'codeRegion':
        alignResult, score = run_lastz_result(
            [alignmentRegion], homoeologGeneInfo, genomeObject)
        matchResult.append(alignResult)
        alignmentRegion=np.append(alignmentRegion,[pairdRegion[3], pairdRegion[7], Regiontype, score])
        matchRegion.append(alignmentRegion)
    if Regiontype == "PromoterRegion":
        ##################################
        # 进行三次比对，筛选其中得分最高的比对
        # 1. 正常的不低
        # 2. At不变，Dt往右边挪动500bp
        # 3. Dt不变，At往右边挪动500bp
        ##################################
        # 1.正常情况下
        alignResult1, score1 = run_lastz_result(
            [alignmentRegion], homoeologGeneInfo, genomeObject)
        # 2.At挪动位置，保证不会得到负数
        alignmentRegion2 = alignmentRegion
        homoeologGeneInfo2 = homoeologGeneInfo
        if Atstand == "+":
            alignmentRegion2[1] = alignmentRegion2[1]-500
            if alignmentRegion2[1] <= 0:
                # 防止减出负数
                alignmentRegion2[1] = 1
            else:
                pass
            homoeologGeneInfo2.iloc[0, 1] = alignmentRegion2[1]
        else:
            alignmentRegion2[2] = alignmentRegion2[2]+500
            homoeologGeneInfo2.iloc[0, 2] = alignmentRegion2[2]
        alignResult2, score2 = run_lastz_result(
            [alignmentRegion2], homoeologGeneInfo2, genomeObject)
        ####################################
        # Dt 挪动位置
        ####################################
        alignmentRegion3 = alignmentRegion
        homoeologGeneInfo3 = homoeologGeneInfo
        if Dtstand == "+":
            alignmentRegion3[4] = alignmentRegion3[4]-500
            if alignmentRegion3[4] <= 0:
                # 防止减出负数
                alignmentRegion3[4] = 1
            else:
                pass
            homoeologGeneInfo3.iloc[0, 6] = alignmentRegion3[4]
        else:
            alignmentRegion3[5] = alignmentRegion3[5]+500
            homoeologGeneInfo3.iloc[0, 7] = alignmentRegion3[5]
        alignResult3, score3 = run_lastz_result(
            [alignmentRegion3], homoeologGeneInfo3, genomeObject)
        ######################################
        # score 相差超过10000得分
        
        ######################################
        if not f(score3, score1) and not f(score2, score1):
            matchResult.append(alignResult1)
            alignmentRegion =np.append(alignmentRegion,[pairdRegion[3],
                                pairdRegion[7], Regiontype, score1])
            matchRegion.append(alignmentRegion)
        elif f(score3, score1) and not f(score2, score1):
            matchResult.append(alignResult3)
            alignmentRegion3 =np.append(alignmentRegion3, [pairdRegion[3],
                                 pairdRegion[7], Regiontype, score3])
            matchRegion.append(alignmentRegion3)
        elif not f(score3, score1) and f(score2, score1):
            matchResult.append(alignResult2)
            alignmentRegion2 =np.append(alignmentRegion2, [pairdRegion[3],
                                 pairdRegion[7], Regiontype, score2])
            matchRegion.append(alignmentRegion2)
        else:
            # 两边延伸效果都可以
            # 不延伸
            matchResult.append(alignResult1)
            alignmentRegion =np.append(alignmentRegion, [pairdRegion[3],
                                pairdRegion[7], Regiontype, score1])
            matchRegion.append(alignmentRegion)
    ###################################################################
    # downRegion
    ###################################################################
    if Regiontype == "downRegion":
        ##################################
        # 进行三次比对，筛选其中得分最高的比对
        # 1. 正常的不低
        # 2. At不变，Dt往右边挪动500bp
        # 3. Dt不变，At往右边挪动500bp
        ##################################
        # 1.正常情况下
        alignResult1, score1 = run_lastz_result(
            [alignmentRegion], homoeologGeneInfo, genomeObject)
        # 2.At挪动位置，保证不会得到负数
        alignmentRegion2 = alignmentRegion
        homoeologGeneInfo2 = homoeologGeneInfo
        if Atstand == "-":
            alignmentRegion2[1] = alignmentRegion2[1]-500
            if alignmentRegion2[1] <= 0:
                # 防止减出负数
                alignmentRegion2[1] = 1
            else:
                pass
            homoeologGeneInfo2.iloc[0, 1] = alignmentRegion2[1]
        else:
            alignmentRegion2[2] = alignmentRegion2[2]+500
            homoeologGeneInfo2.iloc[0, 2] = alignmentRegion2[2]
        alignResult2, score2 = run_lastz_result(
            [alignmentRegion2], homoeologGeneInfo2, genomeObject)
        ####################################
        # Dt 挪动位置
        ####################################
        alignmentRegion3 = alignmentRegion
        homoeologGeneInfo3 = homoeologGeneInfo
        if Dtstand == "-":
            alignmentRegion3[4] = alignmentRegion3[4]-500
            if alignmentRegion3[4] <= 0:
                # 防止减出负数
                alignmentRegion3[4] = 1
            else:
                pass
            homoeologGeneInfo3.iloc[0, 6] = alignmentRegion3[4]
        else:
            alignmentRegion3[5] = alignmentRegion3[5]+500
            homoeologGeneInfo3.iloc[0, 7] = alignmentRegion3[5]
        alignResult3, score3 = run_lastz_result(
            [alignmentRegion3], homoeologGeneInfo3, genomeObject)
        ######################################
        # score 相差超过10000得分
        ######################################
        if not f(score3, score1) and not f(score2, score1):
            matchResult.append(alignResult1)
            alignmentRegion = np.append(alignmentRegion,[pairdRegion[3],
                                pairdRegion[7], Regiontype, score1])
            matchRegion.append(alignmentRegion)
        elif f(score3, score1) and not f(score2, score1):
            matchResult.append(alignResult3)
            alignmentRegion3 =np.append(alignmentRegion3, [pairdRegion[3],
                                 pairdRegion[7], Regiontype, score3])
            matchRegion.append(alignmentRegion3)
        elif not f(score3, score1) and f(score2, score1):
            matchResult.append(alignResult2)
            alignmentRegion2 =np.append(alignmentRegion2, [pairdRegion[3],
                                 pairdRegion[7], Regiontype, score2])
            matchRegion.append(alignmentRegion2)
        else:
            # 两边延伸效果都可以
            # 不延伸
            matchResult.append(alignResult1)
            alignmentRegion =np.append(alignmentRegion, [pairdRegion[3],
                                pairdRegion[7], Regiontype, score1])
            matchRegion.append(alignmentRegion)

with open(sys.argv[2], 'w') as File:
    for line in matchResult:
        File.write(line+"\n")
with open(sys.argv[3], 'w') as File:
    for line in matchRegion:
        File.write("\t".join([str(i) for i in line])+"\n")
