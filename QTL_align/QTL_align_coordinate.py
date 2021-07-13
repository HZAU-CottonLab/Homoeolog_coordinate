'''
Descripttion: 
version: 
Author: zpliu
Date: 2021-07-13 16:39:22
LastEditors: zpliu
LastEditTime: 2021-07-13 21:15:19
@param: 
'''
#########################################################
#!根据同源基因的共线性关系，获取QTL需要比对的区间
#########################################################
import pandas as pd
import sys

chromsomeSize = pd.read_csv(
    "/data/cotton/zhenpingliu/jqyou/TM1/TM1_chr_length.txt", header=None, index_col=0, sep="\t")
'''chromsomeSize File
    Ghir_A01        117757855
    Ghir_A02        108092100
    Ghir_A03        113059412
    Ghir_A04        85149810

'''    
#! all gene Coordinate
geneCoordinate = pd.read_csv(
    "/data/cotton/zhenpingliu/LZP_fiberFullPopulationRNAseq/03express_gene/eQTLgenes_Allstages_340samples/fastQTL_peer_20/Gene_information/All_gene_bed.txt", header=None, index_col=None, sep="\t")
'''gene Coordinate
Ghir_A01        80323913        80324566        Ghir_A01G013980 +
Ghir_A01        80322152        80323048        Ghir_A01G013970 +
Ghir_A01        60753987        60754357        Ghir_A01G013100 +
Ghir_A01        59980270        59983471        Ghir_A01G013050 +
Ghir_A01        59140877        59148861        Ghir_A01G013020 +
Ghir_A01        60294264        60295580        Ghir_A01G013060 +
'''
QTL_flank_gene = pd.read_csv(
    sys.argv[1], header=None, index_col=None, sep="\t")
'''QTL_flank_gene
    Ghir_A01        164527  164527  Ghir_A01G000250,Ghir_A01G000220 Ghir_D01G000260,Ghir_D01G000230 Ghir_A01        Ghir_D01        QTL1
    Ghir_A01        255346  255346  Ghir_A01G000410,Ghir_A01G000380 Ghir_D01G000410,Ghir_D01G000380 Ghir_A01        Ghir_D01        QTL2
    Ghir_A01        269159  269159  Ghir_A01G000410,Ghir_A01G000380 Ghir_D01G000410,Ghir_D01G000380 Ghir_A01        Ghir_D01        QTL3
'''
out = []
for value in QTL_flank_gene.values:
    geneAlist = value[3].strip(",").split(",")
    geneBlist = value[4].strip(",").split(",")
    ChromAlist = value[5].strip(",").split(",")
    ChromBlist = value[6].strip(",").split(",")
    QTLsite = value[1]
    QTLId = value[7]+"*"+"+"
    ##################################
    #! flank length of lead SNP
    ##################################
    flanklength = 500
    if QTLsite >= flanklength:
        QTLflagstart = QTLsite-flanklength
    else:
        QTLflagstart = 1
    if QTLsite+flanklength <= chromsomeSize.loc[ChromAlist[0]][1]:
        QTLflagend = QTLsite+flanklength
    else:
        QTLflagend = chromsomeSize.loc[ChromAlist[0]]
    ######################################
    # target region
    ######################################
    if len(geneAlist) == 1:
        # only one adjacentgene
        chromsomeBsize = chromsomeSize.loc[ChromBlist[0]][1]
        chromsomeB, startB, endB = geneCoordinate.loc[geneCoordinate[3]
                                                      == geneBlist[0]].iloc[0, 0:3]
        startB, endB = sorted([startB, endB])
        if startB < chromsomeBsize-startB:
            #!gene on left of chromosome
            tmpData = [ChromAlist[0], QTLflagstart, QTLflagend,
                       QTLId, chromsomeB, 1, startB, geneBlist[0]+"*"+"+"]
            out.append([str(i) for i in tmpData])
        else:
            # gene on right of chromosome
            tmpData = [ChromAlist[0], QTLflagstart, QTLflagend, QTLId,
                       chromsomeB, endB, chromsomeBsize, geneBlist[0]+"*"+"+"]
            out.append([str(i) for i in tmpData])

    else:
        #! two gene with one chromsome or two chromsomes
        if len(ChromBlist) == 1:
            #!adjacent homoeologene in the same chromsome
            chromsomeB1, startB1, endB1 = geneCoordinate.loc[geneCoordinate[3]
                                                             == geneBlist[0]].iloc[0, 0:3]
            chromsomeB2, startB2, endB2 = geneCoordinate.loc[geneCoordinate[3]
                                                             == geneBlist[1]].iloc[0, 0:3]
            start, tmp1, tmp2, end = sorted([startB1, endB1, startB2, endB2])
            tmpData = [ChromAlist[0], QTLflagstart, QTLflagend, QTLId,
                       chromsomeB1, start, end, "-".join(geneBlist)+"*"+"+"]
            out.append([str(i) for i in tmpData])
        else:
            #!adjacent homoeologene in the differenc chromsome
            chromsomeA1, startA1, endA1 = geneCoordinate.loc[geneCoordinate[3]
                                                             == geneAlist[0]].iloc[0, 0:3]
            chromsomeA2, startA2, endA2 = geneCoordinate.loc[geneCoordinate[3]
                                                             == geneAlist[1]].iloc[0, 0:3]
            start, tmp1, tmp2, end = sorted([startA1, endA1, startA2, endA2])
            #!QTL in this interval 
            IntervalLength = end-start+1
            for geneId in geneBlist:
                chromsomeB1, startB1, endB1 = geneCoordinate.loc[geneCoordinate[3]
                                                                 == geneId].iloc[0, 0:3]
                # chromsomeB1length=chromsomeSize.loc[chromsomeB1]
                startB1, endB1 = sorted([startB1, endB1])
                if startB1-IntervalLength <= 0:
                    tmpData = [ChromAlist[0], QTLflagstart, QTLflagend, QTLId,
                               chromsomeB1, 1, endB1+IntervalLength, geneId+"*"+"+"]
                    out.append([str(i) for i in tmpData])
                else:
                    tmpData = [ChromAlist[0], QTLflagstart, QTLflagend, QTLId, chromsomeB1,
                               startB1-IntervalLength, endB1+IntervalLength, geneId+"*"+"+"]
                    out.append([str(i) for i in tmpData])
with open(sys.argv[2], 'w') as File:
    for item in out:
        File.write("\t".join(item)+"\n")
'''QTL align coordinate
    Ghir_A01        164027  165027  QTL1*+  Ghir_D01        145367  162194  Ghir_D01G000260-Ghir_D01G000230*+
    Ghir_A01        254846  255846  QTL2*+  Ghir_D01        243712  268247  Ghir_D01G000410-Ghir_D01G000380*+
    Ghir_A01        268659  269659  QTL3*+  Ghir_D01        243712  268247  Ghir_D01G000410-Ghir_D01G000380*+
    Ghir_A01        276630  277630  QTL4*+  Ghir_D01        243712  268247  Ghir_D01G000410-Ghir_D01G000380*+
'''        
