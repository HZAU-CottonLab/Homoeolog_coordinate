'''
Descripttion: 
version: 
Author: zpliu
Date: 2021-07-04 09:56:04
LastEditors: zpliu
LastEditTime: 2021-07-06 15:36:37
@param: 
'''
from tempfile import NamedTemporaryFile
import pysam
import os
import pandas as pd
###################################################
#! run lastz
###################################################
# from tempfile import NamedTemporaryFile
# import os
# import redis


def run_lastz_result(geneRegionPairList, homoeologGeneCoordinate: list, genomeObject):
    '''
    run lastz to get match region
    args:
        -geneRegionPairList: the paired region @list
            example:[['Ghir_A01', 70606, 76180, 'Ghir_D01', 39911, 45620],
                     ['Ghir_A01', 68889, 70605, 'Ghir_D01', 38194, 39910],
                     ['Ghir_A01', 76181, 76883, 'Ghir_D01', 45621, 46323]]
        -homoeologGeneCoordinate: pd.DataFrame
            example:  [	Ghir_A01,68889,76180,Ghir_A01G000040,-,Ghir_D01,38781,46323,Ghir_D01G000060,+] 

        -genomeObject: genome sequence file(pysam)           
    return:
        -matchRegion: the match region between homoeolog genes
            example
        -align score: all fragement score
    '''
    out = ''
    totalscore=0
    for pairedRegion in geneRegionPairList:
        ##
        AtFile = NamedTemporaryFile(mode='w+t', encoding='utf-8')
        DtFile = NamedTemporaryFile(mode='w+t', encoding='utf-8')
        if pairedRegion[1] >= pairedRegion[2] or pairedRegion[4] >= pairedRegion[5]:
            #! the Interval only 1bp
            continue
        try:
            ###################index start with 0
            AtFile.write(
                ">At\n"+genomeObject.fetch(start=pairedRegion[1]-1, end=pairedRegion[2], region=pairedRegion[0]))
            DtFile.write(
                ">Dt\n"+genomeObject.fetch(start=pairedRegion[4]-1, end=pairedRegion[5], region=pairedRegion[3]))
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
                "/data/cotton/zhenpingliu/LZP_fiberFullPopulationRNAseq/software/lastz-1.04.03/lastz-distrib/bin/lastz --strand=both --step=20  --hspthresh=5000 --format=cigar {} {}".format(AtFile.name, DtFile.name)).read()
        except:
            #! if region out of chromsome the cigra will be ''
            print('the region out of chromsome', end="\t")
            print(pairedRegion)
            continue
        tmpout,tmpscore= parse_lastz_cigra(pairedRegion, cigraFlag,
                                 homoeologGeneCoordinate)
        totalscore+=int(tmpscore)                         
        out+=tmpout
        # return result
    return (out.strip("\n"),str(totalscore))


def parse_lastz_cigra(pairedRegion, cigraFlagstr, homoeologGeneCoordinate):
    '''
    parse the lastz out to get conserved region
    args:
        - pairedRegion: ther paired region of homoeolog gene
            example@list ['Ghir_A01', 70606, 76180, 'Ghir_D01', 39911, 45620]
        - cigraFlagstr: @str 'cigar: Dt 0 703 + At 0 702 + 60573 M 85 D 1 M 124 I 2 M 492\n'
        -homoeologGeneCoordinate: pd.DataFrame
            example:  [	Ghir_A01,68889,76180,Ghir_A01G000040,-,Ghir_D01,38781,46323,Ghir_D01G000060,+]   
    return:
        - matchRegion Ghir_A01 75171 76180 Ghir_D01 44594 45620 Ghir_A01G000040*- Ghir_D01G000060*+ 82608 reverseDirection
        - total score value 
    '''
    # print(pairedRegion)
    # print(cigraFlagstr)
    # At gene
    AtgeneId = homoeologGeneCoordinate.iloc[0, 3]
    Atstand = homoeologGeneCoordinate.iloc[0, 4]
    # Dt gene
    DtgeneId = homoeologGeneCoordinate.iloc[0, 8]
    Dtstand = homoeologGeneCoordinate.iloc[0, 9]
    out = []
    totalscore=0
    print(cigraFlagstr)
    if cigraFlagstr:
        # paired region with multiple alginment flagment
        cigraList = cigraFlagstr.strip("\n").split("\n")
        for flag in cigraList:
            #! conserved flag to coordinate
            Dtflatstart, Dtflagend, Dtflagstand, tmp, Atflagstart, Atflagend, Atflagstand, score = flag.split(" ")[
                2:10]
            totalscore+=int(score)
            Dtstart = int(Dtflatstart)+pairedRegion[4]
            # # Dtend = max(int(Dtflatstart), int(Dtflagend))+pairedRegion[4]
            Atstart = int(Atflagstart)+pairedRegion[1]
            # Atend = max(int(Atflagstart), int(Atflagend))+pairedRegion[1]
            #! :: avoid the last one
            cigarFlag = flag.split(" ")[10::2]
            cigarBase = flag.split(" ")[11::2]
            ###coordinate begain at 1 
            # cigarBase[0]=int(cigarBase[0])-1
            cigarBase=[int(i)-1 for i in cigarBase]
            # print(cigarBase)
            # print(cigarFlag)
            # print(cigarBase)
            if (Dtstand == Atstand and Dtflagstand == Atflagstand) or (Dtstand != Atstand and Dtflagstand != Atflagstand):
                # flag with same order
                # out.append("\t".join([
                #     pairedRegion[0],str(Atstart),str(Atend),
                #     pairedRegion[3],str(Dtstart),str(Dtend),
                #     AtgeneId+"*"+Atstand, DtgeneId+"*"+Dtstand,
                #     score,'sameDirection'
                # ]))
                # Accurate to each base
                if Dtflagstand == "+" and Atflagstand == "+":
                    inputData = [pairedRegion[0], pairedRegion[3], Atstart, Dtstart,
                                 AtgeneId, Atstand, DtgeneId, Dtstand,
                                 1, 1, 'sameDirection']
                    out+=calcuate_coordinate(cigarFlag, cigarBase, inputData)
                elif Dtflagstand == "+" and Atflagstand == "-":
                    ##! At with reverse order
                    inputData = [pairedRegion[0], pairedRegion[3], Atstart, Dtstart,
                                 AtgeneId, Atstand, DtgeneId, Dtstand,
                                 1, -1, 'sameDirection']
                    out+=calcuate_coordinate(cigarFlag, cigarBase, inputData)
                    
                elif Dtflagstand == "-" and Atflagstand == "+":
                    ##! Dt with reverse order 
                    inputData = [pairedRegion[0], pairedRegion[3], Atstart, Dtstart,
                                 AtgeneId, Atstand, DtgeneId, Dtstand,
                                 -1, 1, 'sameDirection']
                    out+=calcuate_coordinate(cigarFlag, cigarBase, inputData)
                else:
                    print("unfortunately",)
            ############################################################
            # with reverse Direction
            ############################################################
            else:
                if Dtflagstand == "+" and Atflagstand == "+":
                    inputData = [pairedRegion[0], pairedRegion[3], Atstart, Dtstart,
                                 AtgeneId, Atstand, DtgeneId, Dtstand,
                                 1, 1, 'reverseDirection']
                    out+=calcuate_coordinate(cigarFlag, cigarBase, inputData)
                elif Dtflagstand == "+" and Atflagstand == "-":
                    ##! At with reverse order
                    inputData = [pairedRegion[0], pairedRegion[3], Atstart, Dtstart,
                                 AtgeneId, Atstand, DtgeneId, Dtstand,
                                 1, -1, 'reverseDirection']
                    out+=calcuate_coordinate(cigarFlag, cigarBase, inputData)
                    
                elif Dtflagstand == "-" and Atflagstand == "+":
                    ##! Dt with reverse order 
                    inputData = [pairedRegion[0], pairedRegion[3], Atstart, Dtstart,
                                 AtgeneId, Atstand, DtgeneId, Dtstand,
                                 -1, 1, 'reverseDirection']
                    out+=calcuate_coordinate(cigarFlag, cigarBase, inputData)
                else:
                    print("unfortunately",)
                # reverse aligment order
                # out.append("\t".join([
                #     pairedRegion[0], str(Atstart), str(Atend),
                #     pairedRegion[3], str(Dtstart), str(Dtend),
                #     AtgeneId+"*"+Atstand, DtgeneId+"*"+Dtstand,
                #     score, 'reverseDirection'
                # ]))
    else:
        # sequence so diversed and result in low HSP score
        score = '0'
        out.append("\t".join([
                   pairedRegion[0], str(pairedRegion[1]), str(pairedRegion[2]),
                   pairedRegion[3], str(pairedRegion[4]), str(pairedRegion[5]),
                   AtgeneId+"*"+Atstand, DtgeneId+"*"+Dtstand,'Diverse','-'
                   ]))
    # return "\n".join(out),str(totalscore)

    return  "\n".join(out),str(totalscore)


def calcuate_coordinate(cigarFlag, cigarBase, inputData):
    '''
    ## when At and Dt with reverse order; use the different multiple
    '''
    out = []
    chromAt, chromDt, Atstart, Dtstart, AtgeneId, Atstand, DtgeneId, Dtstand, Dtmultiple, Atmultiple, Direction = inputData
    for matchflag, matchlength in zip(cigarFlag, cigarBase):
        if matchflag == "M":
            ###############
            # add same time
            ###############
            At=[Atstart,Atstart+matchlength*Atmultiple]
            Dt=[Dtstart,Dtstart+matchlength*Dtmultiple]
            At.sort() 
            Dt.sort() 
            out.append("\t".join([
                chromAt, str(At[0]), str(At[1]),
                chromDt, str(Dt[0]), str(Dt[1]),
                AtgeneId+"*"+Atstand, DtgeneId+"*"+Dtstand, Direction, 'Match'
            ]))
            Atstart += (matchlength+1)*Atmultiple
            Dtstart += (matchlength+1)*Dtmultiple
        elif matchflag == "I":
            ##################
            # Dt special base
            ##################
            print(Dtstart)
            Dt=[Dtstart,Dtstart+matchlength*Dtmultiple]
            Dt.sort()
            out.append("\t".join([
                chromAt, str(Atstart), str(Atstart+matchlength*0),
                chromDt, str(Dt[0]), str(Dt[1]),
                AtgeneId+"*"+Atstand, DtgeneId+"*"+Dtstand, Direction, 'Insert'
            ]))
            Atstart += 0
            Dtstart += (matchlength+1)*Dtmultiple
        else:
            #################
            # At special base
            #################
            At=[Atstart,Atstart+matchlength*Atmultiple]
            At.sort()
            out.append("\t".join([
                chromAt, str(At[0]), str(At[1]),
                chromDt, str(Dtstart), str(Dtstart+matchlength*0),
                AtgeneId+"*"+Atstand, DtgeneId+"*"+Dtstand, Direction, 'Delte'
            ]))
            Atstart += (matchlength+1)*Atmultiple
            Dtstart += 0
    return out


if __name__ == "__main__":
    #! test run_lastz_result
    homoeologGeneInfo = pd.DataFrame(
        [['Ghir_A01', 771697, 772409, 'Ghir_A01G000070', '+', 'Ghir_D01', 695574, 696287, 'Ghir_D01G000110', '+']])
    genomeFile = '/data/cotton/MaojunWang/WMJ_fiberFullPopulationRNAseq/MappingFPKM/Ghir_Genome_Index/Ghirsutum_genome.fasta'
    genomeObject = pysam.FastaFile(genomeFile)
    result ,score= run_lastz_result(
        [['Ghir_A01', 771697, 772409, 'Ghir_D01', 695574, 696287],], homoeologGeneInfo, genomeObject)
    print(result,score)
# [
#     ['Ghir_A01', 70606, 76180, 'Ghir_D01', 39911, 45620],
#     ['Ghir_A01', 68889, 70605, 'Ghir_D01', 38194, 39910],
#     ['Ghir_A01', 76181, 76883, 'Ghir_D01', 45621, 46323]
# ]
