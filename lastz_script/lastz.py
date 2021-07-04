'''
Descripttion: 
version: 
Author: zpliu
Date: 2021-07-04 09:56:04
LastEditors: zpliu
LastEditTime: 2021-07-04 10:09:02
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
def run_lastz_result(geneRegionPairList,homoeologGeneCoordinate:list,genomeObject):
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
        matchRegion: the match region between homoeolog genes
            example 
    '''
    out=''
    for pairedRegion in geneRegionPairList:
        ##
        AtFile=NamedTemporaryFile(mode='w+t',encoding='utf-8')
        DtFile=NamedTemporaryFile(mode='w+t',encoding='utf-8')
        if pairedRegion[1]>=pairedRegion[2] or pairedRegion[4]>=pairedRegion[5]:
            #! the Interval only 1bp
            continue
        AtFile.write(">At\n"+genomeObject.fetch(start=pairedRegion[1],end=pairedRegion[2],region=pairedRegion[0]))
        DtFile.write(">Dt\n"+genomeObject.fetch(start=pairedRegion[4],end=pairedRegion[5],region=pairedRegion[3]))
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
        cigraFlag=os.popen("/data/cotton/zhenpingliu/LZP_fiberFullPopulationRNAseq/software/lastz-1.04.03/lastz-distrib/bin/lastz --strand=both --step=20  --hspthresh=5000 --format=cigar {} {}".format(AtFile.name,DtFile.name)).read()
        out+=parse_lastz_cigra(pairedRegion,cigraFlag,homoeologGeneCoordinate)+"\n"
        # return result
    return out.strip("\n") 
        



def parse_lastz_cigra(pairedRegion,cigraFlagstr,homoeologGeneCoordinate):
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
    '''
    #At gene
    AtgeneId=homoeologGeneCoordinate.iloc[0,3]
    Atstand=homoeologGeneCoordinate.iloc[0,4]
    #Dt gene
    DtgeneId=homoeologGeneCoordinate.iloc[0,8]
    Dtstand=homoeologGeneCoordinate.iloc[0,9]
    out=[]
    if cigraFlagstr:
        #paired region with multiple alginment flagment
        cigraList=cigraFlagstr.strip("\n").split("\n")
        for flag in cigraList:
            #! conserved flag to coordinate
            Dtflatstart,Dtflagend,Dtflagstand,tmp,Atflagstart,Atflagend,Atflagstand,score=flag.split(" ")[2:10]
            Dtstart=min(int(Dtflatstart),int(Dtflagend))+pairedRegion[4]
            Dtend=max(int(Dtflatstart),int(Dtflagend))+pairedRegion[4]
            Atstart=min(int(Atflagstart),int(Atflagend))+pairedRegion[1]
            Atend=max(int(Atflagstart),int(Atflagend))+pairedRegion[1]
            if (Dtstand==Atstand and Dtflagstand==Atflagstand)or(Dtstand!=Atstand and Dtflagstand!=Atflagstand):
                # flag with same order 
                out.append("\t".join([
                    pairedRegion[0],str(Atstart),str(Atend),
                    pairedRegion[3],str(Dtstart),str(Dtend),
                    AtgeneId+"*"+Atstand, DtgeneId+"*"+Dtstand,
                    score,'sameDirection'
                ]))
            else:
                # reverse aligment order 
                out.append("\t".join([
                    pairedRegion[0],str(Atstart),str(Atend),
                    pairedRegion[3],str(Dtstart),str(Dtend),
                    AtgeneId+"*"+Atstand, DtgeneId+"*"+Dtstand,
                    score,'reverseDirection'
                ]))
    else:
        #sequence so diversed and result in low HSP score
         score='0'
         out.append("\t".join([
                    pairedRegion[0],str(pairedRegion[1]),str(pairedRegion[2]),
                    pairedRegion[3],str(pairedRegion[4]),str(pairedRegion[5]),
                    AtgeneId+"*"+Atstand, DtgeneId+"*"+Dtstand,
                    score,'Diverse'
            ]))
    return "\n".join(out)

if __name__=="__main__":
#! test run_lastz_result
    homoeologGeneInfo=pd.DataFrame([['Ghir_A01',68889,76180,'Ghir_A01G000040','-','Ghir_D01',38781,46323,'Ghir_D01G000060','+']])
    genomeFile='/data/cotton/MaojunWang/WMJ_fiberFullPopulationRNAseq/MappingFPKM/Ghir_Genome_Index/Ghirsutum_genome.fasta'
    genomeObject=pysam.FastaFile(genomeFile)
    result=run_lastz_result(
            [['Ghir_A01', 70606, 76180, 'Ghir_D01', 39911, 45620],
            ['Ghir_A01', 68889, 70605, 'Ghir_D01', 38194, 39910],
            ['Ghir_A01', 76181, 76883, 'Ghir_D01', 45621, 46323]],homoeologGeneInfo,genomeObject)
    print(result)


