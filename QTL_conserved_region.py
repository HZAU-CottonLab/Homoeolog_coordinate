'''
Descripttion: 
version: 
Author: zpliu
Date: 2021-07-13 08:56:35
LastEditors: zpliu
LastEditTime: 2021-07-13 11:15:00
@param: 
'''
from bnMapper import index_chainFile
from bnMapper import getTargetLocation
import logging
import pandas as pd
import numpy as np
import pybedtools
import argparse
logger = logging.getLogger()
logger.setLevel(logging.INFO)
logger.info("loading module...")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-input", nargs='+', default="None")
    parser.add_argument("-alignment", default="None")

    parser.add_argument("-f", '--format', default="BED4",)
    parser.add_argument("-o", '--output', default='stdout',)
    parser.add_argument("-t", '--threshold',
                        metavar="FLOAT", default=0., type=float,)
    parser.add_argument("-s", '--screen', default=False,)
    parser.add_argument('-g', '--gap', type=int, default=-1,)
    parser.add_argument('-v', '--verbose', type=str,  default='info',)
    parser.add_argument("-k", '--keep_split',
                        default=False, action='store_true',)
    parser.add_argument("-i", default="BED",)

    #! only need this thwo arguments
    ###############################
    parser.add_argument("-QTL", default="BED",) 
    '''QTL region file 
    + fild1 Ghir_A01	
    + fild2 586316
    + fild3 586316
    + fild4 Ghir_A01G000740-Ghir_D01G000710 (eGene)
    + fild5 IntraChrAt (QTLtype:cis,trans)
    + fild6 0DPA (stage of development)
    + fild7 BiasQTL
    + fild8 Ghir_A01:586033-586241 (upstream region !import)
    + fild9 Ghir_A01:586309-587347 (intersect with QTL region !import )
    + fild10 Ghir_A01:587356-589248 (downstream region !import)
    '''
    parser.add_argument("-myOut", default="BED",)
    opt = parser.parse_args()
    # * chain psl file
    EPO, Tree = index_chainFile(opt.alignment, opt)
    # * search the conserced region
    # seqrceRegionExample = ('Ghir_D12', 55928807, 55928907, 'conserved')
    # tmpArray=getTargetLocation(seqrceRegionExample,EPO,Tree,opt)
    QTLIntersectBed=pd.read_csv(opt.QTL,header=None,index_col=None,sep="\t")
    Upregion,conservedRegion,downregion=zip(*QTLIntersectBed[[7,8,9]].values)
    # print(Upregion)
    # print(conservedRegion)
    # print(downregion)
    upOut=[]
    conservedOut=[]
    downOut=[]
    ########################################
    #! up region
    ########################################
    for item in Upregion:
        if item=="None":
            #! donw have record
            upOut.append('None')
        else:
            #!Ghir_A01:586309-587347
            chromsome=item.split(":")[0]
            start=item.split("-")[0].split(":")[1]
            end=item.split("-")[1]
            seqrceRegionExample=(chromsome,int(start),int(end),'up')
            tmpArray=getTargetLocation(seqrceRegionExample,EPO,Tree,opt)
            tmpstr=''
            if tmpArray:
                for outitem in tmpArray:
                    #* one to many
                    tmpstr+=outitem[0]+":"+str(outitem[1])+"-"+str(outitem[2])+","
                upOut.append(tmpstr)
            else:
                upOut.append("None")
    #######################################
    #! conserved region
    #######################################
    for item in conservedRegion:
        if  item=="None":
            #! donw have record
            conservedOut.append('None')
        else:
            chromsome=item.split(":")[0]
            start=item.split("-")[0].split(":")[1]
            end=item.split("-")[1]
            seqrceRegionExample=(chromsome,int(start),int(end),'up')
            tmpArray=getTargetLocation(seqrceRegionExample,EPO,Tree,opt)
            tmpstr=''
            if(tmpArray):
                for outitem in tmpArray:
                    tmpstr+=outitem[0]+":"+str(outitem[1])+"-"+str(outitem[2])+","
                conservedOut.append(tmpstr)
            else:
                #! without conserved chain
                conservedOut.append("None")
                
    #################################################
    #! down region
    #################################################
    for item in downregion:
        if item=="None":
            #! donw have record
            downOut.append('None')
        else:
            chromsome=item.split(":")[0]
            start=item.split("-")[0].split(":")[1]
            end=item.split("-")[1]
            seqrceRegionExample=(chromsome,int(start),int(end),'up')
            tmpArray=getTargetLocation(seqrceRegionExample,EPO,Tree,opt)
            tmpstr=''
            if tmpArray:
                for outitem in tmpArray:
                    tmpstr+=outitem[0]+":"+str(outitem[1])+"-"+str(outitem[2])+","
                downOut.append(tmpstr)
            else:
                downOut.append("None")
    QTLIntersectBed['upConserved']=upOut
    QTLIntersectBed['Conserved']=conservedOut
    QTLIntersectBed['down']=downOut
    # print(QTLIntersectBed)
    QTLIntersectBed.to_csv(opt.myOut,header=False,index=False,sep="\t")
                

            
    
    

