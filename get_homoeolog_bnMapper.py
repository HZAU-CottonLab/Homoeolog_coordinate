'''
Descripttion:
version:
Author: zpliu
Date: 2021-07-01 23:09:44
LastEditors: zpliu
LastEditTime: 2021-07-02 16:28:59
@param:
'''
from getCoordinate.bnMapper import index_chainFile
from getCoordinate.bnMapper import getTargetLocation
import logging
import pandas as pd
import numpy as np
import pybedtools
import argparse
logger = logging.getLogger()
logger.setLevel(logging.INFO)
logger.info("loading module...")
###########################
# get request bin for chain
###########################


def getBedtoolsObject(interval_list):
    '''
    # get Bedtools form list
    args:
        interval_list - one line bed
    return: BedTools object
    '''
    s = "\t".join([str(j) for j in interval_list])
    return pybedtools.BedTool(s, from_string=True)


def getRequestBinFromIntersect(intersectList):
    '''
    get chain request Interval by intersect with gene region
    args:
        BedToolsObject: intersect outData
    return:
        request chain: @list
    '''
    out = []
    for item in intersectList:
        try:
            start1 = int(item[1])
            end1 = int(item[2])
            start2 = int(item[6])
            end2 = int(item[7])
        except ValueError:
            #! no intersect with request chain
            continue
        if start1 <= start2 and end1 >= end2:
            out.append((item[0], start2, end2, item[3]+"*"+item[4]))
        elif start1 < start2 and end1 < end2:
            out.append((item[0], start2, end1, item[3]+"*"+item[4]))
        elif start1 > start2 and end1 > end2:
            out.append((item[0], start1, end2, item[3]+"*"+item[4]))
    return out


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
    ###############################
    #! only need this thwo arguments
    ###############################
    parser.add_argument("-homoeolog", default="BED",)
    parser.add_argument("-myOut", default="BED",)
    opt = parser.parse_args()
    homoleologBed = pd.read_csv(
        opt.homoeolog, header=None, index_col=None, sep="\t")
    At = homoleologBed[[0, 1, 2, 3, 4]].values
    Dt = homoleologBed[[5, 6, 7, 8, 9]].values

    ##############################
    # read request chain
    ##############################
    At_request_net = pybedtools.BedTool("../Blastz/out2/At_Dt.rbest.net.bed")
    Dt_request_net = pybedtools.BedTool("../Blastz/out2/Dt_At.rbest.net.bed")

    # bnMapper API
    # index the chain file for getting faster
    ###########################################
    AtEPO, AtTree = index_chainFile('../Blastz/out2/At_Dt.rbest.chain.gz', opt)
    DtEPO, DtTree = index_chainFile('../Blastz/out2/Dt_At.rbest.chain.gz', opt)

    out = []  # log of gene pairs
    for Atitem, Dtitem in zip(At, Dt):
        # iteral homoeolog gene pairs
        At_interactive = getBedtoolsObject(
            Atitem).intersect(At_request_net, loj=True)
        Dt_interactive = getBedtoolsObject(
            Dtitem).intersect(Dt_request_net, loj=True)
        # print(At_interactive)
        # print(Dt_interactive)
        ##############################
        # get the target location
        ##############################
        At_request = getRequestBinFromIntersect(At_interactive)
        Dt_request = getRequestBinFromIntersect(Dt_interactive)

        for item in At_request:
            tmpArray = getTargetLocation(item, AtEPO, AtTree, opt)
            # print(tmpArray)
            if(tmpArray):
                #! request matp to more than one target region
                out.append("\t".join([item[0], str(item[1]), str(item[2]), item[3]]) +
                        "\t"+"\t".join([i[0]+":"+str(i[1])+"-"+str(i[2]) for i in tmpArray]))
            else:
                #! no target region mapper
                out.append(
                    "\t".join([item[0], str(item[1]), str(item[2]), item[3]])+"\t"+"None")
        ########################################
        #! Dt request  to get At target mapping
        ########################################
        for item in Dt_request:
            tmpArray = getTargetLocation(item, DtEPO, DtTree, opt)
            # print(tmpArray)
            if(tmpArray):
                #! request matp to more than one target region
                out.append("\t".join(
                    [item[0], str(item[1]), str(item[2]), item[3]])+"\t"+"\t".join([i[0]+":"+str(i[1])+"-"+str(i[2]) for i in tmpArray]))
            else:
                #! no target region mapper
                out.append(
                    "\t".join([item[0], str(item[1]), str(item[2]), item[3]])+"\t"+"None")
    with open(opt.myOut, 'w') as File:
        for line in out:
            File.write(line+"\n")

