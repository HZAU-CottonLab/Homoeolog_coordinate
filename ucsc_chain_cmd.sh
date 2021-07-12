###
# @Descripttion:
# @version:
# @Author: zpliu
# @Date: 2021-07-09 17:22:17
 # @LastEditors: zpliu
 # @LastEditTime: 2021-07-12 10:05:44
# @@param:

#todo :1.提取基因组中每条染色体的序列信息
for i in 01 02 03 04 05 06 07 08 09 10 11 12 13; do
    samtools faidx /public/home/jqyou/data/genome_Garb.HAU/Lachesis_assembly_changed.fa Chr${i} >A2/Chr${i}.fa
    samtools faidx /public/home/jqyou/data/genome_Grai.HAU/Lachesis_assembly_changed.fa Chr${i} >D5/Chr${i}.fa
done

#todo: 2.排序比对过程中重复序列对序列比对的影响
for i in $(ls D5); do
    bsub -q normal -n 1 -e test.err -o test.out -J ${i} -R span[hosts=1] "
    module load  ucsc_kentUtils/v389;
    module load  TRF/4.0.9;
    trfBig D5/${i} trf/D5_${i} -l=4
    "
done
for i in $(ls A2); do
    bsub -q normal -n 1 -e test.err -o test.out -J ${i} -R span[hosts=1] "
    module load  ucsc_kentUtils/v389;
    module load  TRF/4.0.9;
    trfBig A2/${i} trf/A2_${i} -l=4
    "
done

#todo: 3.将输出文件转换为*nib格式
#* */trf/
for i in 01 02 03 04 05 06 07 08 09 10 11 12 13; do
    bsub -q normal -n 1 -e test.err -o test.out -J ${i} -R span[hosts=1] "
    module load BLAT/3.5
    module load  ucsc_kentUtils/v389
    faToNib trf/A2_Chr${i}.fa  A2/Chr${i}.nib 
    faToNib trf/D5_Chr${i}.fa  D5/Chr${i}.nib 
"
done

#todo: 4.进行染色体水平的比对

for i in $(ls D5/ | grep nib); do
    for j in $(ls /data/cotton/zhenpingliu/LZP_fiberFullPopulationRNAseq/04homologExpressionPattern/sub-genome-chain/Blastz/Dt/ | grep nib); do
        out=$(echo ${i}-${j} | sed 's/.nib//g')
        bsub -q normal -M 40G -n 1 -R span[hosts=1] -J ${out} -e test.err -o test.out "
        /data/cotton/zhenpingliu/LZP_fiberFullPopulationRNAseq/software/lastz-1.04.03/lastz-distrib/bin/lastz  D5/${i}  /data/cotton/zhenpingliu/LZP_fiberFullPopulationRNAseq/04homologExpressionPattern/sub-genome-chain/Blastz/Dt/${j}  ‑‑notransition --step=20 --strand=both --format=lav --allocate:traceback=200M  --hspthresh=8000  --nogapped >D5_Dt_lastz/${out}.lav "
    done
done

#todo: 5.将lav 格式转化为psl格式
#* 1. 运行脚本必须在 */D5_Dt_lastz
#* 2. 运行脚本必须在 */A2_At_lastz
for i in *.lav; do
    bsub -q normal -M 40G -n 1 -R span[hosts=1] -J ${out} -e test.err -o test.out "
    module load  ucsc_kentUtils/v389
    lavToPsl $i $(basename $i .lav).psl 
"
done

##todo: 6.将psl转换为chain文件
#* 1. At nib文件路径 /data/cotton/zhenpingliu/LZP_fiberFullPopulationRNAseq/04homologExpressionPattern/sub-genome-chain/Blastz/At
#* 2. Dt nib文件路径 /data/cotton/zhenpingliu/LZP_fiberFullPopulationRNAseq/04homologExpressionPattern/sub-genome-chain/Blastz/Dt

for i in $(ls A2_At_lastz | grep psl); do
    bsub -q normal -n 1 -e test.err -o test.out -R span[hosts=1] "
module load  ucsc_kentUtils/v389
axtChain A2_At_lastz/${i} A2 /data/cotton/zhenpingliu/LZP_fiberFullPopulationRNAseq/04homologExpressionPattern/sub-genome-chain/Blastz/At chain/$(basename $i .psl).chain  -linearGap=medium -psl
"
done
for i in $(ls D5_Dt_lastz | grep psl); do
    bsub -q normal -n 1 -e test.err -o test.out -R span[hosts=1] "
module load  ucsc_kentUtils/v389
axtChain D5_Dt_lastz/${i} D5 /data/cotton/zhenpingliu/LZP_fiberFullPopulationRNAseq/04homologExpressionPattern/sub-genome-chain/Blastz/Dt chain/$(basename $i .psl).chain  -linearGap=medium -psl
"
done

##todo: 7.合并chain同时过滤一些chain
bsub -q normal -n 1 -e test.err -o test.out -R span[hosts=1] "
module load  ucsc_kentUtils/v389
chainMergeSort chain/*Ghir_A*.chain >A2_At.chain
chainPreNet A2_At.chain A2/A2.sizes   /data/cotton/zhenpingliu/LZP_fiberFullPopulationRNAseq/04homologExpressionPattern/sub-genome-chain/Blastz/At/At.sizes  A2_At.pre.chain "


##todo: 8.构建互惠Reciprocal 的最佳比对
## todo: 根据合并后的chain构建net
Atsize=/data/cotton/zhenpingliu/LZP_fiberFullPopulationRNAseq/04homologExpressionPattern/sub-genome-chain/Blastz/At/At.sizes
Dtsize=/data/cotton/zhenpingliu/LZP_fiberFullPopulationRNAseq/04homologExpressionPattern/sub-genome-chain/Blastz/Dt/Dt.sizes
bsub -q high -n 1 -M 300G -e test2.err -o test2.out -R span[hosts=1] "
module load  ucsc_kentUtils/v389 
chainNet -minSpace=1 -minScore=0  A2_At.pre.chain  A2/A2.sizes ${Atsize}  stdout /dev/null  |   netSyntenic stdin  stdout|gzip -c >A2_At.net.gz
"
##todo: 提取net中包含的chain，并且合并小的片段作为最近chain
bsub -q high -n 1 -M 300G -e test3.err -o test3.out -R span[hosts=1] "
module load  ucsc_kentUtils/v389 
#! chainStitchId这一步可能会导致内存不够用
netChainSubset A2_At.net.gz A2_At.pre.chain stdout|chainStitchId stdin stdout |gzip -c A2_At.best.chain.gz 
#! 不对chain进行合并，需要去除注释信息
# netChainSubset A2_At.net.gz A2_At.pre.chain stdout|chainSort stdin stdout|grep \"^#\" -v |gzip -c A2_At.best.chain.gz 
"
##todo: 反方向获取对应的最佳chain
bsub -q high -n 1 -M 300G -e test2.err -o test2.out -R span[hosts=1] "
module load  ucsc_kentUtils/v389 
chainSwap A2_At.best.chain.gz stdout | chainSort stdin stdout|gizp -c > At_A2.chain.gz
"

##todo: 根据反向的chain文件提取net文件
bsub -q high -n 1 -M 300G -e test2.err -o test2.out -R span[hosts=1] "
module load  ucsc_kentUtils/v389 
chainNet -minSpace=1 -minScore=0  At_A2.chain  ${Atsize} A2/A2.sizes stdout /dev/null | netSyntenic stdin At_A2.net
"



#! 批量提交任务################A2 At
A=A2
B=At
preChain=A2_At.pre.chain 
Asize=A2/A2.sizes
Bsize=/data/cotton/zhenpingliu/LZP_fiberFullPopulationRNAseq/04homologExpressionPattern/sub-genome-chain/Blastz/At/At.sizes  
bsub -q high -n 1 -M 300G -e test3.err -o test3.out -R span[hosts=1] "
module load  ucsc_kentUtils/v389 
chainNet -minSpace=1 -minScore=0  ${preChain}  ${Asize} ${Bsize}  stdout /dev/null  |   netSyntenic stdin  stdout|gzip -c >${A}_${B}.net.gz
netChainSubset ${A}_${B}.net.gz ${preChain} stdout|chainSort stdin stdout|grep \"^#\" -v|gzip -c >${A}_${B}.best.chain.gz 
chainSwap ${A}_${B}.best.chain.gz stdout | chainSort stdin stdout|gzip -c > ${B}_${A}.chain.gz 
chainNet -minSpace=1 -minScore=0  ${B}_${A}.chain.gz  ${Bsize} ${Asize} stdout /dev/null  | netSyntenic stdin   stdout |gzip -c >${B}_${A}.net.gz
"

#! ####################D5 vs Dt ###
A=D5
B=Dt
preChain=D5_Dt.pre.chain 
Asize=D5/D5.sizes
Bsize=/data/cotton/zhenpingliu/LZP_fiberFullPopulationRNAseq/04homologExpressionPattern/sub-genome-chain/Blastz/Dt/Dt.sizes  
bsub -q high -n 1 -M 300G -e test3.err -o test3.out -R span[hosts=1] "
module load  ucsc_kentUtils/v389 
chainNet -minSpace=1 -minScore=0  ${preChain}  ${Asize} ${Bsize}  stdout /dev/null  |   netSyntenic stdin  stdout|gzip -c >${A}_${B}.net.gz
netChainSubset ${A}_${B}.net.gz ${preChain} stdout|chainSort stdin stdout|grep \"^#\" -v|gzip -c >${A}_${B}.best.chain.gz 
chainSwap ${A}_${B}.best.chain.gz stdout | chainSort stdin stdout|gzip -c > ${B}_${A}.chain.gz 
chainNet -minSpace=1 -minScore=0  ${B}_${A}.chain.gz  ${Bsize} ${Asize} stdout /dev/null  | netSyntenic stdin   stdout |gzip -c >${B}_${A}.net.gz
"

A=At 
B=Dt
preChain=At_Dt.pre.chain 
Asize=/data/cotton/zhenpingliu/LZP_fiberFullPopulationRNAseq/04homologExpressionPattern/sub-genome-chain/Blastz/At/At.sizes 
Bsize=/data/cotton/zhenpingliu/LZP_fiberFullPopulationRNAseq/04homologExpressionPattern/sub-genome-chain/Blastz/Dt/Dt.sizes 

##########################################################
#todo :提取保守片段
#* 1. 使用bnMApper构建索引
#* 2. 将net文件转换为bed文件
##########################################################
#! 将染色体按照request分开，按染色体进行检索提高检索速度
chainSplit splitFile At_Dt.best.chain.gz -lump=13  
#! 转化为bed文件
netToBed -maxGap=1 At_Dt.net.gz At_Dt.net.bed 
#* 按照单个染色体构建索引文件
bnMapper test.bed 001.chain -k 
