<!--
 * @Descripttion: 
 * @version: 
 * @Author: zpliu
 * @Date: 2021-07-04 10:19:51
 * @LastEditors: zpliu
 * @LastEditTime: 2021-07-13 21:34:25
 * @@param: 
-->
## Align the homoeolog gene Coordinate 

### required 
1. bnMapper 
2. pysam
3. pybedtools

### Processes

1. Input file format

```bash
Ghir_A01        71162   72994   Ghir_A01G000040*-       Ghir_D01        42201   43973   Ghir_D01G000060*+       codeRegion
Ghir_A01        96241   99719   Ghir_A01G000070*+       Ghir_D01        79012   82604   Ghir_D01G000110*+       codeRegion
Ghir_A01        88884   89741   Ghir_A01G000080*+       Ghir_D01        67888   68745   Ghir_D01G000100*+       codeRegion
Ghir_A01        101839  104938  Ghir_A01G000100*+       Ghir_D01        84695   87883   Ghir_D01G000120*+       codeRegion
Ghir_A01        105797  107808  Ghir_A01G000110*+       Ghir_D01        88728   90735   Ghir_D01G000130*+       codeRegion
Ghir_A01        112915  113995  Ghir_A01G000130*+       Ghir_D01        93835   97091   Ghir_D01G000150*+       codeRegion
Ghir_A01        121080  121961  Ghir_A01G000160*+       Ghir_D01        104104  104984  Ghir_D01G000170*+       codeRegion
Ghir_A01        124197  125271  Ghir_A01G000170*-       Ghir_D01        106623  107585  Ghir_D01G000180*-       codeRegion
# command 
python ./homoeolog_alignment_lastz_v2.py inputData matchRegion.txt align_region_extend.txt
```

2. find diversity region
   > get diversity region in homoeolog pairs
```bash
python ./finde_diversity.py align_region_extend.txt matchRegion.txt diversity.txt
```

#### output file format

> match region File

1. chromosome1
2. start    (samtools faidx extract fasta sequence with start-1)
3. end 
4. chromoesom2
5. start    (samtools faidx extract fasta sequence with start-1)
6. end
7. gene Id 1 (gene*stand)
8. gene Id 2 (gene*stand)
9. alignment direction
10. alignment type (Match, Insert, Delete)


## align the QTL region in two genome by ucsc chain

```bash
gunzip testData/Ghir_D01.chain.gz
python QTL_conserved_region.py -QTL  QTL_region.txt   -alignment testData/Ghir_D01.chain  -myOut test_out
```

### input File

1. genome chain file 

2. QTL region file :
   1. the chromosome of QTL 
   2. lead SNP start
   3. lead SNP end 
   4. QTL regulated eGene pairs (eGene)
   5. QTL type (cis, trans, IntraChrAt, IntraChrDt, Inter)
   6. development stage
   7. QTL type (Bias QTL, eQTL)
   8. upstream of QTL region
   9. QTL region
   10. dowmstream of QTL region

### outFile:

1. the chromosome of QTL 
2. lead SNP start
3. lead SNP end 
4. QTL regulated eGene pairs (eGene)
5. QTL type (cis, trans, IntraChrAt, IntraChrDt, Inter)
6. development stage
7. QTL type (Bias QTL, eQTL)
8. upstream of QTL region
9. QTL region
10. dowmstream of QTL region
11. upstream QTL region correspond in another genome
12. QTL region correspond in another genome
13. downstream QTL region correspond in another genome

***
![](https://img.shields.io/badge/QTL%20align-Collinear-green)
## align QTL in two genome by gene collinear

```bash
   #! 1. get flank homoeolog geneID of QTL
   python QTL_flank_homoeolog.py QTL_coordinate.txt QTL_flank_homoeologId.txt
   #! 2. set a Id for a QTL
   sort -k1,1 -k2,2n QTL_flank_homoeologId.txt |awk '{print $0"\tQTL"NR}' >1
   mv 1 QTL_flank_homoeologId.txt
   #! 3. get align coordinate of QTL
   python QTL_align_coordinate.py QTL_flank_homoeologId.txt All_QTL_aligin_region.txt
   #! 4. align the region using lastz
   python QTL_align.py All_QTL_aligin_region.txt QTL_match_align.txt QTL_match_score.txt 
```

### the form of input file 

1. lead SNP of QTL coordinate File
   1.1 chromsome
   1.2 start
   1.3 end
```
    Ghir_A01        100011098       100011098
    Ghir_A01        10005837        10005837
    Ghir_A01        100086763       100086763
```
2. homoeolog gene Coordinate file
   2.1 chromsome
   2.2 start
   2.3 end
   2.4 homoeolog gene pair

'''bash
    Ghir_A01        70889   74180   Ghir_A01G000040-Ghir_D01G000060
    Ghir_D01        40781   44323   Ghir_D01G000060-Ghir_A01G000040
    Ghir_A01        87326   100744  Ghir_A01G000070-Ghir_D01G000110
    Ghir_D01        78908   83636   Ghir_D01G000110-Ghir_A01G000070
```
3. flank homoeolog geneId of QTL
   3.1 chromosome
   3.2 start
   3.3 end
   3.4 left-right(adjacent homoeolog gene)
   3.5 left-right(adjacent homoeolog gene in another genome)
   3.6 chromosomes
   3.7 chromosomes in another genome (may be Ghir_D01,Ghir_D02)
```bash
Ghir_A01        164527  164527  Ghir_A01G000250,Ghir_A01G000220 Ghir_D01G000260,Ghir_D01G000230 Ghir_A01        Ghir_D01
Ghir_A01        255346  255346  Ghir_A01G000410,Ghir_A01G000380 Ghir_D01G000410,Ghir_D01G000380 Ghir_A01        Ghir_D01
Ghir_A01        269159  269159  Ghir_A01G000410,Ghir_A01G000380 Ghir_D01G000410,Ghir_D01G000380 Ghir_A01        Ghir_D01
Ghir_A01        277130  277130  Ghir_A01G000410,Ghir_A01G000380 Ghir_D01G000410,Ghir_D01G000380 Ghir_A01        Ghir_D01
```
4. chromsomeSize File

'''bash
    Ghir_A01        117757855
    Ghir_A02        108092100
    Ghir_A03        113059412
    Ghir_A04        85149810
''' 

5. gene Coordinate File
'''bash
Ghir_A01        80323913        80324566        Ghir_A01G013980 +
Ghir_A01        80322152        80323048        Ghir_A01G013970 +
Ghir_A01        60753987        60754357        Ghir_A01G013100 +
Ghir_A01        59980270        59983471        Ghir_A01G013050 +
Ghir_A01        59140877        59148861        Ghir_A01G013020 +
Ghir_A01        60294264        60295580        Ghir_A01G013060 +
'''

6. QTL align coordinate

'''bash
    Ghir_A01        164027  165027  QTL1*+  Ghir_D01        145367  162194  Ghir_D01G000260-Ghir_D01G000230*+
    Ghir_A01        254846  255846  QTL2*+  Ghir_D01        243712  268247  Ghir_D01G000410-Ghir_D01G000380*+
    Ghir_A01        268659  269659  QTL3*+  Ghir_D01        243712  268247  Ghir_D01G000410-Ghir_D01G000380*+
    Ghir_A01        276630  277630  QTL4*+  Ghir_D01        243712  268247  Ghir_D01G000410-Ghir_D01G000380*+
''' 