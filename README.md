<!--
 * @Descripttion: 
 * @version: 
 * @Author: zpliu
 * @Date: 2021-07-04 10:19:51
 * @LastEditors: zpliu
 * @LastEditTime: 2021-07-06 16:02:51
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
