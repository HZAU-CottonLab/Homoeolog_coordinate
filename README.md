<!--
 * @Descripttion: 
 * @version: 
 * @Author: zpliu
 * @Date: 2021-07-04 10:19:51
 * @LastEditors: zpliu
 * @LastEditTime: 2021-07-04 10:27:21
 * @@param: 
-->
## Align the homoeolog gene Coordinate 

### required 
1. bnMapper 
2. pysam
3. pybedtools

### Processes

1. get bnMapper result in homoeolog region

```bash
python get_homoeolog_bnMapper.py -homoeolog homoeologGeneRegionfile -myOut MapperRegion.txt 
```

2. filter bnMapper result 

```bash
python filterBnMapper.py homoeolog_bnMapper.txt filter.txt 
```

3. align the homoeolog region

```bash
python homoeolog_alignment_lasstz.py  homoeologGeneRegionfile outfile
```

### results

> **homoeolog_gene_conserved_flagment.txt**

+ chromosomeId A
+ start A
+ end A
+ chromeomeId B
+ start B
+ end B
+ geneId A
+ geneId B
+ lastz score Huge difference in gene sequence with 0
+ align direction of alignment