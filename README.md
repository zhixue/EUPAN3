# EUPAN3
(The project is under development.)
## 1. Introduction
EUPAN3 is a EUkaryote Pan-genome ANalysis toolkit for 3rd generation sequencing.

The overview of pipeline is as follows:
![pipeline](https://raw.githubusercontent.com/zhixue/EUPAN3/master/pics/pipeline.jpg)

## 2. Install
### Dependency
* [Quast](https://github.com/ablab/quast) (>=v5.0.2)
* [Minimap2](https://github.com/lh3/minimap2) (>=v2.17)
* [Cd-hit](https://github.com/weizhongli/cdhit) (>=v4.8.1)
* [Gclust](https://github.com/niu-lab/gclust) (>=v1.0)
* [Blast](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST) (>=2.11.0)
* [Bedtools](https://github.com/arq5x/bedtools2) (>=2.29.2)
### 

## 3. Usage
```
usage: eupan3.py [-h] [-v]
                 {assemsta,unalnbseq,mergebseq,rmredundant,rmctm,fastasta,ptpg,genecov,elecov}
                 ...

EUPAN3: EUkaryote Pan-genome ANalysis toolkit for 3rd generation sequencing
Version 0.1.0

optional arguments:
  -h, --help            show this help message and exit
  -v, --version         show program's version number and exit

sub commands:
  {assemsta,unalnbseq,mergebseq,rmredundant,rmctm,fastasta,ptpg,genecov,elecov}
    assemsta            Use Quast to get assemblies statistics
    unalnbseq           Get partial/full unaligned block sequences
    mergebseq           Merge the block sequences from many files
    rmredundant         Cluster the sequences and remove redundant sequences
    rmctm               Detect and discard the contaminated sequences
    fastasta            Calculate statistics of fasta file (DNA)
    ptpg                Extract the longest transcript elements from genes
    genecov             Compute gene coverage and depth
    elecov              Compute genomic element coverage and depth
```

## 4. Quick start
Under development

## 5. Documentation
The documentation can be viewed at [wiki](https://github.com/zhixue/EUPAN3/wiki).