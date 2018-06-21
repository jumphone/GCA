# GCA: Gene Correlation Analysis in Single-cell RNA-seq data

Author: Feng Zhang

Email: 15110700005@fudan.edu.cn

Date: 2018.06.12

# Demo: 

Human Oligodendrogliomas: http://omics.fudan.edu.cn/static/demo/GCA/GSE70630/index.html

# User Guide:

Requirements:

    Operation system: Mac or Linux, python=2.7, R=3.5
    
    Python packages: scipy
    
R packages: stringr, pastecs, mixtools, igraph

Command options:

i) Expression index building (Python script)

python | step1_BuildIndex.py | NETWORK | EXP | OUTPUT_INDEX | RATE

    NETWORK: a list of tab-delimited gene pairs
    
    EXP: normalized expression matrix (without quotation marks)
    
    OUTPUT_INDEX: the index file name, given by the user
    
    RATE: the cutoff for the non-zero rate. For each gene pair, non-zero rate equals to the proportion of cells of which two genes both are detected. 

ii) Cell-specific network building (Python script)

python | step2_CalZmat.py | OUTPUT_INDEX | OUTPUT_ZMAT | CPU
    
    OUTPUR_INDEX: the index file generated in the first step
    
    OUTPUT_ZMAT: the z-value matrix file name, given by the user
    
    CPU: the number of threads to run GCA 

iii) Mixture models analyzing and result generating (R script).

Rscript | step3_deMix.R | EXP | OUTPUT_ZMAT | OUTPUT | CPU | SEED | CUTOFF
    
    EXP: the expression matrix used in the first step

    OUTPUT_ZMAT: the z-value matrix generated in the second step

    OUTPUT: the name of the output directory, given by the user

    CPU: the number of threads to run GCA 

    SEED: seed for random function in R

    CUTOFF: the cutoff of edge score to draw the gene graph
