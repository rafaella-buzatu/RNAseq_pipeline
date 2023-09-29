#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 20 11:19:40 2023

@author: rafaella
"""

import os
import utils.preprocessing as prep
import sys
import pandas as pd
import scanpy as sc
import pydeseq2
sc.settings.set_figure_params(dpi=100, facecolor='white')


#Initialize global variables
nCores = 40 # max 40
nRAM = 100 # max 200


pathToFastqFiles = 'data/ytTutorial/raw'

#FASTQC
prep.doFastQC (pathToFastqFiles, nCores)

####Quality control
pathToRefAdapters = 'data/resources/adapters.fa'
pathToRefPhix = 'data/resources/phix174_ill.ref.fa.gz'

#Do pre-processing:
    #Adapter trimming
    #Remove phage spike in
    #Trim based on quality score
    #Remove reads shorter than 30bp
pathToProcessed = prep.bulkPreprocessing (pathToFastqFiles, pathToRefAdapters, 
                                          pathToRefPhix, trimQuality = 10, 
                                          minLen = 30, nRAM= 100)

###Alignment

#Create Genome Indices
pathToGenomeReference  = 'data/resources/GCF_000001405.40_GRCh38.p14_genomic.fna'
pathToAnnotations = 'data/resources/GCF_000001405.40_GRCh38.p14_genomic.gtf'
overhang = 50
pathToGenomeDir = 'data/resources/genomeIndices'

prep.createGenomeIndices (pathToGenomeReference, pathToAnnotations, pathToGenomeDir, overhang)

#Run alignment
pathToAlignment = prep.runSTARaligner (pathToProcessed, pathToGenomeDir, nCores = 40)

#Extract count matrix after alignment
countMatrix = prep.extractCountMatrixBulk (pathToAlignment)
#Read count matrix from file



### DIFFERENTIAL EXPRESSION
pathToCountMatrix = pathToFastqFiles.split('raw')[0]+'countMatrix.csv'
countMatrix = pd.read_csv(pathToCountMatrix, index_col = 0)

#Remove genes where counts are 0 for all samples
countMatrix = countMatrix[countMatrix.sum(axis = 1) >0]
#Transpose
countMatrix = countMatrix.T

#metadata dataframe
#index: sample name
#metadata = 

#Read count matrix into AnnData object
dds = pydeseq2.dds.DeseqDataSet(counts = countMatrix,
                                clinical = metadata,
                                design_factors= "condition") #column name from from the metadata table 
#Run differential expression 
dds.deseq2()

#Extract data
statRes = pydeseq2.ds.DeseqStats(dds, n_cpus = nCores, 
                                 contrast = ('condition',
                                             'condition1', 'condition2')) #values of column
statRes.summary()
#Get results into a dataframe
DEgenes = statRes.results_df

#Extract differentially expressed genes
DEgenes = DEgenes[(DEgenes.padj < 0.05) & (abs(DEgenes.log2FoldChange) > 0.5) &(DEgenes.baseMean>=10)]

#Do PCA
sc.tl.pca(dds)
sc.pl.pca(dds, color = 'condition')
