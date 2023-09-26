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
sc.settings.set_figure_params(dpi=100, facecolor='white')


#Initialize global variables
nCores = 40 # max 40
nRAM = 100 # max 200


pathToFastqFiles = 'data/bulk_RNAseq/raw'

#FASTQC
prep.doFastQC (pathToFastqFiles, nCores)

#Quality control
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

#Alignment

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
pathToCountMatrix = pathToFastqFiles.split('raw')[0]+'countMatrix.csv'
countMatrix = pd.read_csv(pathToCountMatrix, index_col = 0)

#Read count matrix into AnnData object
adata= sc.AnnData(countMatrix.T)

#Basic filtering
adata = prep.filter (adata, minCellsFilter = 1, maxCountsFilter = None)
#Normalize
adata = prep.normalize(adata)
#Keep only highly variable
adata = prep.highlyVarGenes(adata)

