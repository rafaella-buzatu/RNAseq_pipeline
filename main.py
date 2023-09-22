#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 20 11:19:40 2023

@author: rafaella
"""

import os
import utils.preprocessing as prep

#Initialize global variables
nCores = 40 # max 40
nRAM = 100 # max 200

#FASTQC
pathToFastqFiles = 'data/Sample_MUC30433/raw'
prep.doFastQC (pathToFastqFiles, nCores)


### Single cell

#CELLRANGER 
pathToFastqFiles = 'data/Sample_MUC30433/raw'
pathToRefTranscriptome = 'data/resources/refdata-gex-GRCh38-2020-A'
expectCells = 10000
prep.doCellRangerSC (pathToFastqFiles, pathToRefTranscriptome, expectCells,
                     nCores=40, nRAM=200)

####  BULK

pathToFastqFiles = 'data/bulk_RNAseq/raw'

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
countMatrix = prep.extractCountMatrix (pathToAlignment)

