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
pathToFastqFiles = 'data/rnaseq/Sample_MUC30433'
prep.doFastQC (pathToFastqFiles)


### Single cell

#CELLRANGER 
pathToFastqFiles = 'data/rnaseq/Sample_MUC30433'
pathToRefTranscriptome = '../software/refdata-gex-GRCh38-2020-A'
expectCells = 10000
prep.doCellRangerSC (pathToFastqFiles, pathToRefTranscriptome,
                                   expectCells,nCores=40, nRAM=200)

####  BULK

pathToFastqFiles = 'data/rnaseq/bulk_RNAseq'
pathToRefAdapters = 'data/resources/adapters.fa'
pathToRefPhix = 'data/resources/phix174_ill.ref.fa.gz'

pathToProcessed = prep.bulkPreprocessing (pathToFastqFiles, pathToRefAdapters, 
                                          pathToRefPhix, trimQuality = 10, 
                                          minLen = 30, nRAM= 100)
print (os.listdir(pathToProcessed))



