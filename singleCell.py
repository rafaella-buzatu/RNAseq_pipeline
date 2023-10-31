#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 25 16:25:50 2023

@author: rafaella
"""


import utils.preprocessing as prep

# Initialize global variables
nCores = 40  # max 40
nRAM = 100  # max 200

# Path to reads
pathToFastqFiles = 'data/Sample_MUC30433/raw'

# FASTQC
prep.doFastQC(pathToFastqFiles, nCores)

# CELLRANGER
pathToRefTranscriptome = 'data/resources/refdata-gex-GRCh38-2020-A'
expectCells = 10000
pathToCellRangerResults = prep.doCellRangerSC(pathToFastqFiles, pathToRefTranscriptome,
                                              expectCells, nCores=40, nRAM=200)


#opening bam file
import pysam
samfile = pysam.AlignmentFile("data/Sample_MUC30433/cellranger/MUC30433/outs/possorted_genome_bam.bam", "rb")
iterator = samfile.fetch("chr22", 32507820, 33058381)