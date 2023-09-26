#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 25 16:25:50 2023

@author: rafaella
"""

import os
import utils.preprocessing as prep
import sys
import scanpy as sc
sc.settings.set_figure_params(dpi=100, facecolor='white')

#Initialize global variables
nCores = 40 # max 40
nRAM = 100 # max 200

#Path to reads
pathToFastqFiles ='data/Sample_MUC30433/raw'

#FASTQC
prep.doFastQC (pathToFastqFiles, nCores)

#CELLRANGER 
pathToRefTranscriptome = 'data/resources/refdata-gex-GRCh38-2020-A'
expectCells = 10000
pathToCellRangerResults = prep.doCellRangerSC (pathToFastqFiles, pathToRefTranscriptome, 
                                                expectCells, nCores=40, nRAM=200)

#pathToCellRangerResults = pathToFastqFiles.split('raw')[0]+'cellranger/MUC30433/outs/filtered_feature_bc_matrix'
#Load into scanpy AnnData object
adata = sc.read_10x_mtx(pathToCellRangerResults, var_names='gene_symbols')

#Basic filtering
adata = prep.filter (adata)

#Normalize
adata = prep.normalize(adata)
#Keep only highly variable
adata = prep.highlyVarGenes(adata)
#ComBat batch correction
sc.pp.combat(adata)
adata.corrected = adata


#PCA
sc.tl.pca(adata, svd_solver='arpack')
sc.pl.pca(adata, color='CST3')