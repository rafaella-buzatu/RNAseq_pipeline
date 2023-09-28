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
import numpy as np
import diffxpy.api as de
import pandas as pd
sc.settings.set_figure_params(dpi=100, facecolor='white')

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

pathToCellRangerResults = pathToFastqFiles.split(
    'raw')[0]+'cellranger/MUC30433/outs/filtered_feature_bc_matrix'
# Load into scanpy AnnData object
adata = sc.read_10x_mtx(pathToCellRangerResults, var_names='gene_symbols')


#Create folder and file to store results
pathToResults = os.path.join ('results', pathToFastqFiles.split('/')[1])
if not os.path.exists(pathToResults):
    os.makedirs(pathToResults)
pathToResultsFile =  os.path.join (pathToResults, 'analysisResults.h5ad')

# Basic filtering
adata = prep.filterData(adata)
# Normalize
adata = prep.normalize(adata)
# Keep only highly variable
adata = prep.highlyVarGenes(adata)


pathToCellCycleGenes = 'data/resources/regev_lab_cell_cycle_genes.txt'

adata = prep.correction(adata, regressTotalCounts=True, regressCellCycle=False,
                       pathToCellCycleGenes=pathToCellCycleGenes, regressMit=True)

# Scale the data to unit variance; Clip values where stdev> max_value
sc.pp.scale(adata, max_value=10)

# do PCA
sc.tl.pca(adata, svd_solver='arpack')
#Plot PCs
sc.pl.pca(adata, color='CST3')
#Plot %variance explained
sc.pl.pca_variance_ratio(adata, log=True)
elbowPCs = 10


### CLUSTERING
#Compute neighbourhood graph
sc.pp.neighbors(adata, n_pcs=elbowPCs)

#Plot umaps coloured by gene expression
sc.tl.umap(adata)
listGenesToPlot = adata.var.highly_variable[:10].index.tolist()
sc.pl.umap(adata, color=listGenesToPlot, use_raw=False)

#Use Leiden for graph clustering
sc.tl.leiden(adata)
sc.pl.umap(adata, color=['leiden'], legend_loc = "on data")

#Find marker genes per cluster
sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')
sc.pl.rank_genes_groups(adata, n_genes=20, sharey=False)
#Extract clustering results in dataframe
markers = sc.get.rank_genes_groups_df(adata, None)
#Keep genes with adj pvals < 0.05 and logfold cahanges > 0.5
markers  = markers[(markers.pvals_adj < 0.05) & (markers.logfoldchanges > 0.5)]
#Get no of clusters
noClusters = np.unique(markers.group).tolist()

### TOOLS TO ANNOTATE
#Can find cell type markers here
#https://panglaodb.se/markers.html?cell_type=%27choose%27
#Plot genes
sc.pl.umap(adata, color=['leiden', 'CD68'], legend_loc = "on data")
#Check logfoldchange of certain gene
markers[markers['names'] == 'CD68']
#Find genes with highest logfoldchanges
markers[markers['group'] == '1'].sort_values(['logfoldchanges'], ascending = False)

#Annotation: (manual)
dictClusterAnnotation = {key: None for key in noClusters}
#Update dictionary:
dictClusterAnnotation['0'] = 'annotation'
dictClusterAnnotation['1'] = 'annotation2'


#Map annotation to samples
adata.obs['cell_type'] = adata.obs.leiden.map(dictClusterAnnotation)
sc.pl.umap(adata, color=['cell_type'], legend_loc = "on data")


adata.write(pathToResultsFile)
adata = sc.read(pathToResultsFile)

#### DIFFERENTIAL EXPRESSION

#Subset cell types
cellType1 = 'annotation2'
cellType2 = 'annotation'
subset = adata[adata.obs['cell_type'].isin([cellType1, cellType2])].copy()
#Convert back to normalized data
subset = subset.raw.to_adata()
subset.X = subset.X.toarray()
#Perform differential expression analysis based on cell type
res = de.test.wald(data = subset, 
                   formula_loc = '~1 + cell_type',
                   factor_loc_totest = 'cell_type')
DEdf = res.summary