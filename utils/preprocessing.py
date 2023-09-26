#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 21 14:08:36 2023

@author: rafaella
"""

import subprocess
import os
import numpy as np
import tqdm
from collections import defaultdict
import pandas as pd
import csv
import gzip
import scipy.io
import scanpy as sc
sc.settings.set_figure_params(dpi=100, facecolor='white')




def doFastQC (pathToFastqFiles, nCores = 40):
    '''
    Runs a FastQC assessment on the fastq files in the input folder.

    Parameters
    ----------
    pathToFastqFiles : A string specifying the path to the folder containing the
    fastq files.
                     
    nCores : An int specifying the number of threads to use. Default is 20

    Returns
    -------
    No return. The results of the FastQC will be saved in a folder named 'fastqc'
    outside of the 'raw' directory.
    '''
    #Create the directory in which to store the fastqc results
    pathToOutputDir = os.path.join (pathToFastqFiles.split('raw')[0], 'fastqc')
    if not os.path.exists(pathToOutputDir):
        os.makedirs(pathToOutputDir)
    
    #Define bash command to call the fastqc software
    command = "fastqc -o {} -t {} {}".format(pathToOutputDir, str(nCores),
                                             pathToFastqFiles + '/*.gz')
    #Send command to console
    subprocess.check_call(command, shell=True)

################ SINGLE CELL

def doCellRangerSC (pathToFastqFiles, pathToRefTranscriptome, expectedCells,
                    nCores=40, nRAM=200):
    '''
    Runs cellRanger on the single cell fastq files of the input directory.

    Parameters
    ----------
    pathToFastqFiles : A string specifying the path to the folder containing the
    fastq files.
    
    pathToRefTranscriptome : A string specifying the path to the folder containing 
    the reference transcriptome.
    
    expectedCells : An int specifying the number of expected cells.
    
    nCores : An int specifying the number of threads to use. Default is 40.
    
    nRAM : An int specifying the amount of RAM memory to use. Default is 200.
    
    raw: A boolean spcifying whether to output the path to the folder with the
    raw or filtered counts.

    Returns
    -------
    The path to the folder containing the count matrix elements.

    '''
    
    #Create output directory
    pathToOutputDir = os.path.join (pathToFastqFiles.split('raw')[0], 'cellranger')
    if not os.path.exists(pathToOutputDir):
        os.makedirs(pathToOutputDir)
        
    #Extract sample name from file name
    sampleName = [file for file in os.listdir (pathToFastqFiles) if file.endswith('.gz')][0].split('_')[0]
   
    #Call cellranger
    command = "cellranger count --id={}\
                                --transcriptome={}\
                                --sample={} \
                                --localcores={} \
                                --expect-cells={} \
                                --localmem={} \
                                --fastqs={}".format (sampleName, 
                                                     pathToRefTranscriptome,
                                                     sampleName,
                                                     nCores,
                                                     expectedCells,
                                                     nRAM,
                                                     pathToFastqFiles
                                )
    subprocess.check_call(command, shell=True)

    #Move to data folder
    command ="mv {} {}".format(sampleName, pathToOutputDir)
    #Define path to the output feature matrix
    pathToMatrices = os.path.join(pathToOutputDir, sampleName, 'outs', 'filtered_feature_bc_matrix')
    
    return pathToMatrices
    
### COUNT MATRIX PROCESSING
def extractCountMatrixSC (pathToCellRangerResults, save = True):
    '''
    Creates a read count matrix from the cellranger output data.

    Parameters
    ----------
    pathToCellRangersResults : A string specifying the path to the folder
    containing the counts 
        
    save: A boolean specifying whether to save the resulting matrix

    Returns
    -------
    The read count matrix.

    '''
   
    #Read count matrix 
    matrix = scipy.io.mmread(os.path.join(pathToCellRangerResults, "matrix.mtx.gz"))

    #Define path to features table
    featuresPath = os.path.join(pathToCellRangerResults, "features.tsv.gz")
    #Read transcript IDs
    featureIDs = [row[0] for row in csv.reader(gzip.open(featuresPath, mode="rt"), delimiter="\t")]
    #Read gene names
    geneNames = [row[1] for row in csv.reader(gzip.open(featuresPath, mode="rt"), delimiter="\t")]
    #Read feature types 
    featureTypes = [row[2] for row in csv.reader(gzip.open(featuresPath, mode="rt"), delimiter="\t")]
    
    #Read barcodes
    barcodesPath = os.path.join(pathToCellRangerResults, "barcodes.tsv.gz")
    barcodes = [row[0] for row in csv.reader(gzip.open(barcodesPath, mode="rt"), delimiter="\t")]
    
    # transform table to pandas dataframe and label rows and columns
    countMatrix = pd.DataFrame.sparse.from_spmatrix(matrix)
    countMatrix.columns = barcodes
    countMatrix.insert(loc=0, column="feature_id", value=featureIDs)
    countMatrix.insert(loc=0, column="gene", value=geneNames)
    countMatrix.insert(loc=0, column="feature_type", value=featureTypes)

    # save the count matrix as a csv
    if (save):
        pathToResults = pathToCellRangerResults.split('cellranger')[0]
        countMatrix.to_csv(pathToResults + "countMatrix.csv", index=False)
    
    return countMatrix

############ BULK

def bulkPreprocessing (pathToFastqFiles, pathToRefAdapters, pathToRefPhix, 
                       trimQuality = 30, minLen = 30, nRAM= 200):
    '''
    Runs the following preprocessing steps on the input bulk RNAseq:
        1.Trim sequencing adapters
        2.Remove phage spike in
        3.Trim based on quality score
        4.Remove reads shorter than 30bp

    Parameters
    ----------
    pathToFastqFiles : A string specifying the path to the folder containing the
    fastq files.
    
    pathToRefAdapters : A string specifying the path to the .fa (Fasta) file
    containing the adapter sequences.
    
    pathToRefPhix : A string specifying the path to the .fa (Fasta) file
    containing the Phix reference genome.
    
    trimQuality : An int specifying the quality score below which to trim. The 
    default is 30.
    
    minLen : An int specifying the minimum length of a transcript to be kept. 
    The default is 30.
        
    nRAM : An int specifying the amount of RAM memory to use. Default is 200.


    Returns
    -------
    The path to the folder in which the processed bulk RNAseq fastq files are 
    stored.

    '''
    
    #Create temporary directory to store processed files
    pathToTemp = os.path.join(pathToFastqFiles.split('raw')[0], 'temp')
    if not os.path.exists(pathToTemp):
        os.makedirs(pathToTemp)
    
    #Get sample names from input folder
    bulkRNAseqSampleNames = np.unique([file.split('_')[0] for file in os.listdir (pathToFastqFiles) if file.endswith('.gz')]).tolist()
    
    #Run pre-processing steps per sample
    for sample in tqdm.tqdm(bulkRNAseqSampleNames) :
        #Trim adapter sequences
        doAdapterTrim(pathToFastqFiles, pathToTemp, sample, pathToRefAdapters, nRAM)
        #Trim reads based on quality score and minimum length
        doQualityTrim(pathToTemp, pathToTemp, sample, trimQuality, minLen, nRAM)
        #Remove phage spikein reads
        doPhageTrim(pathToTemp, pathToTemp, sample, pathToRefPhix, nRAM)
    
    #Create folder to store processed files
    pathToProcessed = os.path.join(pathToFastqFiles.split('raw')[0], 'processed')
    if not os.path.exists(pathToProcessed):
        os.makedirs(pathToProcessed)
    
    #Move processed files to a new folder and delete the temporary folder
    command = "mv {}/*_processed.fastq.gz {}; rm -r {}; ".format (
        pathToTemp, pathToProcessed, pathToTemp)
    subprocess.check_call(command, shell=True)
    
    return pathToProcessed
    
### BULK PREPROCESSING FUNCTIONS

def doAdapterTrim(pathToInput, pathToOutput, sample, pathToRefAdapters, nRAM ):
    '''
    Perform adapter trimming on bulk RNAseq reads based on input adapter sequences.

    Parameters
    ----------
    pathToInput :  string specifying the path to the folder containing the
    fastq files.
    
    pathToOutput : A string specifying the path to the directory in which the
    processed files will be stored.
    
    sample : A string specifying the name of the sample to process.
    
    pathToRefAdapters : A string specifying the path to the .fa (Fasta) file
    containing the adapter sequences.
        
    nRAM : An int specifying the amount of RAM memory to use.
    Returns
    -------
    Nothing is returned. The processed files are saved in the pathToOutput directory.

    '''
    #Generate command to call bbduk
    command = "bbduk.sh -Xmx{}g in1={} in2={} out1={} out2={} \
        ktrim=r k=21 hdist=1 mink=6 ref={} ordered=t".format (
        nRAM,
        '/'.join([pathToInput, sample+'_1.fastq.gz']), 
        '/'.join([pathToInput, sample+'_2.fastq.gz']),
        '/'.join([pathToOutput, sample+'_1_trimmed.fastq']),
        '/'.join([pathToOutput, sample+'_2_trimmed.fastq']),
        pathToRefAdapters
        )
    #Send command to console
    subprocess.check_call(command, shell=True)
    
    
def doQualityTrim(pathToInput, pathToOutput, sample, trimQuality, minLen, nRAM):
    '''
    

    Parameters
    ----------
    pathToInput :  string specifying the path to the folder containing the
    fastq files.
    
    pathToOutput : A string specifying the path to the directory in which the
    processed files will be stored.
    
    sample : A string specifying the name of the sample to process.
    
    trimQuality : An int specifying the quality score below which to trim. 
    
    minLen : An int specifying the minimum length of a transcript to be kept. 
        
    nRAM : An int specifying the amount of RAM memory to use.


    Returns
    -------
    Nothing is returned. The processed files are saved in the pathToOutput directory.
    
    '''
    #Generate command to call bbduk
    command = "bbduk.sh -Xmx{}g in1={} in2={} out1={} out2={} trimq={} minlen={}".format(
        nRAM,
        '/'.join([pathToInput, sample+'_1_trimmed.fastq']),
        '/'.join([pathToInput, sample+'_2_trimmed.fastq']),
        '/'.join([pathToOutput, sample+'_1_quality.fastq']),
        '/'.join([pathToOutput, sample+'_2_quality.fastq']),
        trimQuality,
        minLen)
    #Send command to console
    subprocess.check_call(command, shell=True)
    
    
def doPhageTrim(pathToInput, pathToOutput, sample, pathToRefPhix, nRAM):
    '''
    

    Parameters
    ----------
    pathToInput :  string specifying the path to the folder containing the
    fastq files.
    
    pathToOutput : A string specifying the path to the directory in which the
    processed files will be stored.
    
    sample : A string specifying the name of the sample to process.
    
    pathToRefPhix : A string specifying the path to the .fa (Fasta) file
    containing the Phix reference genome.
    
    nRAM : An int specifying the amount of RAM memory to use.


    Returns
    -------
    Nothing is returned. The processed files are saved in the pathToOutput directory.

    '''
    #Generate command to call bbduk
    command = "bbduk.sh -Xmx{}g in1={} in2={} out1={} out2={} \
               outm1={} outm2={} ref={} k=31 hdist=1 stats={}".format(
        nRAM,
        '/'.join([pathToInput,  sample+'_1_quality.fastq']),
        '/'.join([pathToInput, sample+'_2_quality.fastq']),
        '/'.join([pathToOutput, sample+'_1_processed.fastq.gz']),
        '/'.join([pathToOutput, sample+'_2_processed.fastq.gz']),
        '/'.join([pathToOutput, sample+'_1_matched.fastq']),
        '/'.join([pathToOutput, sample+'_2_matched.fastq']),
        pathToRefPhix,
        '/'.join([pathToOutput, 'stats.txt'])
        )
    #Send command to console
    subprocess.check_call(command, shell=True)


### ALIGNMENT SCRIPTS

def createGenomeIndices (pathToGenomeReference, pathToAnnotations, pathToGenomeDir,
                         overhang = 100, nCores = 40):
    '''
    Generates Genome Indices folder as required by the STAR Aligner.

    Parameters
    ----------
    pathToGenomeReference : A string specifying the path to the .fna genome 
    reference file.
    
    pathToAnnotations : A string specifying the path to the .gtf genome annotations 
    file.
    
    pathToGenomeDir : A string specifying the path to the directory in which the
    resulting files will be saved.
    
    overhang : An int equal to ReadLength-1 The default is 100.
        
    nCores : An int specifying the number of threads to use. The default is 40.

    Returns
    -------
    Nothing is returned. 

    '''
    
    #Create directory to store genome indices
    if not os.path.exists(pathToGenomeDir):
        os.makedirs(pathToGenomeDir)
    
    #Generate genome indices in specified folder
    command = "STAR --runThreadN {}\
                    --runMode genomeGenerate\
                    --genomeDir {}\
                    --genomeFastaFiles {}\
                    --sjdbGTFfile {}\
                    --sjdbOverhang {}".format(
                    nCores,
                    pathToGenomeDir,
                    pathToGenomeReference,
                    pathToAnnotations,
                    overhang)
    
    subprocess.check_call(command, shell=True)
    subprocess.check_call("rm Log.out", shell=True)
    

def runSTARaligner (pathToProcessed, pathToGenomeDir, nCores = 40):
    '''
    Runs the STAR aligner on each sample in the input directory.

    Parameters
    ----------
    pathToProcessed : A string specifying the path to the directory containing
    the fastq files to be aligned.
    
    pathToGenomeDir : A string specifying the path to the genome index directory.
    
    nCores : An int specifying the number of threads to use. The default is 40.

    Returns
    -------
    pathToOutputDir : A string specifying the path to the folder created by the
    STAR aligner in which the results are stored.

    '''
    
    #Create the directory in which to store the star alignment results
    pathToOutputDir = os.path.join (pathToProcessed.split('processed')[0], 'STARaligner')
    if not os.path.exists(pathToOutputDir):
        os.makedirs(pathToOutputDir)
    
    #Extract sample names
    bulkRNAseqProcessed = np.unique([file.split('_')[0] for file in os.listdir (pathToProcessed) if file.endswith('.gz')]).tolist()

    #Load genome indices
    subprocess.check_call("STAR --genomeLoad LoadAndExit --genomeDir {}".format (pathToGenomeDir), shell=True)
    
    for sample in bulkRNAseqProcessed:
        #Generate string of file names to feed into the aligner 
        sampleReads = sorted([os.path.join(pathToProcessed,file) for file in os.listdir(pathToProcessed) if file.startswith(sample)])
        readFilesString =' '.join([file for file in sampleReads])
    
    
        #Call aligner
        command = "STAR --runThreadN {}\
                        --genomeDir {}\
                        --readFilesCommand zcat\
                        --readFilesIn {}\
                        --outFileNamePrefix {}\
                        --genomeLoad LoadAndKeep\
                        --outSAMtype None\
                        --quantMode GeneCounts".format(
                        nCores,
                        pathToGenomeDir,
                        readFilesString,
                        pathToOutputDir + '/' + sample +'_')
        subprocess.check_call(command, shell=True)
    
    #Remove genome indices from memory
    subprocess.check_call("STAR --genomeLoad Remove --genomeDir {}".format (pathToGenomeDir), shell=True)
    
    #Remove temporary files
    #subprocess.check_call("rm *.out*", shell=True)
    #subprocess.check_call("rm -r _STARtmp", shell=True)

    
    return pathToOutputDir

### COUNT MATRIX PROCESSING
def extractCountMatrixBulk (pathToAlignment):
    '''
    Extracts the necessary data and creates a count matrix based on the results
    of the alignment.

    Parameters
    ----------
    pathToAlignment : A string specifying the path to the folder created by
    the STAR aligner.

    Returns
    -------
    The extracted count matrix.
    '''

    #Get filenames of count matrices
    files = [os.path.join (pathToAlignment, file) for file in os.listdir (pathToAlignment) if file.endswith('ReadsPerGene.out.tab')]
    
    #Create empty dictionary to store counts in
    #Format -> [countMatrix[GeneName][Sample] = count]
    countMatrix = defaultdict(lambda: defaultdict(int))
    
    #Iterate over all samples
    for file in files:
        #Extract sample name
        sample = file.split('/')[-1].split('_')[0]
        #Read the file
        with open(file, newline = '') as tableReads:                                                                                          
            csvreader = csv.reader(tableReads, delimiter='\t')
            listOfGenes = []
            #Read lines into lists 
            for line in csvreader:
                listOfGenes.append(line)
            #Remove first lines
            listOfGenes = listOfGenes[4:]
            
            #Extract gene names and both-stranded counts from nested list
            geneNames = [item[0] for item in listOfGenes]
            counts = [item[1] for item in listOfGenes]
            #Append counts to specific gene and sample entry in the dictionary
            for i in range(len(geneNames)):
                countMatrix[geneNames[i]][sample] = counts[i]
            
    #Convert the dictionary into a dataframe
    countMatrixDF = pd.DataFrame.from_dict(countMatrix, orient='index')
    #Save count matrix df
    countMatrixDF.to_csv(pathToAlignment.split('STARaligner')[0] +'countMatrix.csv')

    return countMatrixDF

##### SCANPY

def filter (adata, minGenesFilter = 200, minCellsFilter = 3, 
            maxCountsFilter = 2500, mitCountsFilter = 5):
    '''
    Performs filtering of the count matrix

    Parameters
    ----------
    adata : An AnnData object to filter
    
    minGenesFilter : An int specifying the min amount of genes a cell should 
    have expression values for to not be filtered out. The default is 200.
    
    minCellsFilter : An int specifying the min amount of cells in which a gene
    should be expressed to not be filtered out.The default is 3.
    
    maxCountsFilter : An int specifying the max number of gene counts a cell
    should have to not be filtered out. The default is 2500 for single cell.
    Should be removed for bulk.
    
    mitCountsFilter : An int specifying the max amount of mitochondrial genes a
    cell should express to not be filtered out. The default is 5.

    Returns
    -------
    The filtered AnnData object.
    '''
    
    #Basic filtering of cells with less than 200 genes 
    sc.pp.filter_cells(adata, min_genes=minGenesFilter)
    #Basic filtering of genes expressed in less than 3 cells
    sc.pp.filter_genes(adata, min_cells= minCellsFilter)
    
    #Annotate mitochondiral genes
    adata.var['mt'] = adata.var_names.str.startswith('MT-')
    #Calculate qc metrics 
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
    
    #Remove cells with too many mitochondrial genes -> leakage cytoplasmic mRNA
    adata = adata[adata.obs.pct_counts_mt < mitCountsFilter, :]
    #Remove cells with too many total counts -> doublets
    if (maxCountsFilter is not None):
        adata = adata[adata.obs.n_genes_by_counts < maxCountsFilter, :]
    
    return adata


def normalize (adata):
    '''
    Normalize input AnnData count matrix.
    '''
    
    #Scale raw counts by 10000 and log(x+1) transform
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    #Set the raw 
    adata.raw = adata
    
    return adata

def highlyVarGenes (adata, plot = False):
    '''
    Filters out the genes that are not highly variable


    '''
    #Find highly variable genes
    sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
    if (plot):
        sc.pl.highly_variable_genes(adata)
    #Keep only highly variable genes
    adata = adata[:, adata.var.highly_variable]

    return adata