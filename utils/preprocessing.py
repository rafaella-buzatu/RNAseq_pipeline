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



def doFastQC (pathToFastqFiles, nCores = 20):
    
    pathToOutputDir = os.path.join (pathToFastqFiles, 'fastqc')
    if not os.path.exists(pathToOutputDir):
        os.makedirs(pathToOutputDir)

    command = "fastqc -o {} -t {} {}".format(pathToOutputDir, str(nCores),
                                             pathToFastqFiles + '/*.gz')
    subprocess.check_call(command, shell=True)


def doCellRangerSC (pathToFastqFiles, pathToRefTranscriptome, expectCells,
                    nCores=40, nRAM=200):
    
    #Create output directory
    pathToOutputDir = os.path.join (pathToFastqFiles, 'cellranger')
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
                                                     expectCells,
                                                     nRAM,
                                                     pathToFastqFiles
                                )
    subprocess.check_call(command, shell=True)

    #Move to results folder
    command ="mv {} {}".format(sampleName, pathToOutputDir)
    

def bulkPreprocessing (pathToFastqFiles, pathToRefAdapters, pathToRefPhix, 
                       trimQuality = 10, minLen = 30, nRAM= 100):
    
    #Create temporary directory to store processed files
    pathToTemp = os.path.join(pathToFastqFiles, 'temp')
    if not os.path.exists(pathToTemp):
        os.makedirs(pathToTemp)
    
    #Get sample names from input folder
    bulkRNAseqSampleNames = np.unique([file.split('_')[0] for file in os.listdir (pathToFastqFiles) if file.endswith('.gz')]).tolist()
    
    #Run pre-processing steps per sample
    for sample in tqdm.tqdm(bulkRNAseqSampleNames) :
        #Trim adapter sequences
        doAdapterTrim(pathToFastqFiles, sample, pathToRefAdapters, nRAM)
        #Trim reads based on quality score and minimum length
        doQualityTrim(pathToFastqFiles, sample, trimQuality, minLen, nRAM)
        #Remove phage spikein reads
        doPhageTrim(pathToFastqFiles, sample, pathToRefPhix, nRAM)
    
    #Create folder to store processed files
    pathToProcessed = os.path.join(pathToFastqFiles, 'processed')
    if not os.path.exists(pathToProcessed):
        os.makedirs(pathToProcessed)
    
    #Move processed files to a new folder and delete the temporary folder
    command = "mv {}/*_processed.fastq {}; rm -r {}; ".format (
        pathToTemp, pathToProcessed, pathToTemp)
    subprocess.check_call(command, shell=True)
    
    return pathToProcessed
    
### BULK PREPROCESSING FUNCTIONS

def doAdapterTrim(pathToFastqFiles, sample, pathToRefAdapters, nRAM ):
    
    command = "bbduk.sh -Xmx{}g in1={} in2={} out1={} out2={} \
        ktrim=r k=21 hdist=1 ref={}".format (
        nRAM,
        '/'.join([pathToFastqFiles, sample+'_1.fastq.gz']), 
        '/'.join([pathToFastqFiles, sample+'_2.fastq.gz']),
        '/'.join([pathToFastqFiles, 'temp', sample+'_1_trimmed.fastq']),
        '/'.join([pathToFastqFiles, 'temp', sample+'_2_trimmed.fastq']),
        pathToRefAdapters
        )
    
    subprocess.check_call(command, shell=True)
    
    
def doQualityTrim(pathToFastqFiles, sample, trimQuality, minLen, nRAM):
    
    command = "bbduk.sh -Xmx{}g in1={} in2={} out1={} out2={} trimq={} minlen={}".format(
        nRAM,
        '/'.join([pathToFastqFiles, 'temp', sample+'_1_trimmed.fastq']),
        '/'.join([pathToFastqFiles, 'temp', sample+'_2_trimmed.fastq']),
        '/'.join([pathToFastqFiles, 'temp', sample+'_1_quality.fastq']),
        '/'.join([pathToFastqFiles, 'temp', sample+'_2_quality.fastq']),
        trimQuality,
        minLen)
    
    subprocess.check_call(command, shell=True)
    
    
def doPhageTrim(pathToFastqFiles, sample, pathToRefPhix, nRAM):
    command = "bbduk.sh -Xmx{}g in1={} in2={} out1={} out2={} \
               outm1={} outm2={} ref={} k=31 hdist=1 stats={}".format(
        nRAM,
        '/'.join([pathToFastqFiles, 'temp', sample+'_1_quality.fastq']),
        '/'.join([pathToFastqFiles, 'temp', sample+'_2_quality.fastq']),
        '/'.join([pathToFastqFiles, 'temp', sample+'_1_processed.fastq']),
        '/'.join([pathToFastqFiles, 'temp', sample+'_2_processed.fastq']),
        '/'.join([pathToFastqFiles, 'temp', sample+'_1_matched.fastq']),
        '/'.join([pathToFastqFiles, 'temp', sample+'_2_matched.fastq']),
        pathToRefPhix,
        '/'.join([pathToFastqFiles, 'temp', 'stats.txt'])
        )
    subprocess.check_call(command, shell=True)

