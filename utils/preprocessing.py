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

    Returns
    -------


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

    #Move to results folder
    command ="mv {} {}".format(sampleName, pathToOutputDir)
    

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
@ERR862690.15583333 HS21_14478:5:1214:19704:JIGGIIIJJJIIIIIJJIIIJHDBBACCDCCC9ACDDDDDDDDDDD:<C@CADDEE

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
def extractCountMatrix (pathToAlignment):
    '''
    

    Parameters
    ----------
    pathToAlignment : TYPE
        DESCRIPTION.

    Returns
    -------
    countMatrixDF : TYPE
        DESCRIPTION.

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