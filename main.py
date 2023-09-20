#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 20 11:19:40 2023

@author: rafaella
"""

import os
import subprocess


#FASTQC
pathToFastqFiles = 'data/rnaseq/Sample_MUC30433'
nCores = 40 # max 40

pathToOutputDir = './results/fastqc'
if not os.path.exists(pathToOutputDir):
    os.makedirs(pathToOutputDir)

command = "fastqc -o {} -t {} {}".format(pathToOutputDir, str(nCores), pathToFastqFiles + '/*.gz')
subprocess.check_call(command, shell=True)


### Single cell

#CELLRANGER 

localMem = 200
expectCells = 1000
pathToRefTranscriptome = '../software/refdata-gex-GRCh38-2020-A'

sampleName = [file for file in os.listdir (pathToFastqFiles) if file.endswith('.gz')][0].split('_')[0]
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
                                                 localMem,
                                                 pathToFastqFiles
                            )
subprocess.check_call(command, shell=True)


#### bulk
pathToFastqFiles = 'data/rnaseq/bulk_RNAseq'

#TRIM ADAPTERS

pathToRefAdapters = 'data/resources/adapters.fa'

if not os.path.exists(os.path.join(pathToFastqFiles, 'temp')):
    os.makedirs(os.path.join(pathToFastqFiles, 'temp'))

bulkRNAseqFastqFiles = [file for file in os.listdir (pathToFastqFiles) if file.endswith('.gz')]
command = "bbduk.sh in1={} in2={} out1={} out2={} ref={}".format (
    '/'.join([pathToFastqFiles, bulkRNAseqFastqFiles[0]]), 
    '/'.join([pathToFastqFiles, bulkRNAseqFastqFiles[1]]),
    '/'.join([pathToFastqFiles, 'temp', bulkRNAseqFastqFiles[0].split('.')[0]+'_trimmed.fastq']),
    '/'.join([pathToFastqFiles, 'temp', bulkRNAseqFastqFiles[1].split('.')[0]+'_trimmed.fastq']),
    pathToRefAdapters
    )
subprocess.check_call(command, shell=True)

#QUALITY TRIM

command = "bbduk.sh -Xmx1g in1={} in2={} out1={} out2={} trimq=10".format(
    '/'.join([pathToFastqFiles, 'temp', bulkRNAseqFastqFiles[0].split('.')[0]+'_trimmed.fastq']),
    '/'.join([pathToFastqFiles, 'temp', bulkRNAseqFastqFiles[1].split('.')[0]+'_trimmed.fastq']),
    '/'.join([pathToFastqFiles, 'temp', bulkRNAseqFastqFiles[0].split('.')[0]+'_quality.fastq']),
    '/'.join([pathToFastqFiles, 'temp', bulkRNAseqFastqFiles[1].split('.')[0]+'_quality.fastq'])
    )
subprocess.check_call(command, shell=True)

#REMOVE PHIX SPIKEIN
pathToRefPhix = 'data/resources/phix174_ill.ref.fa.gz'

command = "bbduk.sh -Xmx1g in1={} in2={} out1={} out2={} \
           outm1={} outm2={} ref={} stats={}".format(
    '/'.join([pathToFastqFiles, 'temp', bulkRNAseqFastqFiles[0].split('.')[0]+'_quality.fastq']),
    '/'.join([pathToFastqFiles, 'temp', bulkRNAseqFastqFiles[1].split('.')[0]+'_quality.fastq']),
    '/'.join([pathToFastqFiles, 'temp', bulkRNAseqFastqFiles[0].split('.')[0]+'_unmatched.fastq']),
    '/'.join([pathToFastqFiles, 'temp', bulkRNAseqFastqFiles[1].split('.')[0]+'_unmatched.fastq']),
    '/'.join([pathToFastqFiles, 'temp', bulkRNAseqFastqFiles[0].split('.')[0]+'_matched.fastq']),
    '/'.join([pathToFastqFiles, 'temp', bulkRNAseqFastqFiles[1].split('.')[0]+'_matched.fastq']),
    pathToRefPhix,
    '/'.join([pathToFastqFiles, 'temp', 'stats.txt'])
    )
subprocess.check_call(command, shell=True)
 