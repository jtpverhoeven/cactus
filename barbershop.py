#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import sys
import re
import matplotlib.pyplot as plt
from cactusUtils import *
from collections import Counter
from shutil import copyfile
from subprocess import call
from collections import defaultdict
import logging


""" barbershop.py: Simple trimming on raw fastq reads """

class barberShopObject:

    def __init__(self, readGroups, outputPath, args):
        self.args = args
        self.readGroups = readGroups
        self.outputPath = outputPath
        self.barberOutput = outputPath + self.args.get('directories', 'barbershopOutput')
        checkDirOrCreate(self.barberOutput)

        self.barberUnivec =  outputPath + self.args.get('directories', 'barbershopOutput') + ('/univec')
        self.barberGraphs =  outputPath + self.args.get('directories', 'barbershopOutput') + ('/graphs')

        checkDirOrCreate(self.barberUnivec)
        checkDirOrCreate(self.barberGraphs)

        self.barberUnivecFasta =  outputPath + self.args.get('directories', 'barbershopOutput') + ('/univec/univec.fasta')
        self.barberUnivecDb =  outputPath + self.args.get('directories', 'barbershopOutput') + ('/univec/cactusunivec')

        self.logger = logging.getLogger('cactus')
        self.logger.info('Barbershop object created ')

        if self.args.getboolean('barbershop', 'enableUnivec') is True:
            self.checkAndCreateCustomUnivec()

    """ Run QC """
    def performQC(self):

        self.logger.info('Starting QC')

        for readgroup in self.readGroups:
            readgroupName = readgroup[0]
            readgroupFile = readgroup[1]

            self.logger.info('Processing ' + str(readgroupName))

            #check if this file exists
            if os.path.isfile(readgroupFile) is False:
                self.logger.error('Could not find sequence file as defined in readgroups:' + readgroupFile)
                raise Exception('Could not find sequence file as defined in readgroups:' + readgroupFile)

            #copy to staging area
            outputPath = self.barberOutput + '/' + readgroupName + '.original.fastq'
            copyfile(readgroupFile, outputPath)

            workFile = outputPath

            if self.args.getboolean('barbershop', 'enabled') is True:

                if self.args.getboolean('barbershop', 'enableTrimmer') is True:
                    self.logger.debug('Trimmer')
                    self.qcChecker(readgroupName,workFile)
                    workFile =  self.barberOutput + '/' + readgroupName + '.trimmed.fastq'

                if self.args.getboolean('barbershop', 'enableCutAdapt') is True:
                    self.logger.debug('CutAdapt screening')
                    self.primerTrimming(readgroupName,workFile)
                    workFile = self.barberOutput + '/' + readgroupName + '.cutadapt.fastq'

                if self.args.getboolean('barbershop', 'enableUnivec') is True:
                    self.logger.debug('Univec screening')
                    self.uniVecScreening(readgroupName,workFile)
                    workFile = self.barberOutput + '/' + readgroupName + '.univec.fastq'

                if self.args.getboolean('barbershop', 'enableCompression') is True:
                    self.logger.debug('Compressing')
                    self.compressReadgroup(readgroupName,workFile)
                    workFile = self.barberOutput + '/' + readgroupName + '.compressed.fastq'
            else:
                self.logger.warn('barbershop not enabled, skipping for:' + str(readgroupName))

            outputPath = self.barberOutput + '/' + readgroupName + '.fastq'
            copyfile(workFile, outputPath)

        return

    def compressReadgroup(self, readgroupName, readgroupInput):

        nullHandle = open('/dev/null', 'w')

        fastqOut = self.barberOutput + '/' + readgroupName + '.compressed.fastq'
        vsearchOut = self.barberOutput + '/' + readgroupName + '.compression.result'

        call([self.args.get('bin', 'vsearch'),
            '--quiet',
            '--cluster_fast', readgroupInput,
            '--consout', vsearchOut,
            '--sizeout',
            '--threads', '32',
            '--id', self.args.get('barbershop', 'compressionId'),
        ], stdout=nullHandle)


        centroids = defaultdict(dict)
        for name, seq in fastaIter(vsearchOut):
            realName = re.findall ( '\=(.*?)\;', name)
            try:
                centroids[realName[0]] = name
            except IndexError:
                print(name)
                print(realName)
                sys.exit()

        outputHandle = open(fastqOut, 'w')
        readgroupHandle = open(readgroupInput, 'r')

        self.logger.debug('Compressor number of centroids: ' + str(len(centroids)))

        for name, seq, qual in readfq(readgroupHandle):
            if name in centroids:
                outputHandle.write('@' + name + '\n')
                outputHandle.write(seq + '\n')
                outputHandle.write('+' + '\n')
                outputHandle.write(qual + '\n')

        nullHandle.close()
        outputHandle.close()
        return

    def checkAndCreateCustomUnivec(self):

        addition = self.args.get('barbershop', 'customUnivecAddition')
        original = self.args.get('barbershop', 'univecdb')

        if addition == 'False':
            self.logger.info('Creating univec database as-is')
            originalHandle = open(original, 'r')
            destinationHandle = open(self.barberUnivecFasta,'w')
            shutil.copyfileobj(originalHandle , destinationHandle, 1024*1024*10)
        else:
            univecDb = self.args.get('barbershop', 'univecDb')
            self.logger.info('Creating extended univec database')

            originalHandle = open(original, 'r')
            destinationHandle = open(self.barberUnivecFasta,'w')
            shutil.copyfileobj(originalHandle , destinationHandle, 1024*1024*10)
            destinationHandle.close()

            with open(self.barberUnivecFasta,'a') as customUnivec, open(addition, 'r') as univecAdditions:
                customUnivec.write(univecAdditions.read())

            customUnivec.close()
            univecAdditions.close()

        #else:

        call([self.args.get('bin', 'makeblastdb'),
            '-in', self.barberUnivecFasta,
            '-dbtype', 'nucl',
            '-input_type', 'fasta',
            '-title', 'cactusunivec',
            '-out', self.barberUnivecDb,
        ])


    def primerTrimming(self, readgroupName, readgroupInput):

        thisOut = self.barberOutput + '/' + readgroupName + '.cutadapt.fastq'
        log = self.barberOutput + '/' + readgroupName + '.cutadapt.report.txt'
        logAhandle = open(log, 'w')

        call([self.args.get('bin', 'cutadapt'),
                '-b', 'file:' + self.args.get('barbershop', 'cutAdaptPrimers'),
                '-o', thisOut,
                '-e', '0.1',
                readgroupInput,
            ], stdout=logAhandle)

        return


    def uniVecScreening(self, readgroupName, readgroupInput ):

        nullHandle = open('/dev/null', 'w')
        fastaOut = self.barberOutput + '/' + readgroupName + '.fasta'
        univecOut = self.barberOutput + '/' + readgroupName + '.txt'

        call([self.args.get('bin', 'vsearch'),
            '--quiet',
            '--fastx_filter', readgroupInput,
            '--fastq_qmax', '80',
            '--fastaout', fastaOut,
        ], stdout=nullHandle)

        call([self.args.get('bin', 'blastnbin'),
            '-task', 'blastn',
            '-reward', '1',
            '-penalty', '-4',
            '-gapopen', '3',
            '-gapextend', '3',
            '-dust', 'yes',
            '-soft_masking', 'true',
            '-evalue', '700',
            '-searchsp', '1750000000000',
            '-db', self.barberUnivecDb,
            '-query', fastaOut,
            '-num_threads', '32',
            '-out', univecOut,
            '-outfmt', '10 qseqid sseqid qstart qend length pident sstrand score',
            '-max_target_seqs', '1'
        ], stdout=nullHandle)

        #parse univec results
        univecHits = defaultdict(dict)
        with open(univecOut) as univecResults:
            for i, line in enumerate(univecResults):
                lineSplit = line.split(',')
                univecHits[lineSplit[0]] = defaultdict(dict)
                univecHits[lineSplit[0]]['begin'] = lineSplit[2]
                univecHits[lineSplit[0]]['end'] = lineSplit[3]

        outputPath = self.barberOutput + '/' + readgroupName + '.univec.fastq'
        rejectPath = self.barberOutput + '/' + readgroupName + '.rejected.fastq'
        rejectHandle = open(rejectPath, 'a')

        outputHandle = open(outputPath, 'w')
        inputHandle = open(readgroupInput, 'r')
        trimCount = 0
        trimDiscarded = 0

        for name, seq, qual in readfq(inputHandle):

            #univec found?
            if name in univecHits:
                endPoint = int(univecHits[name]['begin'])
                fromBeginning = int(univecHits[name]['begin'])
                fromEnd = len(seq) - int(univecHits[name]['end'])

                if fromBeginning <= 25 or fromEnd <= 25:
                    seq = seq[:endPoint]
                    qual = seq[:endPoint]
                    trimCount += 1

            thisSeqLen = len(seq)

            if thisSeqLen <= self.args.getint('barbershop', 'minLengthForUnivecTrim'):
                rejectHandle.write('@' + name + ' UNIVEC-FRAG-TOO-SMALL \n')
                rejectHandle.write(seq + '\n')
                rejectHandle.write('+' + '\n')
                rejectHandle.write(qual + '\n')
                trimDiscarded += 1
            else:
                outputHandle.write('@' + name + '\n')
                outputHandle.write(seq + '\n')
                outputHandle.write('+' + '\n')
                outputHandle.write(qual + '\n')

        self.logger.debug('Univec based trims: ' + str(trimCount))
        self.logger.debug('Univec based rejects: ' + str(trimDiscarded))

        outputHandle.close()
        inputHandle.close()
        return

    def qcChecker(self, readgroupName, readgroupInput):
        outputPath = self.barberOutput + '/' + readgroupName + '.trimmed.fastq'
        rejectPath = self.barberOutput + '/' + readgroupName + '.rejected.fastq'

        perBaseQuality = []
        perQualityLength = []

        inputHandle = open(readgroupInput, 'r')
        outputHandle = open(outputPath, 'w')
        rejectHandle = open(rejectPath, 'w')
        qualFile = open(self.barberOutput + '/' + readgroupName + '.quals.txt', 'w')

        readsDropped = 0
        readsCompressionRatioDropped = 0
        readsSalvaged = 0
        readsWritten = 0
        #n, slen, qlen = 0, 0, 0
        for name, seq, qual in readfq(inputHandle):


            if self.args.getboolean('barbershop', 'checkCompressionRatio') is True:
                thisCompressionRatio = compressionRatio(seq)
                if ( thisCompressionRatio <= self.args.getfloat('barbershop', 'compressionRatioCutOff')):
                    rejectHandle.write('@' + name + 'TRIMMER-COMPRESSION-RATIO (' + str(thisCompressionRatio) +')\n')
                    rejectHandle.write(seq + '\n')
                    rejectHandle.write('+' + '\n')
                    rejectHandle.write(qual + '\n')
                    readsDropped = readsDropped + 1
                    readsCompressionRatioDropped += 1
                    continue

            #length only
            if self.args.getboolean('barbershop', 'trimOnlyLength') is True:
                if ( len(seq) >=  self.args.getint('barbershop', 'minLength')):
                    outputHandle.write('@' + name + '\n')
                    outputHandle.write(seq + '\n')
                    outputHandle.write('+' + '\n')
                    outputHandle.write(qual + '\n')
                    readsWritten += 1
                else:
                    rejectHandle.write('@' + name + ' TRIMMER-TOO-SHORT \n')
                    rejectHandle.write(seq + '\n')
                    rejectHandle.write('+' + '\n')
                    rejectHandle.write(qual + '\n')
                    readsDropped = readsDropped + 1
            else:

                #get quality score and calcualte average quality
                qualities = [ord(x) - 33 for x in qual.rstrip('\n')]
                average_quality = sum(qualities) / len(qualities)

                c = Counter(qualities)
                dipSum = sum(v for k, v in c.items() if k < self.args.getint('barbershop','dipCheckTreshold'))
                qualFile.write(str(average_quality) + '\t' + str(dipSum) + '\n')

                if (average_quality >= self.args.getint('barbershop', 'minAverageQuality') and
                dipSum < self.args.getint('barbershop', 'maxNumberOfDips') and
                len(seq) >=  self.args.getint('barbershop', 'minLength')):
                    outputHandle.write('@' + name + '\n')
                    outputHandle.write(seq + '\n')
                    outputHandle.write('+' + '\n')
                    outputHandle.write(qual + '\n')
                    readsWritten += 1
                else:

                    if self.args.getboolean('barbershop','salvageEnabled') and len(seq) > self.args.getint('barbershop','salvageLength'):
                        outputHandle.write('@' + name + '\n')
                        outputHandle.write(seq + '\n')
                        outputHandle.write('+' + '\n')
                        outputHandle.write(qual + '\n')
                        readsSalvaged += 1
                        readsWritten += 1
                    else:
                        rejectHandle.write('@' + name + ' TRIMMER-QUALITY-OR-LENGTH LEN: ' + str(len(seq)) + ' AVG_QUAL: ' + str(average_quality) + ' DIPSUM: ' + str(dipSum) + '  \n')
                        rejectHandle.write(seq + '\n')
                        rejectHandle.write('+' + '\n')
                        rejectHandle.write(qual + '\n')
                        readsDropped = readsDropped + 1

        self.logger.debug('Trimmer quality control dropped reads: ' + str(readsDropped))
        self.logger.debug('Trimmer quality control dropped because of compression ratio: ' + str(readsCompressionRatioDropped))
        self.logger.debug('Trimmer quality control length-salvaged reads: ' + str(readsSalvaged))
        self.logger.debug('Trimmer sum of reads written: ' + str(readsWritten))
        return
