#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import sys
import aligntest
from subprocess import call
import pprint
import hashlib
import json
import csv
import time
from cactusUtils import *
from contig import contigObject
from collections import Counter
from collections import defaultdict
from shutil import copyfile
from subprocess import Popen, PIPE
import matplotlib.pyplot as plt
import logging
import math
import numpy

""" spike.py: Phrap powered assembler
    This will fairly slow and unoptimized, but for
    virus discovery purposes this seems to perform the best
    at assembling from low-depth datasets
"""

class spikeObject:

    def __init__(self, readGroups, outputPath, args):
        self.args = args
        self.readGroups = readGroups
        self.outputPath = outputPath
        self.spikeOutput = outputPath + self.args.get('directories', 'spikeOutput')

        self.spikeWork = outputPath + self.args.get('directories', 'spikeOutput') + self.args.get('directories', 'spikeWork')
        self.spikeContigs = outputPath + self.args.get('directories', 'spikeOutput') + str('/contigs')
        self.spikeGraphs = outputPath + self.args.get('directories', 'spikeOutput') + str('/graphs')
        self.barbershopOutput = outputPath + self.args.get('directories', 'barbershopOutput')
        self.seqConnector = defaultdict(dict)
        self.poolInMemory = defaultdict(dict)
        self.seqCounter = 0
        self.currentPass = 0
        self.poolHandles = defaultdict(dict)
        self.poolContent = defaultdict(dict)

        self.contigRegister = defaultdict(dict)

        checkDirOrCreate(self.spikeOutput)
        checkDirOrCreate(self.spikeWork)
        checkDirOrCreate(self.spikeContigs)
        checkDirOrCreate(self.spikeGraphs)

        self.logger = logging.getLogger('cactus')
        self.logger.info('Spike object created ')


    def combineReadGroups(self):

        self.logger.info('Combining read groups')
        outputPath = self.spikeOutput + '/pool.fasta'
        outputHandle = open(outputPath, 'w')

        for readGroup in self.readGroups:
            readGroupName = readGroup[0]
            self.logger.debug('Processing: ' + readGroupName)
            readGroupInputPath = self.barbershopOutput  + '/' + readGroupName + '.fastq'
            inputHandle = open(readGroupInputPath, 'r')
            for name, seq, qual in readfq(inputHandle):
                m = hashlib.sha256()
                m.update(name.encode('utf-8'))
                hashedName = m.hexdigest()

                if hashedName in self.seqConnector:
                    self.logger.warning('Sequence already exists:'  + hashedName + '/' + name)

                outputHandle.write('@' + hashedName + '\n')
                outputHandle.write(seq + '\n')
                outputHandle.write('+' + '\n')
                outputHandle.write(qual + '\n')

                self.seqConnector[hashedName]['originalName'] = name
                self.seqConnector[hashedName]['readGroupOrigin'] = readGroupName
                self.seqConnector[hashedName]['originalSeq'] = seq

        self.logger.debug('Written : ' + str(len(self.seqConnector)) + ' reads')
        outputHandle.close()
        return

    """ Loads pool into memory so we dont have to perform streneous file operations constantly """
    def loadPoolInMemory(self, specificPool = None):

        self.poolInMemory = defaultdict(dict)
        self.seqCounter = 0

        if specificPool is None:
            poolPath = self.spikeOutput + '/pool.fasta'
        else:
            poolPath = self.spikeOutput + '/' + str(specificPool)

            self.logger.debug('Loading pool in memory: ' + str(poolPath))

        for name, seq, qual in fqIter(poolPath):
            if name in self.poolInMemory:
                raise Exception('I was about to overwrite a sequence! namely this one:' + str(name))
                sys.exit()

            self.poolInMemory[name]['name'] = name
            self.poolInMemory[name]['seq'] = seq
            self.poolInMemory[name]['len'] = len(seq)
            self.poolInMemory[name]['fastq'] = qual
            self.poolInMemory[name]['qual'] =  [ord(x) - 33 for x in qual.rstrip('\n')]
            self.seqCounter += 1

        return True


    def updateContigRegister(self, contigName, contents):

        if contigName in self.contigRegister:
            self.logger.warning('Contig name already present')
        else:
            self.contigRegister[contigName] = contents
        return

    def recreatePool(self, attrition = None):

        originalSequences = defaultdict(dict)

        if self.currentPass == 0:
            previousPoolN = '0'
            previousPool = self.spikeOutput + '/pool.fasta'
        else:
            previousPoolN = self.currentPass - 1
            previousPool = self.spikeOutput + '/pool_' + str(previousPoolN)  +'.fastq'

        self.logger.debug('Recreating pool, attrition: ' + str(attrition))
        self.logger.debug('Previous pool: ' + str(previousPool))

        previousPoolCount = 0

        for name, seq, qual in fqIter(previousPool):
            originalSequences[name]['name'] = name.rstrip('\n')
            originalSequences[name]['seq'] = seq.rstrip('\n')
            originalSequences[name]['fastq'] = qual.rstrip('\n')
            originalSequences[name]['qual'] =  [ord(x) - 33 for x in qual.rstrip('\n')]
            previousPoolCount += 1

        poolPath = self.spikeOutput + '/pool_' + str(self.currentPass) + '.fastq'
        poolHandle = open(poolPath, 'w', 8388608)
        writtenInPoolHandle = []

        if attrition is None:
            poolN = self.args.getint('spike', 'numberOfPools')
        else:
            poolN = attrition

        contigCounter = 0
        singleCounter = 0

        for i in range(0, poolN):
            contigLocation = self.spikeWork + '/pool_' + str(i) + '.fasta.contigs'
            contigQualLocation = self.spikeWork + '/pool_' + str(i) + '.fasta.contigs.qual'
            singletsLocation  = self.spikeWork + '/pool_' + str(i) + '.fasta.singlets'
            aceLocation = self.spikeWork + '/pool_' + str(i) + '.fasta.ace'
            contigCompositions = parseAce(aceLocation)

            #store quality values
            contigQualities = defaultdict(dict)

            for name, qual in qualIter(contigQualLocation):
                m = hashlib.sha256()
                m.update(name.encode('utf-8'))
                hashedName = m.hexdigest()
                contigQualities[hashedName] =  [chr(int(x) + 33)  for x in qual]

            for name, seq in fastaIter(contigLocation):
                m = hashlib.sha256()
                m.update(name.encode('utf-8'))
                hashedName = m.hexdigest()

                forwardName = ''.join([name, str(contigCounter), str(self.currentPass)])
                m.update(forwardName.encode('utf-8'))
                forwardHash = m.hexdigest()

                contigNameSplit = name.split('.')
                contigNumberName = contigNameSplit[-1]

                thisContigComposition = contigCompositions[contigNumberName]
                self.updateContigRegister(forwardHash,thisContigComposition)

                thisQuals =  "".join(str(x) for x in contigQualities[hashedName])
                poolToWrite =  "".join(['@', forwardHash, '\n', seq, '\n', '+',  '\n', thisQuals, '\n' ])
                poolHandle.write(poolToWrite)
                contigCounter += 1

            for name, seq in fastaIter(singletsLocation):

                newhasher = 'S' + str(name) + str(i) + str(self.currentPass)
                m = hashlib.sha256()
                m.update(newhasher.encode('utf-8'))
                hashedName = m.hexdigest()
                self.updateContigRegister(hashedName, [name])

                try:
                    poolToWrite = ''.join(['@', hashedName, '\n', seq, '\n', '+', '\n', originalSequences[name]['fastq'] , '\n'])
                    poolHandle.write(poolToWrite)
                except KeyError:
                    print('Could not find something for sequence: "' + name + '"')
                    sys.exit()

                singleCounter += 1


        self.logger.debug('Number of joins: ' + str(contigCounter))
        self.logger.debug('Number of singlets: ' + str(singleCounter))

        poolHandle.close()
        self.loadPoolInMemory('pool_' + str(self.currentPass) + '.fastq')
        return

    def debugFasta(self, name, seq):
        print('>' + str(name) + str('\n') + str(seq) + str('\n'))

    def doCodonCloneAlignment(self):
        self.combineReadGroups()
        self.loadPoolInMemory()
        self.codonCloneAssemble()
        assemblyFile = self.tidyUp()
        figureOutContigs = self.unravelContigs(assemblyFile)
        return assemblyFile

    def tidyUp(self):

        pool =  self.spikeOutput + '/pool_' + str(self.currentPass)  +'.fastq'
        self.logger.debug('Generating end result, from: ' + str(pool))
        result =  self.spikeOutput + '/spike_result.fasta'
        resultSorted =  self.spikeOutput + '/spike_result_sorted.fasta'
        nullHandle = open('/dev/null', 'w')
        lengths = []

        outputHandle = open(result, 'w')

        for name, seq, qual in fqIter(pool):
            outputHandle.write('>' + name + '\n')
            outputHandle.write(seq.strip() + '\n')
            lengths.append(len(seq))

        outputHandle.close()

        call([self.args.get('bin', 'vsearch'),
            '--quiet',
            '--sortbylength', result,
            '--output', resultSorted
        ], stdout=nullHandle)

        return resultSorted

    def unravelContigs(self, assemblyFile):

        contigFile = open(self.spikeOutput + '/contigregister.txt', 'w')
        for name, seq in fastaIter(assemblyFile):

            toFigureOut = self.contigRegister[name]
            contigFastaFile = self.spikeContigs + '/' + name + '.fasta'
            contigFasta = open(contigFastaFile, 'w')

            while len(toFigureOut) > 0:
                for sequence in toFigureOut:
                    if sequence in self.seqConnector:
                        #update register
                        lineToWrite = name + str('\t') + str(self.seqConnector[sequence]['originalName']) + str('\t') + str(self.seqConnector[sequence]['readGroupOrigin']) + str('\n')
                        contigFile.write(lineToWrite)
                        #update fasta
                        contigFasta.write('>' + str(self.seqConnector[sequence]['originalName']) + '\n')
                        contigFasta.write(str(self.seqConnector[sequence]['originalSeq']) + '\n')
                        toFigureOut.remove(sequence)
                    else:
                        toFigureOut = toFigureOut + self.contigRegister[sequence]
                        toFigureOut.remove(sequence)

            contigFasta.close()

        return


    def codonCloneAssemble(self):

        self.currentPass = 0
        attrition = 32
        self.logger.info('** Strict attrition:' + str(attrition))
        self.createPools(attrition)
        self.alignPools(attrition, True)
        self.recreatePool(attrition)
        self.flushWork()

        # #squeeze strict
        self.currentPass += 1
        attrition = 16
        self.logger.info('** Strict attrition:' + str(attrition))
        self.createPools(attrition)
        self.alignPools(attrition, True)
        self.recreatePool(attrition)
        self.flushWork()

        self.currentPass += 1
        attrition = 8
        self.logger.info('** Strict attrition:' + str(attrition))
        self.createPools(attrition)
        self.alignPools(attrition, True)
        self.recreatePool(attrition)
        self.flushWork()

        self.currentPass += 1
        attrition = 4
        self.logger.info('** Strict attrition:' + str(attrition))
        self.createPools(attrition)
        self.alignPools(attrition, True)
        self.recreatePool(attrition)
        self.flushWork()

        self.currentPass += 1
        attrition = 1
        self.logger.info('** Strict attrition:' + str(attrition))
        self.createPools(attrition)
        self.alignPools(attrition, True)
        self.recreatePool(attrition)
        self.flushWork()

        self.currentPass += 1
        attrition = 4
        self.logger.info('** Loose attrition:' + str(attrition))
        self.createPools(attrition)
        self.alignPools(attrition, None)
        self.recreatePool(attrition)
        self.flushWork()


        self.currentPass += 1
        attrition = 1
        self.logger.info('** Loose attrition:' + str(attrition))
        self.createPools(attrition)
        self.alignPools(attrition, None)
        self.recreatePool(attrition)
        self.flushWork()


        return

    def flushWork(self):

        for workFile in os.listdir(self.spikeWork):
            workFilePath = os.path.join(self.spikeWork, workFile)
            try:
                if os.path.isfile(workFilePath):
                    os.unlink(workFilePath)
                #elif os.path.isdir(file_path): shutil.rmtree(file_path)
            except Exception as e:
                print(e)

    def alignPools(self, attrition = None, broad = None):

        if attrition is None:
            poolN = self.args.getint('spike', 'numberOfPools')
        else:
            poolN = attrition

        self.logger.debug('Pool N:' + str(poolN))
        self.logger.debug('Attrition:' + str(attrition))
        self.logger.debug('Broad:' + str(broad))

        nullHandle = open('/dev/null', 'w')

        if broad == True:

            processes = [Popen([self.args.get('bin', 'phrap'),
                 self.spikeWork + '/pool_' + str(i) + '.fasta',
                '-penalty', '-3',
                '-bandwidth', '4',
                #'-minscore', '40',
                '-minscore', '36',
                #'-maxgap', '30',
                '-maxgap', '35',
                '-ace',
                #'-repeat_stringency', '.95',
                '-repeat_stringency', '.92',
                '-forcelevel', '1',
                '-gap_init', '-3',
                '-gap_ext', '-2',
                '-node_seg', '16',
                '-node_space', '12',
                '-minmatch' , '20',
                '-maxmatch', ' 28',
                '-indexwordsize', '10',
                #'-trim_qual', '20',
                #'-trim_score', '10'
                '-trim_qual', '2',
                '-trim_score', '2'
            ], stdout=nullHandle, stderr=nullHandle, bufsize=1, close_fds=True, universal_newlines=True)
            for i in range(0, poolN)]



        if broad == None:
            processes = [Popen([self.args.get('bin', 'phrap'),
                 self.spikeWork + '/pool_' + str(i) + '.fasta',
                '-penalty', '-3',
                '-bandwidth', '4',
                #'-minscore', '40',
                '-minscore', '30',
                #'-maxgap', '30',
                '-maxgap', '50',
                '-ace',
                #'-repeat_stringency', '.95',
                '-repeat_stringency', '.90',
                '-forcelevel', '1',
                '-gap_init', '-3',
                '-gap_ext', '-2',
                '-node_seg', '16',
                '-node_space', '12',
                '-minmatch' , '18',
                '-maxmatch', ' 26',
                '-indexwordsize', '10',
                #'-trim_qual', '20',
                #'-trim_score', '10'
                '-trim_qual', '2',
                '-trim_score', '2'
            ], stdout=nullHandle, stderr=nullHandle, bufsize=1, close_fds=True, universal_newlines=True)
            for i in range(0, poolN)]

        tLeft = poolN
        while processes:

            for p in processes[:]:
                if p.poll() is not None: # process ended
                    tLeft =  tLeft - 1
                    self.logger.debug('PHRAP threads left: ' + str(tLeft))
                    #print('Threads left' + str(tLeft))
                    processes.remove(p)

            time.sleep(1)

        nullHandle.close()
        time.sleep(1)
        return

    def createPools(self, attrition = None):

        written = []
        lengths = []
        if attrition is None:
            poolN = self.args.getint('spike', 'numberOfPools')
        else:
            poolN = attrition

        for i in range(0,poolN):
            handleLocation = self.spikeWork + '/pool_' + str(i) + '.fasta'
            qualLocation = self.spikeWork + '/pool_' + str(i) + '.fasta.qual'
            self.poolHandles[i]['fasta'] = open(handleLocation, 'w',8388608)
            self.poolHandles[i]['qual'] = open(qualLocation, 'w',8388608)

        currentPool = 0
        sortedList = sorted(self.poolInMemory, key=lambda x: (len(self.poolInMemory[x]['seq'])), reverse=True)
        numberOfSequences = len(sortedList)
        nSeqPerCore = math.floor(numberOfSequences / poolN)

        self.logger.debug('Length of list: ' + str(numberOfSequences) )

        #for sequenceName, sequenceInfo in sortedList.items():
        seqsInThisCore = 0
        seqsWritten = 0

        #some diagnosics
        seqsOver1kb = 0
        seqsOver1kbLengths = []
        seqsOver1KbTotal = 0

        for sequenceName in sortedList:

            sequenceInfo = self.poolInMemory[sequenceName]
            seqLen = len(sequenceInfo['seq'])

            if seqLen > 1000:
                seqsOver1kb += 1
                seqsOver1KbTotal += seqLen
                seqsOver1kbLengths.append(seqLen)

            self.poolHandles[currentPool]['fasta'].write('>' + sequenceName + '\n' + sequenceInfo['seq'] + '\n')
            seqsWritten += 1

            qualAsStr = ' '.join(map(str, sequenceInfo['qual']))
            self.poolHandles[currentPool]['qual'].write( ''.join(['>', sequenceName, '\n', qualAsStr , '\n']))
            seqsInThisCore += 1

            if seqsInThisCore == nSeqPerCore:
                currentPool = currentPool + 1
                if currentPool >= poolN:
                    currentPool -= 1
                else:
                    seqsInThisCore = 0


        self.logger.debug('Sequences written: ' + str(seqsWritten))
        self.logger.debug('Sequences over 1KB : ' + str(seqsOver1kb))
        self.logger.debug('Nucleotides in over 1KB fragments : ' + str(seqsOver1KbTotal))

        avgof1kb = seqsOver1KbTotal / (seqsOver1kb + 1)
        self.logger.debug('Average length of over 1KB fragments : ' + str(avgof1kb))

        if len(seqsOver1kbLengths) > 0:
           self.logger.debug('Min length of over 1KB fragments : ' + str(min(seqsOver1kbLengths)))
           self.logger.debug('Max length of over 1KB fragments : ' + str(max(seqsOver1kbLengths)))

        for poolIdx,pool in self.poolHandles.items():
            pool['fasta'].close()
            pool['qual'].close()
       
        return
