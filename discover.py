#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import sys
from subprocess import call
from subprocess import Popen, PIPE
import pprint
import hashlib
import json

import configparser
import socket
from cactusUtils import *
from collections import Counter
from collections import defaultdict
import argparse
import re
import sys
import csv
import time
import glob,os
import operator
import math
import numpy as np
from shutil import copyfile
import logging
from director import directorObject

""" discover.py: Handles blast/assignment functions
"""

class discoveryObject:

    def __init__(self, args, outputPath, assemblyFasta):
        self.args = args
        self.ec2args = None
        self.outputPath = outputPath
        self.assemblyFasta = assemblyFasta
        self.discoveryOutput = outputPath + self.args.get('directories', 'discoveryWorkDir')
        checkDirOrCreate(self.discoveryOutput)

        self.incoming = self.discoveryOutput + '/incoming'
        self.blastnDir = self.discoveryOutput + '/blastn'
        self.diamondDir = self.discoveryOutput + '/diamond'
        checkDirOrCreate(self.incoming)
        checkDirOrCreate(self.blastnDir)
        checkDirOrCreate(self.diamondDir)

        self.threadStore = self.blastnDir + self.args.get('directories', 'threadStore')
        checkDirOrCreate(self.threadStore)

        self.assemblyCopy =   self.incoming + '/assembly.accept.fasta'
        self.assemblyReject =   self.incoming + '/assembly.reject.fasta'

        self.logger = logging.getLogger('cactus')
        self.logger.debug('Discovery object created ')

        return


    def discovery(self):

        self.incomingCheck()

        if self.args.getboolean('discovery', 'blastn') is True:

            if self.args.getboolean('discovery', 'ec2blast') is True:

                self.logger.info('Running EC2 blast')

                input("Please verify EC2 cloud config correct and instances are ready, press [enter] to continue or [ctrl] + [c] to quit")

                self.ec2args = configparser.RawConfigParser(allow_no_value=True)
                self.ec2args.read(self.args.get('discovery', 'ec2cfg'))

                ec2isup = False

                director = directorObject()

                while ec2isup == False:

                    #do a ping first to see if all the nodes are up and sync'd
                    allOk = True
                    for instanceId,instanceIp in self.ec2args.items('instances'):
                        host = instanceIp
                        instanceNumber = instanceId.split('_')
                        port = self.ec2args.getint('settings', 'port')

                        try:
                            director.connect(host, port)
                            director.sendCommand('ping')
                            director.sendCommand('bye')
                            director.close()
                            self.logger.info('Instance: ' + instanceId + ' is up! ')
                        except socket.timeout:
                            allOk = False
                            self.logger.info('Instance: ' + instanceId + ' was not up, timed out ')
                            break
                        except ConnectionRefusedError:
                            allOk = False;
                            self.logger.info('Instance: ' + instanceId + ' was not up, refused ')


                    if allOk is True:
                        ec2isup = True
                        break
                    else:
                        self.logger.info('Waiting 15 seconds to recheck EC2 cloud')
                        time.sleep(15)

                self.logger.info('EC2 cloud is up! Transmitting manifest and sequences')

                director.connect(self.ec2args.get('instances', 'instance_0'), self.ec2args.getint('settings', 'port'))
                director.sendCommand('ping')

                director.sendFile(self.args.get('discovery', 'ec2cfg'), '/home/ec2-user/ec2received.cfg')
                director.sendCommand('load_config|/home/ec2-user/ec2received.cfg|True')

                #give a secon dto let the config file propagate onto the EFS
                time.sleep(1)

                #clean up working directory
                director.sendCommand('clean')
                time.sleep(1)

                #send and split fasta file
                director.sendFile(self.assemblyCopy, '/home/ec2-user/sequencesreceived.fasta')
                director.sendCommand('split_fasta|/home/ec2-user/sequencesreceived.fasta')

                director.sendCommand('bye')
                director.close()

                self.logger.info('Files succesfully transmited, now starting BLAST on each instance')

                for instanceId,instanceIp in self.ec2args.items('instances'):
                    host = instanceIp
                    instanceNumber = instanceId.split('_')
                    port = self.ec2args.getint('settings', 'port')

                    ec2path = self.ec2args.get('settings', 'cfg') + 'ec2received.cfg'
                    director.connect(host, port)
                    director.sendCommand('load_config|' + str(ec2path))
                    director.sendCommand('start_instance_blast|' + str(instanceNumber[1]))
                    director.sendCommand('bye')
                    director.close()
                    self.logger.debug('Started:' + instanceId)

                cloudReady = False
                while cloudReady == False:
                    ec2path = self.ec2args.get('settings', 'cfg') + 'ec2received.cfg'
                    director.connect(self.ec2args.get('instances', 'instance_0'), self.ec2args.getint('settings', 'port'))
                    director.sendCommand('load_config|' + str(ec2path))
                    poll = director.sendBirectionalCommand('is_instance_ready')
                    director.sendCommand('bye')
                    director.close()

                    if poll == '1':
                        self.logger.info('EC2 cloud responded: Ready!')
                        cloudReady = True
                    else:
                        time.sleep(30)


                self.logger.info('Requesting results from instance 0')

                director.connect(self.ec2args.get('instances', 'instance_0'), self.ec2args.getint('settings', 'port'))
                ec2path = self.ec2args.get('settings', 'cfg') + 'ec2received.cfg'
                director.sendCommand('load_config|' + str(ec2path))
                dataObject = director.receiveResults()
                director.sendCommand('bye')
                director.close()

                self.logger.info('EC2 Results received, writing to file')

                pilePath = self.blastnDir +  '/' + str(self.ec2args.get('settings', 'profile')) + str('.txt')
                pileFaPath = self.blastnDir +  '/' + str(self.ec2args.get('settings', 'profile')) + str('.fasta')

                with open(pileFaPath, 'w') as resultreceivedFasta, open(pilePath, 'w') as resultreceivedTxt:
                        resultreceivedFasta.write(dataObject['resultsFa'])
                        resultreceivedTxt.write(dataObject['results'])
                resultreceivedTxt.close()
                resultreceivedFasta.close()

                self.logger.info('EC2-Blast done')

            if self.args.getboolean('discovery', 'ec2blast') is False:
                files = self.splitInput()
                profile = self.args.get('discovery', 'blastn_profile')
                ntBlastOutput = self.ntBlast(files, profile)

        if self.args.getboolean('discovery', 'diamond') is True:
            profile = self.args.get('discovery', 'diamond_profile')
            diamondOutput = self.diamond(profile)

        return

    def incomingCheck(self):
        fa = fastaIter(self.assemblyFasta)
        assemblyCopyHandle = open(self.assemblyCopy, 'w')
        assemblyRejectHandle = open(self.assemblyReject, 'w')
        count = 0
        for name,seq in fa:

            if self.args.getboolean('discovery', 'trimOutHighCompression') is True:
                bseq = bytes(seq, 'UTF-8')
                complexity = float(len(zlib.compress(bseq)))/len(bseq)
                if complexity < 0.15:
                    count += 1
                    assemblyRejectHandle.write('>' + name.strip() + '\t' + str(complexity)  + '\n')
                    assemblyRejectHandle.write(seq + '\n')
                else:
                    assemblyCopyHandle.write('>' + name.strip() + '\n')
                    assemblyCopyHandle.write(seq + '\n')
            else:
                assemblyCopyHandle.write('>' + name.strip() + '\t' + str(complexity)  + '\n')
                assemblyCopyHandle.write(seq + '\n')

        self.assemblyFasta = self.assemblyCopy
        self.logger.info('Dropped ' + str(count) + ' high compression contigs')

        assemblyCopyHandle.close()
        assemblyRejectHandle.close()
        return

    def diamond(self, profile):
        self.logger.info('Running DIAMOND with profile: ' + str(profile))
        diamondOut = self.diamondDir + '/' + str(profile) + '.daa'
        diamondOutTxt = self.diamondDir + '/' + str(profile) + '.txt'

        #self.assemblyCopy

        call([self.args.get('bin', 'diamond'), 'blastx',
            '--db', self.args.get(profile, 'db'),
            '--query', self.assemblyCopy,
            '--out', diamondOut, '--outfmt', '100',
            '--max-target-seqs', self.args.get(profile, 'maxtargetseqs'),
            '--evalue', self.args.get(profile, 'evalue'),
            '--block-size', self.args.get(profile, 'blocksize'),
            '--index-chunks', self.args.get(profile, 'indexchunks'),
            '--salltitles', '--sallseqid',
            self.args.get(profile, 'additional-flag')
        ])

        call([self.args.get('bin', 'diamond'), 'view',
            '--daa', diamondOut,
            '--outfmt', '0',
            '--out', diamondOutTxt,
        ])

    def ntBlast(self, files, profile):

        self.logger.info('Running BLASTn with profile: ' + str(profile))

        threadShape = self.args.get('discovery', 'threadShape')
        if threadShape == 'linear':
            nThreadsPerWorker = 1
            self.logger.debug('Threadshape linear')
        if threadShape == 'square':
            self.logger.debug('Threadshape square')
            nThreadsPerWorker = self.args.getint('discovery', 'threadShapeFactor')
            self.logger.debug('nThreadsPerWorker ' +str(nThreadsPerWorker))

        processes = [Popen([self.args.get('bin', 'blastnbin'),
                            '-query',  self.threadStore + '/results_thread_' + str(i)  + str('.fasta'),
                            '-db',  self.args.get(profile, 'db'),
                            '-out',self.threadStore + '/results_thread_' + str(i)  + str('.txt'),
                            '-num_threads' , str(nThreadsPerWorker),
                            '-task', 'blastn',
                            '-penalty',  self.args.get(profile, 'penalty'),
                            '-reward', self.args.get(profile, 'reward'),
                            '-gapopen', self.args.get(profile, 'gapopen'),
                            '-gapextend', self.args.get(profile, 'gapextend'),
                            '-evalue', self.args.get(profile, 'evalue'),
                            '-num_descriptions', self.args.get(profile, 'num_descriptions'),
                            '-num_alignments', self.args.get(profile, 'num_alignments'),
                            '-dust', 'yes',
                            '-max_hsps', self.args.get(profile, 'max_hsps'),
                            '-culling_limit', self.args.get(profile, 'culling_limit'),
                            '-soft_masking', 'true',
                            '-outfmt', '0'
                            ], stdout=PIPE,
                       bufsize=1, close_fds=True,
                       universal_newlines=True)
                 for i in range( int(len(files)))]

        #print( str(len(files)) + ' * '  + self.args.get('discovery', 'threadShapeFactor')  + ' threads spun up.')

        self.logger.debug(str(len(files)) + ' * '  + self.args.get('discovery', 'threadShapeFactor')  + ' threads spun up.')
        tLeft = len(files)

        while processes:#
            for p in processes[:]:
                if p.poll() is not None: # process ended
                    tLeft = tLeft - 1
                    #print('Threads left' + str(tLeft))
                    self.logger.debug('Thread groups left: ' + str(tLeft))
                    processes.remove(p)

            time.sleep(1)

        self.resultPileup(profile, files)

    def resultPileup(self, profile, files):
        pilePath = self.blastnDir +  '/' + str(profile) + str('.txt')
        pileFaPath = self.blastnDir +  '/' + str(profile) + str('.fasta')
        threadShape = self.args.get('discovery', 'threadShape')

        if threadShape == 'linear':
            nThreads = self.args.getint('discovery', 'threads')
        if threadShape == 'square':
            nThreads = int(self.args.getint('discovery', 'threads') / self.args.getint('discovery', 'threadShapeFactor'))

        with open(pilePath,'wb') as wfd, open(pileFaPath,'wb') as wfad :
            for n in range(nThreads):
                thisFaPath = self.threadStore + '/results_thread_' + str(n)  + str('.fasta')
                thisTxtPath = self.threadStore + '/results_thread_' + str(n)  + str('.txt')

                with open(thisFaPath,'rb') as fdFa, open(thisTxtPath,'rb') as fdTxt:
                    shutil.copyfileobj(fdTxt, wfd, 1024*1024*10)
                    shutil.copyfileobj(fdFa, wfad, 1024*1024*10)
            #for fId, f in files.items():

    def splitInput(self):
        fasta = self.assemblyFasta
        nSeq = self.countFasta(fasta)
        nThreads = self.args.getint('discovery', 'threads')
        threadShape = self.args.get('discovery', 'threadShape')

        if threadShape == 'linear':
            nSeqPerCore = math.ceil(nSeq / nThreads)
        if threadShape == 'square':
            nThreads = int(self.args.getint('discovery', 'threads') / self.args.getint('discovery', 'threadShapeFactor'))
            nSeqPerCore = math.ceil(nSeq / nThreads)

        seqLib = defaultdict(dict)
        fiter = fastaIter(fasta)
        seqId = 0

        for ff in fiter:
            headerStr, seq = ff
            seqLib[seqId] = dict()
            seqLib[seqId]['name'] = headerStr
            seqLib[seqId]['seq'] = seq
            seqId += 1

        threadFiles = defaultdict(dict)
        threadFilesTxt = defaultdict(dict)
        for n in range(nThreads):
            threadFiles[n] = open(self.threadStore + '/results_thread_' + str(n)  + str('.fasta'), "w")
            threadFilesTxt[n]['query'] =   self.threadStore + '/results_thread_' + str(n) + str('.fasta')
            threadFilesTxt[n]['result'] =   self.threadStore + '/results_thread_' + str(n) + str('.txt')

        seqListSorted = sorted(seqLib, key=lambda x: (len(seqLib[x]['seq'])), reverse=True)

        currThreadFile = 0
        for seqKeyId in seqListSorted:
            seqName = seqLib[seqKeyId]['name']
            seqData = seqLib[seqKeyId]['seq']

            threadFiles[currThreadFile].write('>' + seqName + '\n')
            threadFiles[currThreadFile].write( seqData + '\n')
            currThreadFile += 1

            if(currThreadFile == nThreads):
                currThreadFile = 0

        for n in range(nThreads):
            threadFiles[n].close()

        return threadFilesTxt

    def countFasta(self, fasta):
        nSeq = 0
        fiter = fastaIter(fasta)
        for ff in fiter:
            nSeq += 1

        return nSeq
