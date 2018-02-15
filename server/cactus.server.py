#!/usr/bin/python34
# -*- coding: utf-8 -*-
import math
import socket
import time
import shutil
import os, tarfile
import threading
import struct
import json
import time
import sys
import configparser
import urllib.request
from shutil import copyfile
from collections import defaultdict
from collections import Counter
from itertools import groupby
from subprocess import call
from subprocess import Popen, PIPE

class ncbiUpdater:

    def __init__(self, maxVolumes, prefix, storage):
        self.maxVolumes = range(0, maxVolumes)
        self.volumes = defaultdict(dict)
        self.prefix = prefix
        self.storage = storage
        self.baseLink = "ftp://ftp.ncbi.nlm.nih.gov/blast/db/"
        self.createFileRegister()
        return

    def createFileRegister(self):
        for i in self.maxVolumes:
            fileNamePadded =  self.prefix + "." + str(i).zfill(2) + ".tar.gz"
            fileNamePaddedMd5 =  self.prefix + "." + str(i).zfill(2) + ".tar.gz.md5"
            self.volumes[fileNamePadded] = defaultdict(dict)
            self.volumes[fileNamePadded]['name'] = fileNamePadded
            self.volumes[fileNamePadded]['url'] = self.baseLink + fileNamePadded
            self.volumes[fileNamePadded]['url_md5'] =  self.baseLink +  fileNamePaddedMd5
            self.volumes[fileNamePadded]['savePath'] = self.storage + fileNamePadded
            self.volumes[fileNamePadded]['savePathMd5'] = self.storage + fileNamePaddedMd5
            self.volumes[fileNamePadded]['ready']  = False
            self.volumes[fileNamePadded]['md5Ready']  = False
            self.volumes[fileNamePadded]['checksum_status']  = False

    def md5Check():
        return

    def update(self):
        while 1:
            allReady, grabFile = self.registerCompleted()
            if allReady is True:
                break
            else:
                print('Grabbing file' + str(grabFile))
                self.grabFile(grabFile)
        print('done')
        return

    def registerCompleted(self):
        for thisFileName, thisFileDict in self.volumes.items():
            if thisFileDict['ready'] == False:
                return False, thisFileName
        return True, None


    def grabFile(self, fileName):

        try:
            if(self.volumes[fileName]['md5Ready']== False):
                if os.path.isfile(self.volumes[fileName]['savePathMd5']) == False:
                    urllib.request.urlretrieve(self.volumes[fileName]['url_md5'], self.volumes[fileName]['savePathMd5'])
                self.volumes[fileName]['md5Ready'] = True
                return
            if(self.volumes[fileName]['ready']== False):
                if os.path.isfile(self.volumes[fileName]['savePath']) == False:
                    urllib.request.urlretrieve(self.volumes[fileName]['url'], self.volumes[fileName]['savePath'])
                self.volumes[fileName]['ready'] = True
                return
        except Exception as Ex:
            template = "An exception of type {0} occurred. Arguments:\n{1!r}"
            message = template.format(type(Ex).__name__, Ex.args)
            print(message)
            sys.exit()
        return

class cactusServer:

    def __init__(self):
        #self.host = "localhost"
        self.host = socket.gethostname()
        self.port = 5001
        self.connections = defaultdict(dict)
        self.config = None
        self.allowedCommands = ['receive_file',
                                'load_config',
                                'is_instance_ready',
                                'split_fasta',
                                'bye',
                                'ping',
                                'start_instance_blast',
                                'return_results',
                                'clean']

    def start(self):
        print('* Started')
        threading.Thread(target = self.listen,args = ()).start()

    def listen(self):
        cactusSocket = socket.socket()
        cactusSocket.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
        cactusSocket.bind((self.host,self.port))
        cactusSocket.listen(5)

        while True:
            print('* Waiting for connections')
            clientSocket,addr = cactusSocket.accept()
            thisIp = str(addr[0])
            self.connections[thisIp] = clientSocket
            print('* Connection from:' + str(addr))
            threading.Thread(target = self.clientConnection,args =(thisIp,)).start()

    def clientConnection(self,addr):
        while True:
            if self.connections[addr]:
                data = self.recvall(self.connections[addr], 4)
            else:
                print('* Connection key did not exist')
                data = False

            if not data:
                break

            length, = struct.unpack('!I', data)
            self.receive(addr, length)

        print('* Client connection closed:' + str(addr))

    def receive(self, addr, length ):
        dataReceived = self.recvall(self.connections[addr], length)
        commandObject = json.loads(dataReceived.decode())
        commandList = commandObject['cmd'].split('|')
        command = commandList[0]

        print('Command: ' + command)

        if command not in self.allowedCommands:
            print('* Invalid command received')
            return

        if command == 'receive_file':
            destinationFileName = commandObject['destination_file']
            with open(destinationFileName, 'w') as f:
                f.write(commandObject['contents'])
                f.close()

            self.connections[addr].send("1".encode())

        if command == 'ping':
            print('pong')
            self.connections[addr].send("1".encode())

        if command == 'load_config':
            print('Loading in config:' + commandList[1])
            self.config = configparser.RawConfigParser(allow_no_value=True)
            self.config.read(commandList[1])

            if len(commandList) > 2:
                if commandList[2] == 'True':
                    print('* Copying Config to EC2')
                    pathSplit, fileSplit = os.path.split(commandList[1])
                    manifestCopyPath = self.config.get('settings', 'cfg') + fileSplit
                    with open(commandList[1],'r') as manifestOriginal, open(manifestCopyPath,'w') as manifestCopyPath:
                        shutil.copyfileobj(manifestOriginal, manifestCopyPath, 1024*1024*10)

                manifestOriginal.close()
                manifestCopyPath.close()
            self.connections[addr].send("1".encode())

        if command == 'split_fasta':
            fastaName = commandList[1]
            self.splitFastaOnEFS(fastaName)
            pathSplit, fileSplit = os.path.split(fastaName)
            sequencesCopyPath = self.config.get('settings', 'seqin') + fileSplit
            with open(commandList[1],'r') as seqOriginal, open(sequencesCopyPath,'w') as sequencesCopy:
                shutil.copyfileobj(seqOriginal, sequencesCopy, 1024*1024*10)

            seqOriginal.close()
            sequencesCopy.close()
            self.connections[addr].send("1".encode())
            return

        if command == 'return_results':

            pilePath, pileFaPath = self.resultPileup()

            with open(pilePath, 'r') as f:
                pileToSend = f.read()
                f.close()

            with open(pileFaPath, 'r') as f:
                pileFaToSend = f.read()
                f.close()


            cmdJSON = json.dumps({'cmd' : 'receive_results', 'resultsFa': pileFaToSend ,  'results': pileToSend})
            cmdJsonBytes = bytes(cmdJSON, 'UTF-8')
            length = len(cmdJsonBytes)
            self.connections[addr].sendall(struct.pack('!I', length))
            self.connections[addr].sendall(cmdJsonBytes)

        if command == 'start_instance_blast':
            instanceId = commandList[1]
            threading.Thread(target=self.runBlast,args =(instanceId,)).start()
            self.connections[addr].send("1".encode())

        if command == 'is_instance_ready':
            allReady = True

            for instanceId,instanceIp in self.config.items('instances'):
                instanceNumber = instanceId.split('_')
                readyFile = self.config.get('settings', 'result') + str('ready') + str(instanceNumber[1])
                if os.path.isfile(readyFile) == False:
                    allReady = False
            if allReady:
                self.connections[addr].send("1".encode())
            else:
                self.connections[addr].send("0".encode())


        if command == 'clean':

            foldersToClean = [  self.config.get('settings', 'seqin'),
                                self.config.get('settings', 'seqsplit'),
                                self.config.get('settings', 'seqsplitresult'),
                                self.config.get('settings', 'result'),
                                ]

            for cleanFolder in foldersToClean:
                for workFile in os.listdir(cleanFolder):

                    filePath = os.path.join(cleanFolder, workFile)
                    try:
                        if os.path.isfile(filePath):
                            os.unlink(filePath)
                    except Exception as e:
                        print(e)

            self.connections[addr].send("1".encode())

        if command == 'bye':
            self.connections[addr].send("1".encode())
            del self.connections[addr]

    def resultPileup(self):
        pilePath = self.config.get('settings', 'result')  + str('result.txt')
        pileFaPath = self.config.get('settings', 'result') + str('result.fasta')

        with open(pilePath,'wb') as wfd, open(pileFaPath,'wb') as wfad :

            for instanceId,instanceIp in self.config.items('instances'):
                instanceNumber = instanceId.split('_')
                thisFaPath = self.config.get('settings', 'seqsplit') + 'sequences_instance_' + str(instanceNumber[1])  + str('.fasta')
                thisTxtPath = self.config.get('settings', 'seqsplitresult') + 'sequences_instance_'  +  str(instanceNumber[1])  + str('.txt')

                with open(thisFaPath,'rb') as fdFa, open(thisTxtPath,'rb') as fdTxt:
                    shutil.copyfileobj(fdTxt, wfd, 1024*1024*10)
                    shutil.copyfileobj(fdFa, wfad, 1024*1024*10)
                    fdTxt.close()
                    fdFa.close()
        wfd.close()
        wfad.close()
        return [pilePath, pileFaPath]

    def runBlast(self, instanceId):
        print('* Starting blast for instance:' + str(instanceId))
        profile = self.config.get('settings', 'profile')
        call([self.config.get('settings', 'blastnbin'),
                        '-query',  self.config.get('settings', 'seqsplit') + 'sequences_instance_' + str(instanceId)  + str('.fasta'),
                        '-db', self.config.get('settings', 'ntDb'),
                        '-out', self.config.get('settings', 'seqsplitresult') + 'sequences_instance_' + str(instanceId)  + str('.txt'),
                        '-num_threads' , self.config.get(profile, 'threads'),
                        '-task', 'blastn',
                        '-penalty',  self.config.get(profile, 'penalty'),
                        '-reward', self.config.get(profile, 'reward'),
                        '-gapopen', self.config.get(profile, 'gapopen'),
                        '-gapextend', self.config.get(profile, 'gapextend'),
                        '-evalue', self.config.get(profile, 'evalue'),
                        '-num_descriptions', self.config.get(profile, 'num_descriptions'),
                        '-num_alignments', self.config.get(profile, 'num_alignments'),
                        '-dust', 'yes',
                        '-max_hsps', self.config.get(profile, 'max_hsps'),
                        '-culling_limit', self.config.get(profile, 'culling_limit'),
                        '-soft_masking', 'true',
                        '-outfmt', '0'
                        ], stdout=PIPE,
                   bufsize=1, close_fds=True,universal_newlines=True)

        readyPath = self.config.get('settings', 'result') + '/ready' + str(instanceId)
        readyFile = open(readyPath, 'w')
        readyFile.write(str(time.time()))
        readyFile.close()

        return


    def splitFastaOnEFS(self, inputFasta):

        print('Splitting ')
        fasta = inputFasta
        nSeq = self.countFasta(fasta)

        nInstances = len(self.config.items('instances'))
        self.threadStore = self.config.get('settings', 'seqsplit')
        nSeqPerInstance = math.ceil(nSeq / nInstances)

        seqLib = defaultdict(dict)
        fiter = self.fastaIter(fasta)
        seqId = 0

        for ff in fiter:
            headerStr, seq = ff
            seqLib[seqId] = dict()
            seqLib[seqId]['name'] = headerStr
            seqLib[seqId]['seq'] = seq
            seqId += 1

        threadFiles = defaultdict(dict)
        threadFilesTxt = defaultdict(dict)
        for n in range(nInstances):
            print(self.threadStore + 'sequences_instance_' + str(n)  + str('.fasta'))
            threadFiles[n] = open(self.threadStore + 'sequences_instance_' + str(n)  + str('.fasta'), "w")

        seqListSorted = sorted(seqLib, key=lambda x: (len(seqLib[x]['seq'])), reverse=True)

        currThreadFile = 0
        for seqKeyId in seqListSorted:
            seqName = seqLib[seqKeyId]['name']
            seqData = seqLib[seqKeyId]['seq']

            threadFiles[currThreadFile].write('>' + seqName + '\n')
            threadFiles[currThreadFile].write( seqData + '\n')
            currThreadFile += 1

            if(currThreadFile == nInstances):
                currThreadFile = 0

        for n in range(nInstances):
            threadFiles[n].close()

        return

    def countFasta(self, fasta):
        nSeq = 0
        fiter = self.fastaIter(fasta)
        for ff in fiter:
            nSeq += 1
        return nSeq

    def fastaIter(self, fastaName):
        fh = open(fastaName)
        faiter = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))
        for header in faiter:
            header = header.__next__()[1:].strip()
            seq = "".join(s.strip() for s in faiter.__next__())
            yield header, seq

    def recvall(self, sock, count):
        buf = b''
        bufText = ''
        while count:
            newbuf = sock.recv(count)
            if not newbuf: return None
            buf += newbuf
            count -= len(newbuf)
        return buf

def main():
    #pull in DB
    #change S3 bucket name if needed
    call(["aws", "s3", "sync", "s3://cactusdb/nt/", "/home/ec2-user/scratch/"])
    serverObj = cactusServer()
    serverObj.start()


if __name__ == '__main__':
    main()
