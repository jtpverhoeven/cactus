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
from director import directorObject
import argparse

def cloudBlastInit():

    parser = argparse.ArgumentParser(description="EC2 cloud blast standalone")

    parser.add_argument('--fasta', dest='fasta', action='store',
                    help='Input fasta file', default="input.fastsa")

    parser.add_argument('--cfg', dest='cfg', action='store',
                    help='Configuration file containing ec2 and analysis parameters', default="ec2.cfg")

    args = parser.parse_args()

    #config file valid?
    if not os.path.isfile(args.cfg):
        logger.error('Configuration file was not set')
        raise ValueError('Could not find configuration file')

    config = configparser.RawConfigParser(allow_no_value=True)
    config.read(args.cfg)

    director = directorObject()
    director.connect(config.get('instances', 'instance_0'), config.getint('settings', 'port'))

    director.sendCommand('ping')

    director.sendFile(args.cfg, '/home/ec2-user/ec2received.cfg')
    director.sendCommand('load_config|/home/ec2-user/ec2received.cfg|True')

    #give a secon dto let the config file propagate onto the EFS
    time.sleep(1)

    #clean up working directory
    director.sendCommand('clean')
    time.sleep(1)

    #send and split fasta file
    director.sendFile(args.fasta, '/home/ec2-user/sequencesreceived.fasta')
    director.sendCommand('split_fasta|/home/ec2-user/sequencesreceived.fasta')

    director.sendCommand('bye')
    director.close()

    for instanceId,instanceIp in config.items('instances'):
        host = instanceIp
        instanceNumber = instanceId.split('_')
        port = config.getint('settings', 'port')

        ec2path = config.get('settings', 'cfg') + 'ec2received.cfg'
        director.connect(host, port)
        director.sendCommand('load_config|' + str(ec2path))
        director.sendCommand('start_instance_blast|' + str(instanceNumber[1]))
        director.sendCommand('bye')
        director.close()

    cloudReady = False
    while cloudReady == False:
        ec2path = config.get('settings', 'cfg') + 'ec2received.cfg'
        director.connect(config.get('instances', 'instance_0'), config.getint('settings', 'port'))
        director.sendCommand('load_config|' + str(ec2path))
        poll = director.sendBirectionalCommand('is_instance_ready')
        director.sendCommand('bye')
        director.close()

        if poll == '1':
            cloudReady = True
        else:
            time.sleep(30)

    director.connect(config.get('instances', 'instance_0'), config.getint('settings', 'port'))
    ec2path = config.get('settings', 'cfg') + 'ec2received.cfg'
    director.sendCommand('load_config|' + str(ec2path))
    dataObject = director.receiveResults()
    director.sendCommand('bye')
    director.close()

    pilePath = 'cloudblast_result.txt'
    pileFaPath = 'cloudblast_result.fa'

    with open(pileFaPath, 'w') as resultreceivedFasta, open(pilePath, 'w') as resultreceivedTxt:
            resultreceivedFasta.write(dataObject['resultsFa'])
            resultreceivedTxt.write(dataObject['results'])
    resultreceivedTxt.close()
    resultreceivedFasta.close()

if __name__ == "__main__":
    cloudBlastInit()
