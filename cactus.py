#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
import sys
import configparser
import argparse
import logger
import barbershop
import pprint
from art import *
from cactusUtils import *
from spike import spikeObject
from discover import discoveryObject
import logging


"""cactus.py: Main entry point into the Cactus workflow"""

__author__ = "Joost Verhoeven"
__copyright__ = "Copyright 2017, Memorial University of Newfoundland"
__credits__ = ["Dr. Suzanne Dufour" , "Dr. Marta Canuti"]

__license__ = "GPL"
__version__ = "1.0"
__maintainer__ = "Joost Verhoeven"
__email__ = "jverhoeven@mun.ca"
__status__ = "Development"


""" The lobby boy receives and dispatches requests based upon the configuration file settings """
class lobbyBoy(object):
    def __init__(self, args, cargs, logger):
        self.args = args
        self.cargs = cargs
        self.readGroups = []
        self.basePath = os.path.dirname(os.path.realpath(__file__))
        self.outputPath = self.basePath +  self.args.get('directories', 'outputDirectory')
        self.logger = logger

        #run setup functions
        self.numberOfReadGroups = self.storeReadGroups()
        self.checkDirectoryOk()


    """ Is the output directory already there, if so quit  """
    def checkDirectoryOk(self):
        resultPathExisted = checkDirOrCreate(self.outputPath)
        if resultPathExisted is True:
            #if set to overwrite, drop directory drop directory and create again
            if self.args.getboolean('advanced', 'overwriteOnExist') is False:
                raise Exception('Output directory already eixsted')
            else:
                #we now also need to check if we are going to remove it or not based on
                #skip_qc and skip_assembly
                if self.cargs.skip_qc is False and self.cargs.skip_assembly is False and self.cargs.skip_discovery is False:
                    removeDirectoryTree(self.outputPath)
                    resultPathExisted = checkDirOrCreate(self.outputPath)
                    if resultPathExisted is True:
                        raise Exception('An error occured in overwriting the existing results folder')
                else:
                    self.logger.warn('Skipped directory removal because a skip step was found')

        return True

    """ Stores, and returns the number of, all indicated FASTQ files in self.readgroups """
    def storeReadGroups(self):
        self.readGroups = self.args.items('readgroups')
        self.logger.info('Read groups: ' + str(self.readGroups))
        return len(self.readGroups)

    """ Check if we need to run QC, otherwise, copy the unmodified files to the staging location """
    def runQC(self):
        bs = barbershop.barberShopObject(self.readGroups, self.outputPath, self.args)
        bs.performQC()
        return

    """ Run assembly  """
    def assembly(self):
        spike = spikeObject(self.readGroups, self.outputPath, self.args)
        assemblyFasta = spike.doCodonCloneAlignment()
        self.discovery(assemblyFasta)

    def discovery(self, assemblyFasta):
        if assemblyFasta is False:
            assemblyFasta = self.outputPath + self.args.get('directories', 'spikeOutput') + '/spike_result_sorted.fasta'

        disco = discoveryObject(self.args, self.outputPath,assemblyFasta)
        disco.discovery()

    def carbonCopyLogFile(self):
        logFileCopy = self.outputPath + '/runinfo.log'
        with open('runinfo.log','r') as logFile, open(logFileCopy,'w') as logFileCopyHandle:
                shutil.copyfileobj(logFile, logFileCopyHandle, 1024*1024*10)


def initLogger(config):

    logger = logging.getLogger('cactus')
    logger.setLevel(logging.DEBUG)

    ch = logging.StreamHandler()
    ch.setLevel(logging.DEBUG)

    fh = logging.FileHandler('runinfo.log', mode='w')
    fh.setLevel(logging.DEBUG)

    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    ch.setFormatter(formatter)
    fh.setFormatter(formatter)

    logger.addHandler(ch)
    logger.addHandler(fh)

def cactusInit():

    parser = argparse.ArgumentParser(description="Cactus: Discovery pipeline for low-depth NextGen sequence datasets")
    parser.add_argument('--config', dest='cfg', action='store',
                    help='Configuration file containing analysis parameters', default="default.cactus.cfg")
    parser.add_argument('--rebuilddb', dest='rebuilddb', action='store_true',
                    help='Do not run pipeline, rebuild all databases', default=False)

    parser.add_argument('--skip_qc', dest='skip_qc', action='store_true',
                    help='Do not perform QC, will skip output directory removal', default=False)
    parser.add_argument('--skip_assembly', dest='skip_assembly', action='store_true',
                    help='Do not perform Assembly, will skip output directory removal', default=False)
    parser.add_argument('--skip_discovery', dest='skip_discovery', action='store_true',
                help='Do not perform discovery, will skip output directory removal', default=False)

    args = parser.parse_args()

    printCactusIntro()

    #config file valid?
    if not os.path.isfile(args.cfg):
        logger.error('Configuration file was not set')
        raise ValueError('Could not find configuration file')

    #read config file and load our lobbyboy
    config = configparser.RawConfigParser(allow_no_value=True)
    config.read(args.cfg)

    initLogger(config)

    logger = logging.getLogger('cactus')
    logger.info('CACTUS workflow started')

    if args.rebuilddb:
        ncbiObj = ncbi.ncbiObject(args)
    else:
        lb = lobbyBoy(config, args, logger)

        if args.skip_qc is False:
            lb.runQC()

        if args.skip_assembly is False:
            lb.assembly()

        else:
            if args.skip_discovery is False:
                lb.discovery(False)
            else:
                #future reporting implementation here. 
                print('No work')

        logger.info('CACTUS workflow has finished!  ヽ(´▽`)/ ')
        lb.carbonCopyLogFile()

if __name__ == "__main__":
    cactusInit()