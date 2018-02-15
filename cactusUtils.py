#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import shutil
import os
from itertools import groupby
from collections import defaultdict
import zlib

def compressionRatio(seq):
    bseq = bytes(seq, 'UTF-8')
    complexity = float(len(zlib.compress(bseq)))/len(bseq)
    return complexity

def gapClipper(seq):
    maxPosteriorGapLength = 10
    maxPosteriorGapIslandLength = 6

    nucl = [ 'A', 'C', 'T', 'G', 'a', 'c', 't', 'g']
    replacer = '-'
    pos = 0
    onlyGaps = True
    posteriorGapStart = False
    posteriorGapLength = 0
    posteriorGapGapBoundary = False
    posteriorGapIslandLength = 0

    for letter in seq:
        if onlyGaps == True:
            if letter in nucl:
                onlyGaps = False
                posteriorGapStart = pos
        else:
            if letter == '-':
                posteriorGapLength += 1
            else:
                posteriorGapGapBoundary = pos + 1
                posteriorGapIslandLength += 1

            if posteriorGapIslandLength > maxPosteriorGapIslandLength:
                return seq

            if posteriorGapLength > maxPosteriorGapLength:
                islandLength = posteriorGapGapBoundary - posteriorGapStart
                if islandLength == 0:
                    islandLength = 1

                if posteriorGapStart == 0:
                    mutated =  replacer * islandLength + seq[islandLength:]
                else:
                    mutated =  seq[:posteriorGapStart] + replacer * islandLength + seq[posteriorGapGapBoundary:]

                return mutated
        pos += 1

    return seq

def revcomp(dna):
    bases = 'ATGCN-TACGN-'
    complement_dict = {bases[i]:bases[i+6] for i in range(6)}
    dna = reversed(dna)
    result_as_list = None
    result_as_list = [complement_dict[base] for base in dna]
    return ''.join(result_as_list)

def fastaIter(fastaName):
    fh = open(fastaName)
    faiter = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))
    for header in faiter:
        header = header.__next__()[1:].strip()
        seq = "".join(s.strip() for s in faiter.__next__())
        yield header, seq

def qualIter(qualName):
    fh = open(qualName)
    faiter = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))
    for header in faiter:
        header = header.__next__()[1:].strip()
        qual = " ".join(s.strip() for s in faiter.__next__())
        yield header, qual.split(' ')

def fqIter(fqName):
    f = open(fqName)
    lines = f.readlines()
    f.close()
    records = len(lines)
    lastLine = False
    index = 0

    while lastLine == False:
        try:
            name = lines[index]
            seq = lines[index + 1]
            qual = lines[index + 3]
            index = index + 4
        except IndexError:
            print('index error for index' + str(index))
            sys.exit()
        if index >= records:
            lastLine = True
        else:
            #check if empty
            content = lines[index].strip()
            if not content:
                lastLine = True

        yield name[1:].strip(),seq,qual

def parseAce(aceName):

    contigComposition = defaultdict(dict)
    currentContig = False

    with open(aceName) as aceFile:
        for i, line in enumerate(aceFile):
            lineStart = line[:2]

            if lineStart == 'CO':
                lineSplit = line.split(' ')
                currentContig = lineSplit[1]
                contigComposition[currentContig] = []

            if lineStart == 'RD':
                lineSplit = line.split(' ')
                contigComposition[currentContig].append(lineSplit[1])

    return contigComposition

def replaceLeadingAndTrailingGaps(seq, replaceChar = '~'):
    leadingRemoved = seq.lstrip('-actg')
    trailingRemoved = seq.rstrip('-actg')
    leadingN = len(seq) - len(leadingRemoved)
    trailingN = len(seq) - len(trailingRemoved)
    mutatedSeq = seq[:leadingN].replace('-', replaceChar) + seq[leadingN:]
    trailCut = len(mutatedSeq) - trailingN
    mutatedSeq = mutatedSeq[:trailCut] + replaceChar * trailingN
    return mutatedSeq.upper()

def cactusConsensus(seqCol, idx):
    iupac = dict()

    iupac['A'] = 'A'
    iupac['C'] = 'C'
    iupac['T'] = 'T'
    iupac['G'] = 'G'

    iupac['R'] = 'AG'
    iupac['Y'] = 'CT'
    iupac['S'] = 'CG'
    iupac['W'] = 'AT'
    iupac['K'] = 'GT'
    iupac['M'] = 'AC'
    iupac['B'] = 'CGT'
    iupac['D'] = 'AGT'
    iupac['H'] = 'ACT'
    iupac['V'] = 'ACG'
    iupac['N'] = 'ACGT'
    iupac['AG'] = 'R'
    iupac['CT'] = 'Y'
    iupac['CG'] = 'S'
    iupac['AT'] = 'W'
    iupac['GT'] = 'K'
    iupac['AC'] = 'M'
    iupac['CGT'] = 'B'
    iupac['AGT'] = 'D'
    iupac['ACT'] = 'H'
    iupac['ACG'] = 'V'
    iupac['ACGT'] = 'N'

    n = len(seqCol[0])
    allowed = 'TGCAN-'
    profile = { 'T':[0]*n,'G':[0]*n ,'C':[0]*n,'A':[0]*n, 'N':[0]*n, '-':[0]*n }
    consensus = ''

    for seq in seqCol:
        for i, char in enumerate(seq):
            if char in allowed:
                profile[char][i] += 1

    for i in range(n):
        d = {N:profile[N][i] for N in ['T','G','C','A', '-', 'N']}
        m = max(d.values())
        l = [N for N in ['T','G','C','A', '-', 'N'] if d[N] == m]

        if(m == 0):
            continue
        elif len(l) > 1:
            #remove - from here to cooerce nucleotide assignment
            if '-' in l:
                l.remove('-')
            c = ''.join(sorted(l))
            if c in iupac.keys():
                consensus += iupac[c]
            else:
                consensus += 'N'
        else:
            consensus += l[0]

    degappedConsensus = consensus.replace('-', '')
    return degappedConsensus

def readfq(fp):
    last = None
    while True:
        if not last:
            for l in fp:
                if l[0] in '>@':
                    last = l[:-1]
                    break
        if not last: break
        name, seqs, last = last[1:].partition(" ")[0], [], None
        for l in fp:
            if l[0] in '@+>':
                last = l[:-1]
                break
            seqs.append(l[:-1])
        if not last or last[0] != '+':
            yield name, ''.join(seqs), None
            if not last: break
        else:
            seq, leng, seqs = ''.join(seqs), 0, []
            for l in fp:
                seqs.append(l[:-1])
                leng += len(l) - 1
                if leng >= len(seq):
                    last = None
                    yield name, seq, ''.join(seqs);
                    break
            if last:
                yield name, seq, None
                break
""" Checks if dir exists, if not creates it (will return false), if it does, return True """
def checkDirOrCreate(dirPath):
    dirExists = os.path.isdir(dirPath)
    if dirExists is True:
        return True;
    if dirExists is False:
        os.makedirs(dirPath)
        return False;

def removeDirectoryTree(dirPath):
    shutil.rmtree(dirPath)
