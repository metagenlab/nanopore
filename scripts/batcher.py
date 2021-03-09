#!/usr/bin/env python

import sys
import subprocess
from threading import Thread
import shutil
import os
import logging
from argparse import ArgumentParser

class batcher:
    """
    Batcher class to watch a directory for fast5 files and regroup them in batches for basecalling
    parts copied from https://gitlab.com/ModernisingMedicalMicrobiology/CRuMPIT
    """
    def __init__(self,path,outfol,batchsize=1000,batches=0,fileType='fast5'):
        self.job_queue=[]
        self.fileType=fileType
        self.batchsize=int(batchsize)
        self.batches=batches
        self.path=path + '/'
        self.outpath=outfol
        self.stop=False
        self.previousf5s=set()
        self.runNums={}
        self.reruns=0
        logging.basicConfig(level=logging.DEBUG, filename='batcher.log')


    def getPreviousBatches(self):
        """
        This function will check for .batch text files in the output path
        Check if there already is a batch file and add _previous to its name (if the script is rerun) 
        """
        out = subprocess.check_output(["find", f"{self.outpath}","-name","*.batch*"],universal_newlines=True)
        for bf in out.splitlines()[:-1]:
            self.batches+=1
            with open(bf,'rt') as handle:
                for f in handle.read().split('\n'):
                    self.previousf5s.add(f)
            shutil.move(bf, f'{bf}_previous')

    def getFiles(self,path):
        """
        This function will list all fast5s or fastqs in the directory
        """
        maxdepth='10'#Descend at most 10 levels of directories below the command line arguments. 
        try:
            out=subprocess.check_output(["find", f"{path}",'-maxdepth',maxdepth,"-name",f"*.{self.fileType}"],universal_newlines=True)
            out=out.splitlines()
        except subprocess.CalledProcessError as e:
            out=[]
        return out

    def run_batch(self,b):
        if self.fileType=='fast5':
            outfile=f"{self.outpath}/{self.batches}.batch"
        else:
            outfile=f"{self.outpath}/{b[0].split(',')[1]}.batch"
        with open(outfile,'wt') as outf:
            for i in b: outf.write(f"{i}\n") 

    def run(self):
        self.getPreviousBatches()
        self.f5s=self.getFiles(self.path)
        self.npf5s=set(self.f5s)-self.previousf5s
        logging.info(f"run Files {self.npf5s}")
        for f in self.npf5s:
            if len(self.job_queue) >= self.batchsize:
                self.run_batch(self.job_queue)
                self.batches+=1
                self.job_queue=[f]
            else:
                self.job_queue.append(f)
                self.f5s.extend(f)
        self.f5s.extend(self.previousf5s)
        del(self.previousf5s)
        if len(self.job_queue)==0:return
        else:
            self.batches+=1
            self.run_batch(self.job_queue)

    def keepWatch(self):
        while self.stop is False:
            newF5s=self.getFiles(self.path)
            newF5s=set(newF5s)-set(self.f5s)
            for f in newF5s:
                if len(self.job_queue) >= self.batchsize:
                    self.run_batch(self.job_queue)
                    self.batches+=1
                    self.job_queue=[]
                else:
                    self.job_queue.append(f)
                    self.f5s.append(f)
    
    def runWatch(self):
        rcv=Thread(target=self.keepWatch)
        rcv.start()
    
    def stopWatch(self):
        self.stop=True

def batchRun(opts):
    bt=batcher(opts.path, opts.out, opts.batchsize, opts.batches, opts.filetype)
    bt.run()
    if opts.watch==True:
        try:
            bt.runWatch()
        except KeyboardInterrupt:
            bt.stopWatch()

def batchGetArgs(parser):
    parser.add_argument('-p', '--path', required=True,
                                 help='fastq or fast5 input directory path')
    parser.add_argument('-o', '--out', required=True,
                                 help='output path for batch files')
    parser.add_argument('-f', '--filetype', required=False,
                                 help='file type (fast5 or fastq) to batch', default='fast5', type=str)
    parser.add_argument('-bs', '--batchsize', required=False,default=1000,type=int,
                                 help='batch size for basecalling (number of fast5 or fastq files)')
    parser.add_argument('-bn', '--batches', required=False,default=0,type=int,
                                 help='number of previous batch files if restarting')
    parser.add_argument('-w', '--watch', required=False,action='store_true',
                                 help='Option to watch output')
    return parser

if __name__ == "__main__":
    # args
    parser = ArgumentParser(description='CRuMPIT batcher: batch fast5 files or fastq file for nanopore workflow')
    parser=batchGetArgs(parser)
    opts, unknown_args = parser.parse_known_args()
    # run
    batchRun(opts)

