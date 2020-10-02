#!/usr/bin/env python
#  Primary script to run the SuperFold analysis pipeline (in other words RUN THIS SCRIPT)
#
#  - Requires a .map file as an argument (see README for details)
#  - Requires the following programs be executable from any location (i.e. in the PATH):
#        python, Fold, partition, ProbabilityPlot
#  - Take a look at the README for other required modules, installation, and execution help.
#  - Public release 1.0
#  - Copyright Greggory M Rice 2014

#  - Update: Sep 30, 2020: Addressed issues with SHAPE reactivity not displaying on Shannon entropy plot.
#			Altered the Shannon Entropy plot to display Shannon and SHAPE on separate y axes.
#			Fixed np.nan issue with SHAPE plotting.
#			-999 values no longer counted as low SHAPE when determining low SHAPE/low Shannon regions
#			Fixed issue with beginning of sequence being incorrectly included in low SHAPE/low Shannon region
##################################################################################
# GPL statement:
# This file is part of Shapemapper.
#
# SuperFold is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# SuperFold is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with SuperFold.  If not, see <http://www.gnu.org/licenses/>.

# 17 Nov 2014
# all rights reserved
# 1.0 build
##################################################################################
import batchSubmit as batch
import argparse, sys, shlex, os, subprocess, hashlib, time
from RNAtools import dotPlot, CT, padCT, writeSHAPE

# set the plotting environment to be non-interactive
import matplotlib as mtl
mtl.use("Agg")
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

# import arc drawing functions
from drawArcRibbons_simple import writeManyPlot as arcplot
import numpy as np
from PyCircleCompareSF import makeCircle

class shapeMAP:
    def __init__(self,fIN):
        self.seq    = []
        self.ntNum  = []
        self.shape  = []
        self.stdErr = []
        self.zfactor = []
        
        if fIN:
            parsed = self.readFile(fIN)
            
            self.ntNum = map(int,parsed[0])
            self.shape = map(float,parsed[1])
            self.stdErr= map(float,parsed[2])
            self.seq = parsed[3]
            
            if len(parsed.keys()) == 5:
                try:
                    self.zfactor = map(float,parsed[4])
                except:
                    print "formatting error: {0}\nNon float in zfactor column".format(fIN)
                    sys.exit()
            # replace T's with U's
            for i in range(len(self.seq)):
                if self.seq[i] == "T": self.seq[i] = "U"
            
    
    def readFile(self,fIN):
        data = {}
        lineNum = 0
        for line in open(fIN, "rU").readlines():
            x = line.rstrip().split()
            
            # initialize an array obj for the first line
            if lineNum == 0:
                for i in range(len(x)):
                    data[i] = [x[i]]
            
            # populate the array obj
            else:
                for i in range(len(x)):
                    data[i].append(x[i])
            
            lineNum += 1
        
        
        return data


def main():
    
    debug = False
    
    # start the stopwatch
    startTime = time.time()
    
    # check to see if cmd line programs are available
    runCheck()
    
    # parse the input arguments
    args = parseArgs()
    
    
    # set the results directory
    resultsDir = "results_"+args.safeName
    
    print resultsDir
    try:
        os.mkdir(resultsDir)
    except:
        pass
    
    try:
        os.mkdir(resultsDir+"/regions/")
    except:
        pass
    
    # set location of the logfile
    logFile = open("{0}/log_{0}.txt".format(resultsDir),"a")
    sys.stdout = logFile
    print >> sys.stderr, "log file location: {0}/log_{0}.txt".format(resultsDir)
    
    print "\n"*3,"#"*51
    print """#    _____                     ______    _     _  #
#  /  ___|                     |  ___|  | |   | | #
#  \ `--. _   _ _ __   ___ _ __| |_ ___ | | __| | #
#   `--. \ | | | '_ \ / _ \ '__|  _/ _ \| |/ _` | #
#  /\__/ / |_| | |_) |  __/ |  | || (_) | | (_| | #
#  \____/ \__,_| .__/ \___|_|  \_| \___/|_|\__,_| #
#              | |                                #
#              |_|                                #"""                             
    print "#{0: ^49}#".format( "" )
    print "#{0: ^49}#".format( "Superfold ver. Alpha_22-Sept-2014" )
    print "#{0: ^49}#".format( "" )
    print "#{0: ^49}#".format("starting job: " + args.safeName)
    print "#{0: ^49}#".format( time.strftime("%c") )
    print "#{0: ^49}#".format( "" )
    
    print "#"*51,"\n\n# Job Submitted with following attributes:"
    print args, "\n"
    
    # first step is to run the partition function
    print >> sys.stderr, "\nstarting Partition function calculation..."
    partitionPairing = dotPlot()
    
    if not debug:
        partitionPairing = generateAndRunPartition(args.mapObj, args.allConstraints,args.partitionWindowSize, args.partitionStepSize, args.safeName, args.SHAPEslope, args.SHAPEintercept, args.np, args.maxPairingDist)
    
    # write the partition function file
    partitionFileName="{0}/merged_{1}.dp".format(resultsDir, args.safeName)
    
    if not debug:
        partitionPairing.writeDP(partitionFileName)
    
    # debug line
    if debug:
        partitionPairing = dotPlot(partitionFileName)
    
    
    # get the 99% most probable pairs to use as constraints in folding
    probablePairs99 = partitionPairing.requireProb(0.004364805402450088).pairList()
    dsConstraint = {0:[],1:[]}
    for i,j in probablePairs99:
        dsConstraint[0].append(i)
        dsConstraint[1].append(j)
    
    
    # calculate shannon entropy
    bpShannonEntropy = partitionPairing.calcShannon()
    
    # write the shannon entropy in .shape format
    shannonEntropyName = "{0}/shannon_{1}.txt".format(resultsDir, args.safeName)
    writeSHAPE(bpShannonEntropy,shannonEntropyName)
    
    
    
    # generate the folded structure model
    print >> sys.stderr, "starting Fold..."
    
    initialStructure = CT()
    if not debug:
        initialStructure = generateAndRunFold(args.mapObj, args.allConstraints,dsConstraint, args.foldWindowSize, args.foldStepSize, args.safeName,args.SHAPEslope,args.SHAPEintercept,args.np, args.maxPairingDist)
    
    # write the folded structure
    initialStructureFileName = "{0}/merged_{1}.ct".format(resultsDir, args.safeName)
    if not debug:
        initialStructure.writeCT(initialStructureFileName)
    
    #debug
    if debug:
        initialStructure.readCT(initialStructureFileName)
    
    # write final files and figures
    print >> sys.stderr, "drawing figures..."
    
    # add in former pk constraints
    pkPair = []
    for i,j in zip(args.dsConstraints[0],args.dsConstraints[1]):
        pkPair.append((i,j))
    nonPKpairs = initialStructure.pairList()
    finalStructure = CT()
    finalStructure.pair2CT(nonPKpairs+pkPair, seq=initialStructure.seq,name='finalStructure wPKs')
    finalStructureName = "{0}/merged_wPK_{1}.ct".format(resultsDir, args.safeName)
    
    # if there are pk pairs in the file, then write a second ct file with pk's included
    if pkPair != []:
        finalStructure.writeCT(finalStructureName)
    
    # calcualte low shannon/shape regions with expansions based on structure
    
    shannonShapeName = "{0}/shannonShape_{1}.pdf".format(resultsDir, args.safeName)
    
    lowSHAPEregions = mainShannonFunc(args.mapObj.origSHAPE, bpShannonEntropy, shannonShapeName, finalStructure)
    
    
    # plot the arcs along with the shannon/shape reactivities
    arcFileName = "{0}/arcPlot_{1}.pdf".format(resultsDir, args.safeName)
    ensembleRNA_splitPlot(partitionPairing, initialStructure, pk=pkPair, outFile=arcFileName)
    
    # export the structures of the low shannon/shape regions
    maxChar = len(str(lowSHAPEregions[-1][1]))
    
    # file for containing all regions
    ps_comb = "{0}/regions_{1}.ps".format(resultsDir,args.safeName)
    ps_write = open(ps_comb, "w")
    
    if args.noPVclient:
        try:
            import pvclient
        except:
            print "PVclient failed to load"
            args.drawPVclient = False
    
    for i,j in lowSHAPEregions:
        #define file names
        print >> sys.stderr, i,j
        ct_name = "{2}/regions/region_{3}_{0:0>{maxChar}}_{1:0>{maxChar}}.ct".format(i,j,resultsDir,args.safeName,maxChar=maxChar)
        ps_name = "{2}/regions/region_{3}_{0:0>{maxChar}}_{1:0>{maxChar}}.ps".format(i,j,resultsDir,args.safeName,maxChar=maxChar)
        pvclient_name = "{2}/regions/region_{3}_{0:0>{maxChar}}_{1:0>{maxChar}}".format(i,j,resultsDir,args.safeName,maxChar=maxChar)
        
        # write new ct file
        x = finalStructure.cutCT(i,j)
        x.name = "region {0}-{1}, ".format(i,j) + args.name
        x.writeCT(ct_name)
        
        # circle plotting functions
        tmpSHAPE = args.mapObj.origSHAPE[i-1:j]
        tmpZeros = np.zeros_like(tmpSHAPE)
        
        # file lines
        lines = makeCircle(x,x,tmpZeros,tmpSHAPE,{'i':[],'j':[],'correl':[]},[],offset=i)
        
        w = open(ps_name,"w")
        w.write(lines)
        ps_write.write(lines)
        w.close()
        
        if args.noPVclient:
            try:
                pvclient.python_client(x, tmpSHAPE, i, pvclient_name)
            except:
                print "Structure drawing failed Region {0}-{1}".format(i,j)
    ps_write.close()
    
    runtime = "{0:.2f}".format(time.time() - startTime)
    print "\n","#"*51
    print "#{0:^49}#".format("job finished: "+args.safeName)
    print "#{0:^49}#".format( time.strftime("%c"))
    print "#{0:^49}#".format("Total Runtime: " + runtime + " sec.")
    print "#"*51


def ensembleRNA_splitPlot(dpObj, ctObj, pk=None, outFile="arcs.pdf"):
    
    def rgb_int2pct(rgbArr):
        out = []
        for rgb in rgbArr:
            out.append((rgb[0]/255.0, rgb[1]/255.0, rgb[2]/255.0))
        return out
        
    x = dpObj
    y = ctObj
    
    # binning is in log10 scale
    #binning = [0.0,0.09691,0.5228,1.0,2.0]
    binning = [1.5228, 1.0, 0.5228, 0.09691, 0.0]
    # 3% 10% 30% 80% 100%
    
    alphaList = [0.7, 0.7, 0.7, 0.3]
    alphaList = [0.3, 0.7, 0.7, 0.7]
    #alphaList = [1.0, 1.0, 1.0, 1.0]
    #colorList  = ["red", "orange", "yellow", "green","blue", "violet"]
    
    # nat methods palett
    #colorList  = [(215, 25, 28), (253, 174, 97), (171, 221, 164), (43, 131, 186)]
    colorList  = [ (150,150,150), (255,204,0),  (72,143,205) ,(81, 184, 72)  ]
    
    # colorbrewer palett
    #colorList  = [ (43, 131, 186), (171, 221, 164), (253, 174, 97), (215, 25, 28)]
    colorList = rgb_int2pct(colorList)
    
    nucArr = []
    colors = []
    alpha  = []
    
    # bin the pairs by cutoff
    for i in range(0, len(binning)-1):
        
        probPairs = x.requireProb(binning[i],binning[i+1]).pairList()
        
        for pair in probPairs:
            
            temp = np.zeros_like(y.ct)
            temp[pair[0]-1] = pair[1]
            
            #tempCT = RNA.CT()
            #tempCT.pair2CT([pair],y.seq)
        
            #nucArr.append(tempCT.stripCT())
            nucArr.append(temp)
        
            #add a color from the choice list
            colors.append(colorList[i])
            alpha.append(alphaList[i])
    
    if pk:
        for pair in pk:
            temp = np.zeros_like(y.ct)
            temp[pair[0]-1] = pair[1]
            nucArr.append(temp)
            colors.append((0,0,0))
            alpha.append(0.8)
    #nucArr.append(y.stripCT())
    #colors.append("gray")
    
    arcplot(outPath=outFile, pairedNucArr=nucArr, arcColors=colors,seq=y.seq, alpha=alpha, maxDistance=None)


def generateAndRunFold(mapObj, constraints, dsConstraints, windowSize, stepSize, prefix, shapeSlope, shapeIntercept, nprocs, maxDist):
    
    
    #make list of commands, skip those files that have already been calculated
    # run those commands
    # run master model, return structure
    rnaLength = len(mapObj.seq)
    #print rnaLength
    
    # make the directory if it doesn't exist
    dirname = "fold_" + prefix
    try:
        os.mkdir(dirname)
    except:
        pass
    
    # generate the files for the run and the jobs to submit
    jobQueue1 = []
    
    # handle the case where the folding window can cover almost the entire RNA
    if rnaLength-windowSize<200:
        cut_i = 1
        cut_j = rnaLength
        fname = "{0}/{1}_{2}_{3}".format(dirname,prefix,cut_i,cut_j)
        genFiles(mapObj,constraints,dsConstraints,cut_i,cut_j,fname)
        
        foldCMD = "Fold {0}.seq {0}.ct -sh {0}.shape -sm {1} -si {2} -md {3} -C {0}.const -m 100 -w 0".format(fname, shapeSlope,shapeIntercept, maxDist)
        
        jobQueue1.append(shlex.split(foldCMD))
    
    else:
            
        # middle folds
        for i in range(1,rnaLength-windowSize,stepSize):
            cut_i = i
            cut_j = i+windowSize-1
            fname = "{0}/{1}_{2}_{3}".format(dirname,prefix,cut_i,cut_j)
            genFiles(mapObj,constraints,dsConstraints,cut_i,cut_j,fname)
            
            foldCMD = "Fold {0}.seq {0}.ct -sh {0}.shape -sm {1} -si {2} -md {3} -C {0}.const -m 100 -w 0".format(fname, shapeSlope,shapeIntercept, maxDist)
            
            jobQueue1.append(shlex.split(foldCMD))
        
        # 5prime folds
        # 3prime folds
        for i in [-100,-50,50,100]:
            cut5prime_j = windowSize+i
            fname = "{0}/{1}_{2}_{3}".format(dirname,prefix,1,cut5prime_j)
            genFiles(mapObj,constraints,dsConstraints,1,cut5prime_j,fname)
            
            foldCMD = "Fold {0}.seq {0}.ct -sh {0}.shape -sm {1} -si {2} -md {3} -C {0}.const -m 100 -w 0".format(fname, shapeSlope,shapeIntercept, maxDist)
            
            jobQueue1.append(shlex.split(foldCMD))
            
            cut3prime_i = rnaLength-windowSize + i
            fname = "{0}/{1}_{2}_{3}".format(dirname,prefix,cut3prime_i,rnaLength)
            genFiles(mapObj,constraints,dsConstraints,cut3prime_i,rnaLength,fname)
            
            foldCMD = "Fold {0}.seq {0}.ct -sh {0}.shape -sm {1} -si {2} -md {3} -C {0}.const -m 100 -w 0".format(fname, shapeSlope,shapeIntercept, maxDist)
            
            jobQueue1.append(shlex.split(foldCMD))
        
        
    # run the generated jobs
    batch.batchSubmit(jobQueue1,nprocs)
    
    
    # run the master modeler
    
    
    # generate a dummy RNA to align to
    targetRNA = CT()
    targetRNA.pair2CT([],"".join(mapObj.seq))
    
    targetFolderRNAs, baseCount = MasterModel_readAndRenumberAll(targetRNA,dirname, prefix)
    
    pairs = MasterModel_findOverlapPairs(targetFolderRNAs, baseCount)

    masterModelStructure = CT()
    masterModelStructure.pair2CT(pairs,targetRNA.seq,'ConsensusModel')
    
    # return the ct structure
    return masterModelStructure
    

def generateAndRunPartition(mapObj, constraints, windowSize, stepSize, prefix, shapeSlope, shapeIntercept, nprocs, maxDist):
    
    #make list of commands, skip those files that have already been calculated
    #run partition function
    #make list of probability plot commands and run
    #run assembleDP program
    #return dotplot file from assemble
    
    rnaLength = len(mapObj.seq)
    #print rnaLength
    
    # make the directory if it doesn't exist
    dirname = "partition_" + prefix
    try:
        os.mkdir(dirname)
    except:
        pass
    
    # generate the files for the run and the jobs to submit
    jobQueue1 = []
    jobQueue2 = []
    
    # store the names of the files to delete pfs after it has been processed to txt
    pfsNames = []
    
    
    # if the length of the RNA is near the window size, just fold the whole thing at once
    if rnaLength-windowSize < 200:
        cut_i = 1
        cut_j = rnaLength
        fname = "{0}/{1}_{2}_{3}".format(dirname,prefix,cut_i,cut_j)
        genFiles(mapObj,constraints,{0:[],1:[]},cut_i,cut_j,fname)
        
        foldCMD = "partition {0}.seq {0}.pfs -sh {0}.shape -sm {1} -si {2} -md {3} -C {0}.const".format(fname, shapeSlope,shapeIntercept, maxDist)
        
        parseFold = "ProbabilityPlot {0}.pfs {0}.dp -t".format(fname)
        jobQueue1.append(shlex.split(foldCMD))
        jobQueue2.append(shlex.split(parseFold))
        
        pfsNames.append(fname)
    
    # else statement from ln 331, fold the RNA in windows using the options
    else:
        # middle folds
        for i in range(1,rnaLength-windowSize,stepSize):
            cut_i = i
            cut_j = i+windowSize-1
            fname = "{0}/{1}_{2}_{3}".format(dirname,prefix,cut_i,cut_j)
            genFiles(mapObj,constraints,{0:[],1:[]},cut_i,cut_j,fname)
            
            foldCMD = "partition {0}.seq {0}.pfs -sh {0}.shape -sm {1} -si {2} -md {3} -C {0}.const".format(fname, shapeSlope,shapeIntercept, maxDist)
            
            parseFold = "ProbabilityPlot {0}.pfs {0}.dp -t".format(fname)
            jobQueue1.append(shlex.split(foldCMD))
            jobQueue2.append(shlex.split(parseFold))
            
            pfsNames.append(fname)
        
        # 5prime folds
        # 3prime folds
        for i in [-100,-50,50,100]:
            # 5' end
            cut5prime_j = windowSize+i
            fname = "{0}/{1}_{2}_{3}".format(dirname,prefix,1,cut5prime_j)
            genFiles(mapObj,constraints,{0:[],1:[]},1,cut5prime_j,fname)
            
            foldCMD = "partition {0}.seq {0}.pfs -sh {0}.shape -sm {1} -si {2} -md {3} -C {0}.const".format(fname, shapeSlope,shapeIntercept, maxDist)
            parseFold = "ProbabilityPlot {0}.pfs {0}.dp -t".format(fname)
            
            jobQueue1.append(shlex.split(foldCMD))
            jobQueue2.append(shlex.split(parseFold))
            
            pfsNames.append(fname)
            
            # 3' end
            cut3prime_i = rnaLength-windowSize + i
            fname = "{0}/{1}_{2}_{3}".format(dirname,prefix,cut3prime_i,rnaLength)
            genFiles(mapObj,constraints,{0:[],1:[]},cut3prime_i,rnaLength,fname)
            
            foldCMD = "partition {0}.seq {0}.pfs -sh {0}.shape -sm {1} -si {2} -md {3} -C {0}.const".format(fname, shapeSlope,shapeIntercept, maxDist)
            parseFold = "ProbabilityPlot {0}.pfs {0}.dp -t".format(fname)
            
            jobQueue1.append(shlex.split(foldCMD))
            jobQueue2.append(shlex.split(parseFold))
            
            pfsNames.append(fname)
        
    
    # run the generated jobs
    batch.batchSubmit(jobQueue1,nprocs)
    batch.batchSubmit(jobQueue2,nprocs)
    
    # clean up the fnames
    for fname in pfsNames:
        try:
            os.remove(fname+".pfs")
        except:
            pass
    
    
    # run the assemble DP routine
    dpObject = mainAssemble(dirname,trim=300)
    
    return dpObject
    
    
def runCheck():
    """
    check to see if necessary commands are available to call
    """
    neededCmds = ["Fold", "partition", "ProbabilityPlot"]
    count = 0
    
    for each in neededCmds:
        try:
            out, err = "", ""
            subprocess.call([each], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            count += 1
        except OSError:
            print "Program {0} not found in the path".format(each)
    
    if len(neededCmds) != count:
        print "...exiting"
        sys.exit()
    
    try:
        # get the data path and test to see if it is real
        datapath = os.environ.get("DATAPATH")
        os.listdir(datapath)
    except:
        print "DATAPATH is not set. RNAstructure will not run"
        print "...exiting"
        sys.exit()


def parseArgs():
    
    def parseTXT(fIN):
        data = {}
        lineNum = 0
        for line in open(fIN, "rU").readlines():
            x = line.rstrip().split()
            
            try:
                x = map(int,x)
            except:
                print "Unexpected character in {0}...exiting".format(fIN)
                return 0
            
            # initialize an array obj for the first line
            if lineNum == 0:
                for i in range(len(x)):
                    data[i] = [x[i]]
            
            # populate the array obj
            else:
                for i in range(len(x)):
                    data[i].append(x[i])
            
            lineNum += 1
        
        
        return data
    
    def shapeEnergy(x,slope,intercept):
        out = []
        for i in x:
            if i<-500:
                out.append(-999)
                continue
            if -500<i<0:i=0.0
            g = slope*np.log(i+1)+intercept
            out.append(g)
        return np.array(out)
    
    def differentialEnergy(x,a,b,zfactor):
        zfactor = np.array(zfactor)
        
        # set -inf and nan to non significant zfactors
        zfactor[np.isinf(zfactor)] = -999
        zfactor[np.isnan(zfactor)] = -999
        
        out = []
        for i in x:
            if i<-500:
                out.append(-999)
                continue
            #slow fit; linear so easy
            if i >= 0:
                g = a*i+b
            #slow case
            if i < 0:
                g = -999#c*i+d
                #c = abs(i)
                #if c > 0.2:
                #    g = -0.620*np.log(func(c,f1[0],f1[1],f1[2])/func(c,f2[0],f2[1],f2[2]))
                #else:
                #    g = 2.57275425027 * c + 0.0239
            out.append(g)
        return np.array(out)*(zfactor>0)
    
    def fakeSHAPE(shape,diff,slope,intercept):
        #shape_f,diff_fmap=(float(shape)),map(float(diff))
        #print shape
        out = []
        for i,j in zip(shape,diff):
            if i < -500 and j<-500:
                out.append(-999)
                continue
            if i < -500:
                g = j
            if j < -500:
                g = i
            if i >-500 and j>-500:
                g = i + j
            #out.append(g)
            out.append(np.exp( (g-intercept) / slope ) -1.0)
        return out

    arg = argparse.ArgumentParser(description="SuperFold takes a windowing approach to break up the folding of large RNAs. Dividing the folding of a large RNA into smaller segments allows modern multi-core workstations to model RNA structures in a modest amount of clock-time. See README file for further details and file descriptions", epilog="SuperFold v1.0 by Gregg Rice ( gmr@unc.edu )")
    arg.add_argument('mapFile', type=str, help='SHAPE-MaP .map file, trimmed to fold desired sequence')
    arg.add_argument('--ssRegion', type=str, help='file containing forced single stranded regions')
    arg.add_argument('--pkRegion', type=str, help='text file containing pseudoknotted basepairs:\npaired nucleotides on each line seperated by whitespace. e.g.\n\n1 50(newline)2 49(newline)...')
    arg.add_argument('--np', type=int, default=2, help='number of processors to use, default:2')
    arg.add_argument('--SHAPEslope',type=float,default=1.8, help='SHAPE pseudofreeenergy slope, default:1.8')
    arg.add_argument('--SHAPEintercept',type=float,default=-0.6,help='SHAPE pseudofreeenergy intercept, default:-0.6')
    arg.add_argument('--differentialFile',type=str, help='.mapd file containing NMIA-1M6 calculated values')
    arg.add_argument('--differentialSlope',type=float,default=2.1,help='differential SHAPE slope, default:2.1')
    arg.add_argument('--trimInterior',type=int,default=300, help='number of nucleotides to trim to improve negative end effects during windowed folding, default:300')
    arg.add_argument('--partitionWindowSize',type=int,default=1200, help='length of the partition window size, default:1200')
    arg.add_argument('--partitionStepSize',type=int,default=100, help='spacing between partition windows, default:100')
    arg.add_argument('--foldWindowSize',type=int,default=3000, help='length of the Fold window size, default:3000')
    arg.add_argument('--foldStepSize',type=int,default=300, help='spacing between Fold windows, default:300')
    arg.add_argument('--maxPairingDist', type=int, default=600, help='Maximum pairing distance for partition and Fold, default:600')
    arg.add_argument('--noPVclient', action='store_false', help="Don't draw secondary structures using PVclient")
    
    o = arg.parse_args()
    
    o.name = o.mapFile
    # read the shapemap file in
    mapfile = shapeMAP(o.mapFile)
    o.mapObj = mapfile
    
    o.mapObj.origSHAPE = np.array(o.mapObj.shape)
    
    if o.differentialFile:
        o.diffMapFile = shapeMAP(o.differentialFile)
        
        calcSHAPEenergy = shapeEnergy(mapfile.shape,o.SHAPEslope,o.SHAPEintercept)
        calcDiffEnergy  = differentialEnergy(o.diffMapFile.shape,o.differentialSlope,0,o.diffMapFile.zfactor)
        
        calcFakeShape   = fakeSHAPE(calcSHAPEenergy,calcDiffEnergy,o.SHAPEslope,o.SHAPEintercept)
        
        o.mapObj.shape = calcFakeShape
        
        #print o.diffMapFile.zfactor
    
    # read the constraint files
    
    # first the ss constraints
    try:
        ss = parseTXT(o.ssRegion)
        
        # catch the exception
        if ss == 0:
            sys.exit()
    except:
        if o.ssRegion == None:
            ss = {0:[]}
        else:
            sys.exit()
        
    o.ssConstraints = ss
    
    # then the ds constraints
    
    try:
        ds = parseTXT(o.pkRegion)
        
        # catch the exception
        if ds == 0:
            sys.exit()
        
        # check to make sure that the pks have partners
        if len(ds[0]) != len(ds[1]):
            print "pkRegion file incorrectly formatted...exiting"
            sys.exit()
    except:
        # if file not given, fill empty bps
        if o.pkRegion == None:
            ds = {0:[],1:[]}
        else:
            sys.exit()
    
    o.dsConstraints = ds
    
    # concatonate all the constraints
    allConst = ss[0] + ds[0] + ds[1]
    o.allConstraints = allConst
    
    
    # generate a safe short name for the constraint, shape, and sequence files
    # by hashing the names of the input values, changing any input will make a new hash
    m = hashlib.md5()
    m.update(str(o.mapFile))
    m.update(str(o.ssRegion))
    m.update(str(o.pkRegion))
    m.update(str(o.differentialFile))
    m.update(str(o.foldWindowSize))
    m.update(str(o.partitionWindowSize))
    m.update(str(o.maxPairingDist))
    m.update(str(o.partitionStepSize))
    
    o.safeName = o.mapFile[:] +"_"+ m.hexdigest()[:4]
    
    
    return o


def genFiles(mapObj, ssConstraints, dsConstraints, ntStart, ntEnd, fName):
    
    def shapeFile(SHAPEdata, fOUT):
        """
        writes a .SHAPE file compatable with RNAstructure
        """
        w = open(fOUT, "w")
        
        for i in range(len(SHAPEdata)):
            w.write("{0}\t{1}\n".format(i+1,SHAPEdata[i]))
        w.close()
        return fOUT
        
    def seqFile(sequence, fOUT, name=None):
        """
        write a seq file compatable with RNAstructure format
        """
        
        if not name:
            name = str(fOUT)
        w = open(fOUT, "w")
        w.write(";\n\n{0}\n\n".format(name))
        
        for i in range(len(sequence)):
            w.write(sequence[i])
            if (i+1) % 50 == 0:
                w.write("\n")
            elif (i+1) % 10 == 0: w.write(" ")
        w.write("1\n")
        w.close()
        return fOUT
    
    def constraintFile(ssConstraint, dsConstraint, fOUT):
        """
        writes a constraint file compatable with RNAstructure
        """
        
        w = open(fOUT, "w")
        w.write("DS:\n-1\nSS:\n")
        for const in ssConstraint: w.write("{0}\n".format(const))
        w.write("-1\nMod:\n-1\nPairs:\n")
        
        for i in range(len(dsConstraint[0])):
            line = "{0} {1}\n".format( dsConstraint[0][i], dsConstraint[1][i])
            w.write(line)
                
        w.write("-1 -1\nFMN:\n-1\nForbids:\n-1 -1\nMicroarray Constraints:\n0\n")
        w.close()
        
        
        return fOUT
    
    ### renumber constraints
    
    renumSS = []
    renumDS = {0:[],1:[]}
    
    for i in ssConstraints:
        if i >= ntStart and i <=ntEnd:
            renumSS.append(i-(ntStart-1))
    
    for i,j in zip(dsConstraints[0],dsConstraints[1]):
        if i >= ntStart and j >= ntStart:
            if i <= ntEnd and j<=ntEnd:
                renumDS[0].append(i-(ntStart-1))
                renumDS[1].append(j-(ntStart-1))
                
        # if only one of the pairs is within the desired window hold it out
        elif j >= ntStart and j <= ntEnd:
            renumSS.append(j-(ntStart-1))
        elif i >= ntStart and i <= ntEnd:
            renumSS.append(i-(ntStart-1))
        
    
    renumSS.sort()
    
    
    # SHAPE file
    shapeFile(mapObj.shape[ntStart-1:ntEnd],fName+".shape")
    # sequence file
    seqFile(mapObj.seq[ntStart-1:ntEnd], fName+".seq",name=fName)
    # constraint file
    constraintFile(renumSS,renumDS,fName+".const")
    

def mainAssemble(folderPath, trim=300):
    """
    input is the path of a folder of dp files. will return an assembled
    dotPlot object with all of the windows. Trim defines how much of each
    window to cut out of each of the dp files.
    
    Must have the nt start and ending numbers as the last two parts of the
    file name seperated by underscores
    """
    
    targetDP = {}
    
    # load all the dp files into memory
    numFiles = len(filter( lambda x:x[-2:]=="dp",os.listdir(folderPath)  ) )
    
    num = 1
    print "reading files..."
    
    for dpFileName in os.listdir(folderPath):
        #ignore all files that are not dp files
        if dpFileName[-2:] != "dp": continue
        
        #print "reading {0}...".format(dpFileName)
        dp = dotPlot(folderPath + "/" + dpFileName)
        
        start = int(dpFileName.split("_")[-2])
        end   = dpFileName[:-3].split("_")[-1]
        progress(num, numFiles)
        num +=1
        
        targetDP[(int(start), int(end))] = dp
        
    
    # find the first and last dotplot files in sequence
    # we only want to trim one side of those
    firstDP = min([i[0] for i in targetDP.keys()])
    lastDP  = max([i[1] for i in targetDP.keys()])
    
    
    #container for the final dp structure
    finalDP = dotPlot()
    finalDP.name = "assembled dotPlot"
    finalDP.length = lastDP
    
    coverage = []
    num = 1
    
    print "trim and resorting..."
    # average slipped pairs, renumber, and trim the ends
    # the dpKey is the starting number of window
    if len(targetDP.keys()) == 1:
        return targetDP[targetDP.keys()[0]]
    
    for dpKey in targetDP.keys():
        
        #print "resorting window at {0}...".format(dpKey)
        progress(num,numFiles)
        num+=1
        # remove some of the very lowly probable base pairs to simplify later calculations
        # targetDP[dpKey] = targetDP[dpKey].requireProb(4)
        
        # preaverage the slipped base pairs
        #targetDP[dpKey].averageSlippedBPs(predictedOnly=False)
        
        # trim 3' end of the first structure in sequence
        if dpKey[0] == firstDP:
            #print dpKey, "3prime"
            targetDP[dpKey] = targetDP[dpKey].trimEnds(trim, which='3prime')
            coverage.append((dpKey[0],dpKey[1]-(trim-1)))
        # trim 5' end of the first structure in sequence
        elif dpKey[1] == lastDP:
            #print dpKey, "5prime"
            targetDP[dpKey] = targetDP[dpKey].trimEnds(trim, which='5prime')
            coverage.append((dpKey[0]+(trim-1),dpKey[1]))
        # trim middle dplots at both ends
        else:
            #print dpKey, "both", trim
            targetDP[dpKey] = targetDP[dpKey].trimEnds(trim)
            coverage.append((dpKey[0]+(trim-1),dpKey[1]-(trim-1)))
        
        # reset the numbers to the correct global numbers
        targetDP[dpKey].dp['i'] += dpKey[0] - 1
        targetDP[dpKey].dp['j'] += dpKey[0] - 1
        
        # append to the final dp file
        finalDP.dp['i'] = np.append(finalDP.dp['i'], targetDP[dpKey].dp['i'])
        finalDP.dp['j'] = np.append(finalDP.dp['j'], targetDP[dpKey].dp['j'])
        finalDP.dp['logBP'] = np.append(finalDP.dp['logBP'], targetDP[dpKey].dp['logBP'])
        
    # save a checkpoint file
    #pickle.dump(finalDP, open("finaldp1.bin", "wb"))
    
    finalDP = concatonateDP(finalDP, coverage)
    
    
    return finalDP


def concatonateDP(dpObj, coverage):
    """
    merges duplicate entries in a dp file by averaging
    """
    def calcCoverage(i, j, coverage):
        n = 0
        if j-i >= 600: return 500
        for a,b in coverage:
            if a <= i <= b and a<= j <= b:
                n+=1
        return n
        
    #print dpObj
    #print coverage
    dpObj = dpObj.requireProb(2)
    
    # make a shortcut
    dp = dpObj.dp
    
    # maxDist is the max search space
    maxDist = 600
    
    # construct the return object
    outObj        = dotPlot()
    outObj.length = dpObj.length
    outObj.name   = dpObj.name
    
    # grab the potential pairs to limit looping space
    pairs = set(dpObj.pairList())
    
    fullSize = len(pairs)
    outObj.dp['i'] = np.zeros(fullSize)
    outObj.dp['j'] = np.zeros(fullSize)
    outObj.dp['logBP'] = np.zeros(fullSize)
    outObj.dp['coverage'] = np.zeros(fullSize)
    
    n = 0
    
    oldFilter = np.zeros_like(dp['logBP'])
    dp['logBP'] = 10**( -dp['logBP'])
    
    print "merging dotplots..."
    for i,j in pairs:
        i = int(i)
        j = int(j)
        
            
        # remove any long distance base pairs
        #if j-i > maxDist: continue
        dpFilter = ( dp['i']== i ) * ( dp['j'] == j )
        
        
        oldFilter += dpFilter
        # skip nonexisting pairs
        entries = np.sum(dpFilter)
        outObj.dp['i'][n] = i
        outObj.dp['j'][n] = j
        outObj.dp['logBP'][n] = 0.0
        outObj.dp['coverage'][n] = calcCoverage(i, j, coverage)
        if entries == 0: continue
        #elif entries == 1:
            #outObj.dp['i'] = np.append(outObj.dp['i'], i )
            #outObj.dp['j'] = np.append(outObj.dp['j'], j )
            #outObj.dp['logBP'] = np.append(outObj.dp['logBP'], dp['logBP'][dpFilter] )
        else:
            #outObj.dp['i'] = np.append(outObj.dp['i'], i )
            #outObj.dp['j'] = np.append(outObj.dp['j'], j )
            #outObj.dp['logBP'] = np.append(outObj.dp['logBP'], np.average(dp['logBP'][dpFilter]) )
            #outObj.dp['logBP'][n] = np.average(dp['logBP'][dpFilter])
            outObj.dp['logBP'][n] = np.sum(dp['logBP'][dpFilter])/outObj.dp['coverage'][n]
            #print np.average(dp['logBP'][dpFilter]), dp['logBP'][dpFilter]
        
        #remove already searched from full list
        #debugging
        n+=1
        progress(n, fullSize)
        if n%1000 == 0:
            dp['i'] = np.delete(dp['i'],np.nonzero(oldFilter))
            dp['j'] = np.delete(dp['j'],np.nonzero(oldFilter))
            dp['logBP'] = np.delete(dp['logBP'],np.nonzero(oldFilter))
            oldFilter = np.zeros_like(dp['logBP'])
            #print n, len(dp['logBP'])
    
    print "DONE!"
    
    
    outObj.dp['logBP'] = -np.log10(outObj.dp['logBP'])
    #print np.sum(outObj.dp['coverage'] ==2)
    
    return outObj
    

def progress(num,outof):
    num = float(num)
    outof = float(outof)
    width = 30
    line = '['
    
    meter =  ''.join(['=']*int(num/outof*width)) + ''.join([' ']*int((outof-num)/outof*width))
    line += meter[:width-4]
    line += ']'
    
    line += ' %s / %s' % (int(num),int(outof))
    if num == 1:
        sys.stdout.write(line)
    elif num==outof:
        sys.stdout.write('\r'+line)
        sys.stdout.flush()
        sys.stdout.write("\n")
    else:
        sys.stdout.write('\r'+line)
        sys.stdout.flush()


def MasterModel_readAndRenumberAll(targetCT,folder, prefix):
    """go through the folder of ct files
       read and renumber them in the context
       of the big sequence"""
    targetCTfiles = []
    
    #initialize bp overlap object
    #stores how many times particular base is possible
    bpoverlap = np.zeros(len(targetCT.ct))
    
    for i in os.listdir(folder):
        # skip non-ct files
        if i[-2:] != 'ct': continue
        
        # check file prefix
        checkLen = len(prefix)
        if i[:checkLen] != prefix: continue
        
        #load the fold
        loadedCT = CT(folder+"/"+i)
        #renumber it in context of the target genome
        paddedFold,pos = padCT(loadedCT,targetCT,giveAlignment=True)
        
        #add one to each position if it is contained in
        #the aligned ct
        for i in range(pos,pos+len(loadedCT.ct)):
            bpoverlap[i]+=1
        
        #add it to the array
        targetCTfiles.append(paddedFold)
    
    return targetCTfiles, bpoverlap

def MasterModel_findOverlapPairs(ctObjectList, baseCount):
    """
    goes over the list of CT objects and
    calculates the overlap of basepairs
    """
    
    # first step is to go through the CT objects and collect
    # all the common basepairs. Format is a a dict {(i,j):n}
    # where n is how many times the combination i,j has appeared
    
    bpairs = {}
    
    for rna in ctObjectList:
        # for the nts in the ct file
        for nuc in range(len(rna.ct)):
            # if it is not paired, skip it
            if rna.ct[nuc] == 0: continue
            
            # set i and j from the defined nuc variable
            i, j = nuc+1, rna.ct[nuc]
            
            # skip those j greater than i, no duplicates
            if j < i: continue
            
            # initialize the dict pair if it doesnt exist...
            if (i,j) not in bpairs:
                bpairs[(i,j)] = 1
            
            # ...add to it if it does
            else:
                bpairs[(i,j)] += 1
    
    # second step is to go through the generated pairs and if they meet
    # the cutoff (appearing in more than half of the possible windows)
    # then add them to the final 
    outpairs = {}

    for key in bpairs.keys():
        #print key
        # key[0] = i, key[1] = j
        if bpairs[key]/float(min(baseCount[key[0]-1],baseCount[key[1]-1])) > 0.5:
            outpairs[key] = bpairs[key]

    return outpairs.keys()

############################
# SHANNON SHAPE SEARCH FUNCTIONS
############################

def mainShannonFunc(shapeReact, shannonEntropy, saveName, ctStruct):
    
	def findTransitions(react,cutoff):
		transitions = []
		prev = 1
		for i in range(0,len(react)):
			curr = (react[i] > cutoff)
			#Update V1.1: Counting -999 (nan) values as high SHAPE
			#This is to prevent -999 tracks being called as low SHAPE
			if np.isnan(react[i]):
				curr = True				
			if curr != prev:
				transitions.append((i,curr))
				prev = curr
		return transitions

	def cullTransitions(transitions,minLength):
		dist = []
		culled = []
		skip = False
		# dist array is n-1 in length to the transitions
		# it is the sequence distance length
		# transitions = [(2,T),(5,F),(10,T)]
		# dist = [3,5]
		for i in range(len(transitions)-1):
			dist.append(transitions[i+1][0]-transitions[i][0])
		for i in range(len(dist)):
			if skip:
				skip=False
				continue
			if dist[i]<minLength:
				skip=True
				continue
			culled.append(transitions[i])
		culled.append(transitions[-1])
		return culled

	def selectCutsites(trans1, trans2, seqLength, minlength=40):
		def genStretch(shortlist,seqLength,high=True,pstat=False):
			# creates an array containing the potential good cut sites as a list of 1,0
			# based on the high low values of the culled transitions
			# high flag inverses the transition logic
			addnums = 1-abs(shortlist[0][1])
			#addnums = abs(shortlist[0][1])
			if not high:
				addnums = abs(addnums - 1)
			addList = []
			#loop through all nts in the sequence
			#V1.1 Update: loop had started at 1, this skipped the first
			#region and resulted in mis-assigned regions
			for i in range(0,seqLength+1):
				if pstat:print addnums,
				if addnums:
					addList.append(i)
				if not addnums:
					addList.append(0)
				#if the current nt matches one of the transition starts
				#set the addnums flag to match the boolean of the transition
				for j in shortlist:
					if j[0] == i:
						addnums=abs(j[1])
						if not high:
							addnums=abs(1-j[1])
			return addList
		t1 = np.array(genStretch(trans1,seqLength,high=False))
		t2 = np.array(genStretch(trans2,seqLength,high=False,pstat=False))
		
		# truth logic to find union
		combined = (t1*t2)>0
		#return combined
	
		# find the transitions again, cull for the minlength, and regenerate the stretch
		transitions = cullTransitions(findTransitions(combined,0.5),minlength)
		culled_combined = genStretch(transitions,seqLength)
	
		return culled_combined
	
	def movingWindow(dataIn, degree):
		"""
		windowsize is 2*degree+1
		"""
		out = np.zeros(len(dataIn))
		dataIn = np.array(dataIn)
		dataIn[dataIn<-500] = np.nan
		# window = 2*degree + 1
		for i in range(degree, len(dataIn)-degree):
			#### v1.1 update, medians now calculated with nanmedian not median
			out[i] = np.nanmedian(dataIn[i-degree:i+degree+1])
	
		# pad the 5' end
		for i in range(0,degree):
			out[i] = float(out[degree])
	
		# pad the 3' end
		for i in range(len(dataIn)-degree, len(dataIn)):
			out[i] = float(out[len(dataIn)-degree-1])
	
		#for i in out: print i
		return out

	shannonUnscaled = shannonEntropy
	shapeUnscaled = shapeReact

	# get a moving window of the SHAPE and Shannon data
	shape = movingWindow(shapeUnscaled, 25)
	shannon = movingWindow(shannonUnscaled, 25)

	shapeTrans = findTransitions(shape,np.nanmedian(shape)/1)
	#shapeTrans = findTransitions(shape,0.25)
	#shannonTrans = findTransitions(shannon,np.median(shannon)*1)
	shannonTrans = findTransitions(shannon,np.median(shannon))
	
	shape_culled = cullTransitions(shapeTrans,10)
	shannon_culled = cullTransitions(shannonTrans,10)

	# find the union of the two culled transitions
	a = np.array(selectCutsites(shape_culled,shannon_culled,len(shape)))
	b=a>0
	
	print "####### shannonShapeTransitions #######"
	#print findTransitions(b,0.5)
	regions = findTransitions(b,0.5)

	# optional: print out the range to stdout
	#for i in range(len(b)):
	#    print i+1,b[i]*1.0
	#print shannon_culled


	#plotting functions

	#define x axis positions
	x = np.linspace(0,len(shape),len(shape))

	#set up plots and xtick locations
	fig = plt.gcf()
	fig.set_size_inches(10.5,8.1)
	fig.set_size_inches(11.9,3.5)
	fig.set_size_inches(24,7)
	pdf = PdfPages(saveName)
	plt.subplots_adjust(left=0.05,right=0.95,bottom=0.06,top=0.98)
	ax = plt.subplot(111)

	az = fig.add_axes()

	#plot labels
	plt.xlabel('Nucleotide',fontsize=10)
	plt.ylabel('Shannon entropy',fontsize=10)
	#plt.setp(ax.get_xticklabels(),fontsize=8,color='w')
	#plt.setp(ax.get_yticklabels(),fontsize=8,color='w')
	plt.setp(ax.get_xticklabels(),fontsize=8)
	plt.setp(ax.get_yticklabels(),fontsize=8)

	#set minor and major tick locations
	majLoc   = MultipleLocator(500)
	minorLoc = MultipleLocator(100)
	yMajLoc  = MultipleLocator(0.1)
	YminorLoc= MultipleLocator(0.05) 
	yMajLoc2  = MultipleLocator(0.1)
	YminorLoc2= MultipleLocator(0.05)
	ax.xaxis.set_major_locator(majLoc)
	ax.xaxis.set_minor_locator(minorLoc)
	ax.yaxis.set_major_locator(yMajLoc)
	ax.yaxis.set_minor_locator(YminorLoc)


	ax.spines['right'].set_color('none')
	ax.spines['top'].set_color('none')
	ax.spines['left'].set_color('none')

	#ax.get_xaxis().set_visible(False)
	#ax.get_yaxis().set_visible(False)

	#plot range
	plt.ylim(0,2)
	plt.xlim(0,5000)

	#plot the highlighted range
	lowShannonLowSHAPE = plt.fill_between(range(len(b)),b*3.0,0,alpha=0.2,color='blue', label = "Low Shannon low SHAPE")

	#    for i in range(len(b)):
	#        print i+1, float(b[i])
	
	#### V1.1 Update #####
	#plt the shannon entropy and it's median
	plt.plot(x,np.median(shannon)*np.ones(len(x)),color='y')
	#plt.plot(shannon,color='r')
	shannonPlot = plt.fill_between(x,shannon,np.zeros_like(shannon),facecolor='brown',interpolate=True, label = "Shannon")

	#twin the x axis
	ax2 = ax.twinx()
	#plot the shape reactivity above the shannon
	shape_med = np.nanmedian(shape)
	#### V1.1 Update
	#overwrite all np.nan values with median
	#this prevent errors in the plotting
	shapeForPlotting = shape[:]
	for i in range(len(shapeForPlotting)):
		if np.isnan(shapeForPlotting[i]):
			shapeForPlotting[i] = shape_med
	#plt.fill_between(x,shape,shape_med, where=shape>shape_med, facecolor='red',interpolate=True)
	#plt.fill_between(x,shape,shape_med, where=shape<shape_med, facecolor='blue',interpolate=True)
	shapePlot = ax2.fill_between(x,shapeForPlotting,shape_med, facecolor='black',interpolate=True, label = "SHAPE")
	
	plt.ylabel('SHAPE reactivity')
	plt.ylim(-1,1)
	plt.legend([lowShannonLowSHAPE,shapePlot,shannonPlot],["Low SHAPE low Shannon","SHAPE","Shannon"],loc=2)
	ax2.yaxis.set_major_locator(yMajLoc2)
	ax2.yaxis.set_minor_locator(YminorLoc2)
	
	#### End of V1.1 Update #####
	
	#save first 5k
	pdf.savefig(dpi=300,transparent=True)

	#change x range and save next 5k
	plt.xlim(5000,len(x))
	pdf.savefig(dpi=300,transparent=True)

	#change to whole range to plot in the pop up window
	plt.xlim(0,len(x))
	pdf.close()
	plt.savefig(saveName,dpi=300,transparent=True)


	#####
	# expand regions to cover the entire secondary structure
	#####
	n = 0
	regionPair = []

	#print >> sys.stderr, regions
	#print >> sys.stderr, len(regions)

	try:
		if regions[0][1] == False:
			#V1.1 Update, prevent an issue with the first nt
			if regions[0][0] != 0:
				regionPair.append([1, regions[0][0]])
			n = 1
	
		for i in range(n,len(regions)-n,2):
			#print >> sys.stderr, i,n
			regionPair.append([regions[i][0], regions[i+1][0]])
	
		if regions[-1][1] == True:
			regionPair.append([regions[-1][0], len(shape)])
	
	except:
		print "No Defined Regions"
		regionPair.append([1,len(shape)])

	print "Unexpanded regions:"
	for i,j in regionPair:
		print i,"-",j
	
	allPair = ctStruct.pairList()

	expandedRegion = []
	for i,j in regionPair:
		new_i = i
		new_j = j
	
		for k,l in allPair:
		
			if k<i and l>i:
				if k< new_i:
					new_i = k
				if l> new_j:
					new_j = l
			if k < j and l > j:
				if k< new_i:
					new_i = k
				if l> new_j:
					new_j = l
			
		expandedRegion.append((new_i,new_j))

	print "Expanded regions:"
	for i,j in expandedRegion:
		print i,"-",j

	return expandedRegion
    


if __name__ == '__main__':
    main()


