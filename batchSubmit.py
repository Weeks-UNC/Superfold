#!/usr/bin/env python
##################################################################################
# GPL statement:
# This file is part of SuperFold.
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
# Copywrite 2014
# Greggory M Rice
# all rights reserved
# 1.0 build
##################################################################################



# Superfold pipeline
# performs windowed folding and partiton analysis
# copyright 2014, Gregg Rice
# 15 june 2014
# all rights reserved
# pre-alpha build
import subprocess as sub
import shlex, sys, argparse, time, os
from multiprocessing import Pool, Process

def readCommandList(fIN):
    
    out = []
    for line in open(fIN).readlines():
        i = line.rstrip()
        out.append(shlex.split(i))
        
    return out

def submitJob(command):
    pid = os.getpid()
    logname = '.job-{0}.log'.format(pid)
    errname = '.job-{0}.err'.format(pid)
    #print 'running:>% {0}'.format(' '.join(command))
    try:
        log = open(logname,'w')
        err = open(errname,'w')
        x = sub.Popen(command,stdout=log,stderr=err)
        x.wait()
        # remove the log files if the job finishes
        os.remove(logname)
        os.remove(errname)
        return 0
    except:
        print 'Job Failed:>% {0}'.format(' '.join(command))
        return 1
    
def parseArgs():
    arg = argparse.ArgumentParser()
    arg.add_argument('commandList', type=str, help='list of commands in .sh format, each line will be a job')
    arg.add_argument('-n', type=int, default=2, help='number of processors to use, defaults to 2')
    arg.add_argument('--time', action='store_true', help='times how long the job takes to complete')
    o = arg.parse_args()
    return o

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

def updateJobs(currentJobs, jobTrack):
    """
    updates the job tracker. Current jobs holds an array of job objects.
    job track is an array of 01 for job status
    """
    for i in range(len(currJobs)):
        try:
            jobTrack[i] = currJobs[i].is_alive()
        except:
            jobTrack[i] = 0
    return jobTrack


def batchSubmit(jobCommands, nproc=4):
    """
    wrapper for the batch submitter. input is a list of single jobs that have
    been sanitized using the shlex.split utility. The function will wait until
    all jobs have completed before the next command will be processed.
    
    This utility uses the process command. Error files can be tracked for stdout writing.
    These are automatically removed after the command completes but will remain if there
    is a failure.
    
    Command line jobs are assumed to be single threaded
    """

    
    jobs = len(jobCommands)
    
    # set up the job status tracking
    jobTrack = [0] * 24
    # and the job array
    global currJobs # figure out why I need to declare this as a global var. not really important now but may be in the future
    currJobs = [None] * 24
    
    # limit the number of threads that can run concurrently
    maxJobs = nproc
    
    
    while jobs:
        
        
        #update job status
        jobTrack = updateJobs(currJobs, jobTrack)
        runningJobs = sum(jobTrack)
        jobs = len(jobCommands)
        
        # if there are no more jobs wait for running ones to finish
        if not jobs:
            while runningJobs:
                time.sleep(1)
                jobTrack = updateJobs(currJobs, jobTrack)
                runningJobs = sum(jobTrack)
            break
        
        
        #submit jobs if procs available
        if runningJobs < maxJobs:
            a = jobCommands.pop()
            
            for i in range(len(jobTrack)):
                if jobTrack[i] == 0:
                    currJobs[i] = Process(target=submitJob,args=(a,))
                    print 'running:>% {0}'.format(' '.join(a))
                    currJobs[i].start()
                    #continue
                    break
        
        # otherwise sleep this process
        else:
            time.sleep(1)

if __name__ == '__main__':
    
    start = time.time()
    
    # get the command line argument
    args = parseArgs()
    
    # read and parse the commands
    x = readCommandList(args.commandList)
    
    print x
    batchSubmit(x, args.n)
    sys.exit()
    

    
    jobs = len(x)
    totalJobs = int(jobs)
    
    # set up the job status tracking
    jobTrack = [0] * 24
    # and the job array
    currJobs = [None] * 24
    
    # limit the number of threads that can run concurrently
    maxJobs = args.n
    
    
    while jobs:
        
        
        #update job status
        jobTrack = updateJobs(currJobs, jobTrack)
        runningJobs = sum(jobTrack)
        jobs = len(x)
        
        
        # if there are no more jobs wait for running ones to finish
        if not jobs:
            while runningJobs:
                time.sleep(1)
                jobTrack = updateJobs(currJobs, jobTrack)
                runningJobs = sum(jobTrack)
            break
        
        
        #submit jobs if procs available
        if runningJobs < maxJobs:
            a = x.pop()
            
            for i in range(len(jobTrack)):
                if jobTrack[i] == 0:
                    currJobs[i] = Process(target=submitJob,args=(a,))
                    currJobs[i].start()
                    #continue
                    break
        
        # otherwise sleep this process
        else:
            time.sleep(1)
    
    runningJobs = 1
    
    # wait for the rest of the jobs to finish
    
    #while runningJobs:
    #    time.sleep(1)
    #    jobTrack = updateJobs(currJobs, jobTrack)
    #    runningJobs = sum(jobTrack)
        
    
 
    if args.time:
        print 'Total Runtime: {0:.4f} seconds'.format(time.time() - start)
    
