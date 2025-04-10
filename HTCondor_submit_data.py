#!/usr/bin/env python

import os
import optparse
import subprocess
import sys

''' Inputs for the skim code '''
usage = 'usage: %prog [options]'
parser = optparse.OptionParser(usage)
parser.add_option('-i', '--infiles', dest='infiles', help='input list of files (.txt)', default='', type='string')
parser.add_option('-o', '--outfiles', dest='outputfiles', help='output name', default='', type='string')
parser.add_option('-f', '--jobflav', dest='jobflavour', help='job flavour (espresso=20min, microcentury=1h, longlunch=2h, workday=8h, tomorrow=1d, testmatch=3d, nextweek=1w)', default='espresso', type='string')
parser.add_option('-c', '--ncpu', dest='numberofcpu', help='number of cpu requested (2GB per cpu)', default='1', type='int')
parser.add_option('-n', '--njobs', dest='numberofjobs', help='number of jobs to be submitted (integer)', default='1', type='int')
parser.add_option('-s', '--subfiles', dest='subfiles', help='HTCondor submission file', default='HTcondor_sub_data_', type='string')
parser.add_option('-u', '--uncertanties', dest='uncertanties', help='Systematic uncertainties number', default='0', type='int')
parser.add_option('--split', dest='split', help='Split files into multiple parts', action='store_true', default=False)

(opt, args) = parser.parse_args()

inFiles = opt.infiles
outFiles = opt.outputfiles
jobFlavour = opt.jobflavour
nCpu = opt.numberofcpu
nJobs = opt.numberofjobs
subFiles = opt.subfiles
uncerSys = opt.uncertanties
splitFiles = opt.split

''' Read list of files '''
if not os.path.exists(inFiles + '.txt'):
    sys.exit(f"Error: Input file list {inFiles}.txt does not exist.")

listOfFiles = open(inFiles + '.txt', 'r')
Lines = listOfFiles.readlines()
print(f"Number of files: {len(Lines)}")
print(f"Number of jobs: {nJobs}")

# Calculate ratio to split files across jobs
ratio = len(Lines) / nJobs if nJobs != 0 else len(Lines)
if ratio < 1:
    sys.exit("Number of jobs greater than number of files, please reduce the number of jobs.")
ratioint = int(ratio)
print(f"Files per job: {ratio} --> closest integer: {ratioint}")

''' Start the write submission file '''
fsubfile = open(subFiles + ".sub", "w")
command_lines = f'''universe   = vanilla
getenv     = True
executable = htsub.sh
+JobFlavour = "{jobFlavour}"
requirements = ((OpSysAndVer =?= "AlmaLinux9") && (CERNEnvironment =?= "qa"))
RequestCpus = {nCpu}
'''

# If splitting is not enabled, submit one job with all files
if not splitFiles:
    temp = f'''
log        = cond/{subFiles}.log
output     = cond/{subFiles}.out
error      = cond/{subFiles}.err
arguments = {inFiles}.txt {outFiles} 0 0 0 10 5 2.0 0 0 0 {uncerSys}
queue
'''
    command_lines += temp
else:
    ''' Loop over files to create separate input files for each job '''
    for i in range(nJobs):
        outtempfiles = open(f"{inFiles}_part{i}.txt", "w")
        starti = i * ratioint
        endi = (i + 1) * ratioint
        if i == nJobs - 1:
            endi = len(Lines)
        for line in Lines[starti:endi]:
            outtempfiles.write(line)
        outtempfiles.close()

        temp = f'''
log        = cond/{subFiles}_part_{i}.log
output     = cond/{subFiles}_part_{i}.out
error      = cond/{subFiles}_part_{i}.err
arguments = {inFiles}_part{i}.txt {outFiles}_job_{i} 0 0 0 10 5 2.0 0 0 0 {uncerSys}
queue
'''
        command_lines += temp

# Write the complete submission file
fsubfile.write(command_lines)
fsubfile.close()

# Submit the jobs
print(f"Submitting HTCondor jobs using: {subFiles}.sub")
subprocess.call(["condor_submit", subFiles + ".sub"])
