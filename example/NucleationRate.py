import re
import os
import math
import numpy as np
from optparse import OptionParser


def readFileAsLines(filename):
    fileLines = open(filename, 'r').readlines()
    open(filename, 'r').close()
    return(fileLines)
    
path = os.getcwd() + '/'
os.chdir(path)

#folderlist = os.listdir(path)
#for filename in folderlist:
#    if re.search('slurm',filename):  # or rename the jobout file
#        jobout = filename

parser = OptionParser()
parser.add_option('--datafile',type = str,default ='in.data',help = "Name of data file(default: %default)")
parser.add_option('--logfile', type = str,default = 'slurm-7648355.out',help = 'Name of slurm job output file (default: %default)')
(options, args) = parser.parse_args()
dataFileInput = options.datafile
jobout = options.logfile

ffsInputLines = readFileAsLines('ffs.input')
lammpsInputLines = readFileAsLines('lammps.input')
joboutLines = readFileAsLines(jobout)
structureLines = readFileAsLines(dataFileInput)  # input strucutre file
trajectoryLines = readFileAsLines('trajectory.out.txt')



for line in ffsInputLines:
    if re.match('equilibrium ', line):
        equilibriumSteps = float(line.split(' ')[1])
        
    if re.match('lambda ', line):
        lambdaList = line.split(' ')[1:len(line.split(' '))]
        num_interface = (len(line.split(' '))-3)


for line in lammpsInputLines:
    if re.match('timestep ', line):  # timestep 5fs
        timestep = float(line.split(' ')[1])
 
initialCrossing = 0
universe,step,targetLambda,volumes = [],[],[],[]
crossing=[]
for i in range(num_interface):
    crossing.append([0,0])

    
for line in joboutLines:    
    if re.match('\[date=\d*\]\s\[universe=\d*\]\s\[steps=\d*\]\s+:\s\d*\\s+\.\.\.\s+\d*\\s', line):            
        universe.append(re.split('\[|\]|=|\s',line)[6])
        step.append(re.split('\[|\]|=|\s',line)[10])
        targetLambda.append(re.split('\[|\]|=|\s',line)[15])
    if re.match('\s+\d*\s+\_+\s\(\_+\)\s+>==\d*\s+\d*==>\s+\d*\s\(xyz.0\_*\d*\_\d*',line):
        initialCrossing = initialCrossing + 1
    if re.match('\s*(?:\d+\s+)?\d+\s*\(xyz\.(\d+)__\d+_\d+\)\s*>==\s*\d+\s*\d+\s*==>\s*\d+\s*\((.+)\)\s*',line):
        split = re.split('\.|__',line)
        if 'xyz' in split[-3]:
             crossing[int(split[1])][0] += 1 # successful crossing
             crossing[int(split[1])][1] += 1 # total trial runs
        else:
             crossing[int(split[1])][1] += 1
    
        
def volume_bulk(lines):
    for line in lines:
        if re.match('\d*\_*\d*\_\d*\s\((-?\d+)(\.\d+)?,\s(-?\d+)(\.\d+)?,\s(-?\d+)(\.\d+)?',line):
            xLow = float(re.split('\(|\)|,',line)[1])
            yLow = float(re.split('\(|\)|,',line)[2])
            zLow = float(re.split('\(|\)|,',line)[3])
            xHigh = float(re.split('\(|\)|,',line)[5])
            yHigh = float(re.split('\(|\)|,',line)[6])
            zHigh = float(re.split('\(|\)|,',line)[7])
            volume = (xHigh - xLow) * (yHigh - yLow) * (zHigh - zLow)
            volumes.append(volume)
    volumes.sort()
    v1 = math.floor((len(volumes) - 1) / 2)
    v2 = math.floor(len(volumes) / 2)
    meanVolume = ((volumes[v1] + volumes[v2]) /2) * 10 ** (-30)
    return meanVolume

bulkVolume = volume_bulk(joboutLines)
    

simustep = []
for i in range(500):
    simustep.append(0)
    
for i in range(len(targetLambda)):
    if targetLambda[i] == lambdaList[1]:
        simustep[int(universe[i])] = int(step[i]) 

totalStep = 0
for i in simustep:
        if i>equilibriumSteps:
            totalStep = totalStep+i-equilibriumSteps
        else:
            totalStep = totalStep

phi = float(initialCrossing)/(timestep * 10 ** (-12) * totalStep * bulkVolume)   # time length = time step * 5fs
            
growthPossibilitys = []
overall = 1   
overalls = []   
for i in range(len(crossing)):
    if crossing[i][1] == 0:
        break
    growthPossibility = crossing[i][0]/crossing[i][1]
    growthPossibilitys.append(growthPossibility)
    overall *= growthPossibility
    overalls.append(overall)
    


nucleationRate = phi * overall


print('    %-15s %-10d'%('time steps',totalStep))
print('    %-15s %-10.2f'%('time(ns)',totalStep*timestep/1000))
for i in range(num_interface):
    print('    %-4d => %-4d     %-4d / %-7d %-10.2e' %(int(lambdaList[i+1]),int(lambdaList[i+2]),crossing[i][0],crossing[i][1],growthPossibilitys[i]))

print('************************************')

print('    %-15s %-10.2e'%('growth prob',overall))
print('    %-15s %-10.2e'%('volume',bulkVolume))
print('    %-15s %-10.2e'%('flux rate',phi))
print('    %-15s %-10.2e'%('nucleation rate',nucleationRate))
            



