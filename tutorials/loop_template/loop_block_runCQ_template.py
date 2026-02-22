#!/bin/python
import os
import shutil
import subprocess
import time
import re

subprocess.run(['bash', './../environment_variables.sh'])

# CQ root directory
CQROOT = os.getenv('CQROOT')
if (not CQROOT): 
    CQROOT='/Users/lioneltruflandier/CONQUEST-tutorials/bin'
# CQ mpi command
CQMPI = os.getenv('CQMPI')
if (not CQMPI): 
    CQMPI='mpirun -np'
# CQ nprocess
CQNUM = os.getenv('CQNUM')
if (not CQNUM): 
    CQNUM='2'
# CQ nthreads not greater than 4
OMP_NUM_THREADS = os.getenv('OMP_NUM_THREADS')
if ( not  OMP_NUM_THREADS ):
    os.environ['OMP_NUM_THREADS'] = '1'
# Setup Conquest command
CQCMD = CQMPI+' '+CQNUM+' '+CQROOT+'/Conquest'
# List of files to copy
files = ['Si.ion', 'coords.dat', 'Conquest_template', 'block.txt']
# List of cutoff energy values in Ha
values=[1, 2, 4, 6]
# Basename
base = 'kp'
# Pattern to search and replace
pattern = re.compile(r"\$val")
# Start timer
start = time.time()

print('Conquest cmd:',CQCMD)
# Loop over each cutoff energy value
for val in values:
    dirname = str(val)+str(val)+str(val)
    path = base+'_'+dirname
    if not os.path.isdir(path):
        os.makedirs(path)
        
    for file in files:
        shutil.copy2(file,path)

    os.chdir(path)    
    shutil.copy2('Conquest_template', 'Conquest_input')
    with open('block.txt', 'r') as src, open('Conquest_input', 'a') as dst:
        for line in src:
            new_line = pattern.sub(str(val), line)
            dst.write(new_line)
            
    print('Conquest is running in', path,'...')
    subprocess.run(CQCMD+'> Conquest_out', shell=True)
    
    os.chdir('./../')

# End timer
end = time.time()

# Print elapsed time
print('total of %.3f seconds elapsed'%(end -start))   
