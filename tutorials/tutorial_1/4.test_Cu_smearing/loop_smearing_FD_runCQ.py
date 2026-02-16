#!/bin/python
import os
import shutil
import subprocess

os.environ['OMP_NUM_THREADS'] = '1'

# CQ root directory
CQROOT='/Users/lioneltruflandier/CONQUEST-tutorials/bin'
# CQ mpi command
CQMPI='mpirun -np'
# CQ nprocess
CQNUM='2'
# CQ command
CQCMD=CQMPI+' '+CQNUM+' '+CQROOT+'/Conquest'
# Ion files
ions = ['Cu.ion']
# Loop over each temperature in Ha 
values = [0.0001, 0.0005, 0.001, 0.01]
# basename
base = 'fd'

print('Conquest cmd',CQCMD)
# Loop over each cutoff energy value
for val in values:
    dirname=str(f"{val:.4f}")
    path=base+'_'+dirname
    if not os.path.isdir(path):
        os.makedirs(path)
                
    shutil.copy2('./coords.dat',path)
    shutil.copy2('./Cu.ion',path)        
    shutil.copy2('./Conquest_input',path)

    os.chdir(path)    
    with open('Conquest_input','a') as f:
        f.write('\n\ndiag.smearingtype 0')
        f.write('\ndiag.kt '+str(val))
         
    print('Conquest is running in', path,'...')
    subprocess.run(CQCMD+'> Conquest_out', shell=True)
    
    os.chdir('./../')
    
