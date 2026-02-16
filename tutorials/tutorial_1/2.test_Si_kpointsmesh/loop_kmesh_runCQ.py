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
ions = ['Si.ion']
# List of cutoff energy values
values=[1, 2, 4, 6]
# basename
base = 'kp'

print('Conquest cmd',CQCMD)
# Loop over each cutoff energy value
for val in values:
    dirname=str(val)+'x'+str(val)+'x'+str(val)
    path=base+'_'+dirname
    if not os.path.isdir(path):
        os.makedirs(path)
                
    shutil.copy2('./coords.dat',path)
    shutil.copy2('./Si.ion',path)        
    shutil.copy2('./Conquest_input',path)

    os.chdir(path)    
    with open('Conquest_input','a') as f:
        f.write('\n\ndiag.mpmesh  True')
        f.write('\ndiag.mpmeshx '+str(val))
        f.write('\ndiag.mpmeshy '+str(val))
        f.write('\ndiag.mpmeshz '+str(val))
         
    print('Conquest is running in', path,'...')
    subprocess.run(CQCMD+'> Conquest_out', shell=True)
    
    os.chdir('./../')
    
