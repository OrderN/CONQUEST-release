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
print('Conquest cmd:',CQCMD)
# basename
base = 'a'

# Reference lattice constant (in angstroms)
a0_ang=3.5856
ang_to_bohr=1.88973

# Scaling factors: from -5% to +7% around a0
scales=[0.95, 0.98, 1.00, 1.02, 1.05, 1.07]
# Loop over each cutoff energy value
for s in scales:
    a_bohr = a0_ang*ang_to_bohr*s
    a_ang  = a0_ang*s
    cell = [
    str(f"{a_bohr:8.6f}") +' '+str(f"{0:8.6f}")      +' '+str(f"{0:8.6f}")+'\n',
    str(f"{0:8.6f}")      +' '+str(f"{a_bohr:8.6f}") +' '+str(f"{0:8.6f}")+'\n',
    str(f"{0:8.6f}")      +' '+str(f"{0:8.6f}")      +' '+str(f"{a_bohr:8.6f}") +'\n',
    ]

    dirname=str(f"{a_ang:4.3f}")
    path=base+'_'+dirname
    if not os.path.isdir(path):
        os.makedirs(path)
                
    shutil.copy2('./coords.dat',path)
    shutil.copy2('./Conquest_input',path)
    shutil.copy2('./C.ion',path)

    os.chdir(path)    
    with open('coords.dat', 'r') as f:
        pos = f.readlines()            
        
    with open('coords.dat', 'w') as f:
        f.writelines(cell)
        f.writelines(pos)    
         
    print('Conquest is running in', path,'...')
    subprocess.run(CQCMD+'> Conquest_out', shell=True)
    
    os.chdir('./../')
    
