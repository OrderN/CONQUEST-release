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
ions = ['SZ', 'SZP', 'DZP', 'TZTP']
# basename
base = 'bs'

print('Conquest cmd:',CQCMD)
# Loop over each cutoff energy value
for val in ions:
    dirname=str(val)
    path=base+'_'+dirname
    if not os.path.isdir(path):
        os.makedirs(path)
                
    shutil.copy2('./coords.dat',path)
    shutil.copy2('./Conquest_input',path)
    shutil.copy2('./Si_'+str(val)+'.ion',path)

    os.chdir(path)    
    with open('Conquest_input','a') as f:

        f.write('\n\n%block ChemicalSpeciesLabel')
        f.write('\n1   28.085 Si_'+str(val))
        f.write('\n%endblock ChemicalSpeciesLabel')
         
    print('Conquest is running in', path,'...')
    subprocess.run(CQCMD+'> Conquest_out', shell=True)
    
    os.chdir('./../')
    
