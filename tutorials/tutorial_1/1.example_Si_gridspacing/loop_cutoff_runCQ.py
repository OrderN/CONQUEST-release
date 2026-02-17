#!/bin/python
import os
import shutil
import subprocess
import time

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
# List of cutoff energy values in Ha
values=[50, 75, 100, 120, 150]
# basename
base = 'gs'
# Start timer
start = time.time()

print('Conquest cmd:',CQCMD)
# Loop over each cutoff energy value
for val in values:
    dirname=str(val).rjust(3,'0')
    path=base+'_'+dirname
    if not os.path.isdir(path):
        os.makedirs(path)
        #print('directory %s created'%path)
                
    shutil.copy2('./coords.dat',path)
    shutil.copy2('./Si.ion',path)        
    shutil.copy2('./Conquest_input',path)

    os.chdir(path)    
    with open('Conquest_input','a') as f:
        f.write('\n\ngrid.gridcutoff '+str(val))
        
    print('Conquest is running in', path,'...')
    subprocess.run(CQCMD+'> Conquest_out', shell=True)
    
    os.chdir('./../')

# End timer
end = time.time()

# Print elapsed time
print('total of %.3f seconds elapsed'%(end -start))   
