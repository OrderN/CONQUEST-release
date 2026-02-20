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
# List of MSSF range in Bohr 
values=[1.0, 5.0, 8.0, 11.0, 14.0]
# basename
base = 'msr'
# Start timer
start = time.time()

print('Conquest cmd:',CQCMD)
# Loop over each cutoff energy value
for val in values:
    dirname=str("{:0>4.1f}".format(val))
    path=base+'_'+dirname
    if not os.path.isdir(path):
        os.makedirs(path)
        #print('directory %s created'%path)
                
    shutil.copy2('./coords.dat',path)
    shutil.copy2('./Si.ion',path)        
    shutil.copy2('./Conquest_input',path)

    os.chdir(path)    
    with open('Conquest_input','a') as f:
        f.write('\n\n%block Si') 
        f.write('\nAtom.NumberOfSupports 4')
        f.write('\nAtom.MultisiteRange '+str(val))
        f.write('\nAtom.LFDRange '+str(val))
        f.write('\n%endblock')
        
    print('Conquest is running in', path,'...')
    subprocess.run(CQCMD+'> Conquest_out', shell=True)
    
    os.chdir('./../')

# End timer
end = time.time()

# Print elapsed time
print('total of %.3f seconds elapsed'%(end -start))   
