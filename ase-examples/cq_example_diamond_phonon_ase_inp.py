#!/usr/bin/env python
# coding: utf-8

# ### Example for phonon calculation of diamond using finite-displacement with ASE+Conquest
import os
import sys

from ase.build import bulk
from ase.spacegroup import crystal
from ase.phonons import Phonons

from ase.calculators.conquest import Conquest

# #### Directory for storing calculation files
working_directory = 'cq_example_diamond_phonon_ase'

# #### Generate diamond cell ; enforced to be conventional by using `cubic=True`
atoms = bulk('C', crystalstructure='diamond', a=3.572117687, cubic=True)

# #### Get arguments from command line when use as a python script
for i in range(len(sys.argv)):
    print('argument:', i, 'value:', sys.argv[i])
#
CQ_NPROC = int(sys.argv[1])
CQ_ROOT  = sys.argv[2]

# #### Define Conquest environment variables for ASE
#CQ_NPROC= 8
#CQ_ROOT = '/Users/lioneltruflandier/CONQUEST-develop'

os.environ['ASE_CONQUEST_COMMAND'] = 'mpirun -np '+str(CQ_NPROC)+' '+CQ_ROOT+'/bin/Conquest'
os.environ['CQ_PP_PATH']           = CQ_ROOT+'/pseudo-and-pao/'
os.environ['CQ_GEN_BASIS_CMD']     = CQ_ROOT+'/tools/BasisGeneration/MakeIonFiles'

# #### Define Conquest PAO basis
basis = {'C' : { 'gen_basis' : False}}

# #### Setup Conquest as a calculator
calc = Conquest(basis=basis,
                directory    = working_directory,
                scf_tolerance= 1e-8, # default 1e-6
                grid_cutoff  =  120) # default 100 Ha

# #### Setup Phonon calculation with small displacement
N      =  4
ph_dir = 'phonon'
ph     =  Phonons(atoms, calc, supercell=(N, N, N), delta=0.05, name=working_directory+'/'+ph_dir)

# #### Run Phonon calculation
ph.run()

