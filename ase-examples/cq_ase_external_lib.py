#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 29 15:11:52 2022

@author: lioneltruflandier
"""

from numpy import amax, amin, matmul, array, trapz, size, exp, zeros, floor
from numpy import empty

import ase.utils

from ase.geometry import orthorhombic, is_orthorhombic
from ase.constraints import FixAtoms, FixScaled
from ase.units import Hartree, Rydberg, Bohr
from ase import Atoms

from spglib import standardize_cell, get_spacegroup

from inspect import currentframe, getframeinfo

class error(Exception):
    """Base class for exceptions in this module
    """
    pass

class warning(Warning):
    """Base class for warning in this module
    """
    pass

class conquest_err(error):
    """Exceptions related to Conquest I/O

    Attributes
    ----------
    message : explanation of the error
    """
    def __init__(self, message):
        self.message = message

class conquest_warn(warning):

    def __init__(self, message):
        self.message = message
    
        frameinfo = getframeinfo(currentframe().f_back)

        print('## ConquestWarning ##')
        print('file:', frameinfo.filename, ', line:', frameinfo.lineno)
        print('>>>> ', message)


###############################################################################
#
# 
#
###############################################################################
def print_cell_data(atoms, verbose):
    cell = atoms.get_cell()
    cell_param = atoms.get_cell_lengths_and_angles()
    print('Bravais lattice:')
    print('  ', cell.get_bravais_lattice())
    print('Space group:')
    print('  ', get_spacegroup(atoms, symprec=1e-5))
    print('Cell parameters (Ang. and degree):')
    print('   a = {:14.6f}'.format(cell_param[0]))
    print('   b = {:14.6f}'.format(cell_param[1]))
    print('   c = {:14.6f}'.format(cell_param[2]))
    print('   alpha = {:8.4f}'.format(cell_param[3]))
    print('   beta  = {:8.4f}'.format(cell_param[4]))
    print('   gamma = {:8.4f}'.format(cell_param[5]))
    #
    if (verbose):
        atom_positions = atoms.get_positions()
        atom_names = atoms.get_chemical_symbols()
        print('Cartesian atomic positions (Ang.)')
        for i in range(len(atoms)):
            # print formated forces
            print('  {}{:14.6f}{:14.6f}{:14.6f}'.format(atom_names[i],
                                                        atom_positions[i, 0],
                                                        atom_positions[i, 1],
                                                        atom_positions[i, 2]))


###############################################################################
#
# 
#
###############################################################################
def read_conquest(fileobj, fractional=True, atomic_order=[]):
    """
    Read CONQUEST structure files.
    Returns atoms object
    """
    if atomic_order == []:
        raise conquest_err("Atomic order must be specified as a list.")

    if isinstance(fileobj, str):
        fileobj = open(fileobj)

    # The first 3 lines are always the unit cell vectors in Bohr
    lines = []
    for i in range(3):
        lines.append(fileobj.readline().split())
    natoms = int(fileobj.readline().split()[0])
    # CONQUEST always uses orthorhombic cells
    cell = [float(lines[0][0]) * Bohr, float(lines[1][1]) * Bohr,
            float(lines[2][2]) * Bohr]
    # This is for constraining certain atoms
    moveflags = empty((natoms, 3), dtype=bool)
    con = {'t': True, 'f': False}
    n = 0
    species = []
    positions = []
    for line in fileobj.readlines():
        x, y, z, speciesindex, cx, cy, cz = line.split()
        speciesindex = int(speciesindex) - 1
        if fractional:
            position = [float(x) * cell[0], float(y) * cell[1],
                        float(z) * cell[2]]
        else:
            position = [float(x) * Bohr, float(y) * Bohr, float(z) * Bohr]

        moveflags[n] = [not con[cx.lower()], not con[cy.lower()],
                        not con[cz.lower()]]
        n += 1
        species.append(atomic_order[speciesindex])
        positions.append(position)

    atoms = Atoms(symbols=species, positions=positions, cell=cell)
    constraints = []
    indices = []
    for ind, moveflag in enumerate(moveflags):
        if moveflag.any() and not moveflag.all():
            constraints.append(FixScaled(atoms.get_cell(), ind, moveflag))
        elif moveflag.all():
            indices.append(ind)
    if indices:
        constraints.append(FixAtoms(indices))
    if constraints:
        atoms.set_constraint(constraints)
    atoms.set_pbc((True, True, True))
    return atoms

###############################################################################
#
# 
#
###############################################################################
def conquest_orthorhombic_check(atoms, verbose):
    # get cell from ASE Atoms object
    cell = atoms.get_cell()
    # test cell and if not orthorhombic standardize and check again
    if not is_orthorhombic(cell):
        # Print cell type before modification
        conquest_warn('Current cell not orthorhombic:')
        print_cell_data(atoms, verbose=verbose)
        # If not orthorhombic standardize and check

        atoms_ = ase.utils.atoms_to_spglib_cell(atoms)
        cell, scaled_positions, numbers = standardize_cell(atoms_,
                                                           to_primitive=False,
                                                           no_idealize=False,
                                                           symprec=1e-5)
        # New atoms fractional positions
        atoms = Atoms(numbers, positions=scaled_positions,
                      pbc=[True, True, True])
        # Set new cell
        atoms.set_cell(cell, scale_atoms=True)
        # Test again and stop if not orthorhombic
        orthorhombic(cell)
        # Print new cell data
        print('\n')
        conquest_warn('New orthorhombic cell:')
        print_cell_data(atoms, verbose=verbose)
        # Take care "atoms" may have been modified !
    return atoms

def get_gapwind( calc, kpath, n_occ ):

    if ( n_occ / floor(n_occ) != 1.0 ):
        print('n_occ = {} is not an integer'.format(n_occ))
        quit()
    else:
        eig_homo = -1e9   
        eig_lumo = +1e9   
        for kpt in range( len(kpath.kpts) ):
        
            eig_kpt = calc.get_eigenvalues(kpt=kpt)
        
            if ( eig_homo < eig_kpt[int(n_occ[0])-1] ):
                eig_homo  = eig_kpt[int(n_occ[0])-1]
                kpt_homo  = kpt
                        
            if ( eig_lumo > eig_kpt[int(n_occ[0])] ):
                eig_lumo  = eig_kpt[int(n_occ[0])]
                kpt_lumo  = kpt

    return [eig_homo,kpt_homo], [eig_lumo,kpt_lumo]
        #print_occ(eig,fermi_energy,0.0,'aufbau',str(kpt))

###############################################################################
#
# Given file_xsf and number of atoms add ANIMSTEPS n and backup last struct
#
###############################################################################
def cq_postproc_xsf(file_xsf,n_atoms):
    
    f = open(file_xsf, 'r')
    n = len(f.readlines())  # number of lines in traj_file
    f.close()
    n_len = 8 + n_atoms     # length of the repeating unit = 8 + number of atoms
    
    if ( n % n_len ) !=0 :
        f = open(file_xsf, 'r')
        line_test = f.readline()
        if ( 'ANIMSTEPS' in line_test ): 
            print('WARNING: {} already post-processed...'.format(file_xsf))
            tmp = line_test.split()
            n_animxsf = int(tmp[1])

        else:
            print('WARNING: {}/{} probably corrupted ; stop processing...'.format(n,n_len))
            return
            
    # Prepend 'ANIMSTEPS n_animxsf' at the beginning of traj_file
    # n_animxsf = number of optimisation steps in traj_file 
    else:
        n_animxsf = n // n_len
        with open(file_xsf, 'r+') as file: 
            originalContent = file.read() 
            file.seek(0)              
            file.write('ANIMSTEPS {}\n'.format(n_animxsf))
            file.write(originalContent) 
        f.close()
            
    # Extract last structure from traj_file
    # assumed that 'ANIMSTEPS n_animxsf' at the first line of traj_file
    n_last = 1 + n_len*(n_animxsf - 1)
    #print(n_last,n_len,n_atoms,n_animxsf)
    f     = open(file_xsf, 'r')
    lines = f.readlines()[n_last:n_last+n_len] # read last structure
    #
    #print(lines)
    cell_out    = zeros((3,3)) # to store cell
    coords_out  = []           # to store oordinates in Ang.
    symbols_out = ''           # to store atomic symbols (concatenate)
    i = 0
    for line in lines:
        if ( 6 > i > 2 ): 
            tmp = line.split()
            #print(tmp,'1')
            cell_out[i-3,0] = tmp[0]
            cell_out[i-3,1] = tmp[1]
            cell_out[i-3,2] = tmp[2]
            
        elif( i > 7  ):
            tmp         = line.split()
            #print(tmp,'2')
            symbols_out = symbols_out + tmp[0]
            coords_out.append( (float(tmp[1]),float(tmp[2]),float(tmp[3])) )
        i += 1

    return Atoms(symbols=symbols_out,cell=cell_out,positions=coords_out,pbc=[True,True,True]) 
        
###############################################################################
#
# Given eigenvalues, Fermi level and method output occupation numbers
#
###############################################################################
def print_occ(eigenvalues,fermi,kT,method,label):

    occ_num = [ ]
    
    if (method == 'aufbau'):
       for eig in eigenvalues:
            if ( eig <= fermi):
                occ_num.append( 1.0 )
            else:
                occ_num.append( 0.0 )

    elif (method == 'fermi-dirac'):
        for eig in eigenvalues:
            tmp = (eig - fermi)/Rydberg/kT
            if   ( tmp >  400 ):
                occ_num.append( 0.0 )
            elif ( tmp < -400 ):
                occ_num.append( 1.0 )
            else:
                occ_num.append( 1.0/( exp(tmp) + 1.0 ) )
    else:        
        print('{}'.format('Error: method unrocognized... '))
    
    print()
    print('{}'.format('Band energies at the {} point'.format(label)))
    print('{}'.format('# energy in eV (occupation number)'))
    
    occ_tot = 0.0
    for eig, occ in zip(reversed(eigenvalues),reversed(occ_num)):
        print('{:14.5f} ({:5.3f})'.format(eig,occ*2))
        occ_tot = occ_tot + occ*2
    print('{} {:5.3f})'.format('(total'.rjust(21),occ_tot)) 


###############################################################################
#
# Given QEspresso output (ASE Atoms object) print data
#
###############################################################################
def qe_print_output_data(struct,verbose):
    # results are extracted from struct.calc.output !
    dft_energy   = struct.get_potential_energy()
    #fermi_energy = struct.calc.output.get_fermi_level()
    atom_names   = struct.get_chemical_symbols()
    atom_forces  = struct.get_forces()
    stress_tens  = struct.get_stress()
    stress_tens_comp = ['xx', 'yy', 'zz', 'yz', 'xz', 'xy']
    
    print('Total electronic energy = {:14.6f} eV'.format(dft_energy))
    print('Atoms kinetic energy    = {:14.6f} eV'.format(struct.get_kinetic_energy())) # of the atoms (relevant only with MD)
    print('Total energy = {:14.6f} eV'.format(struct.get_total_energy())) # total electronic + atoms-kinetic energy
    #print('Fermi energy = {:14.6f} eV'.format(fermi_energy)) # Fermi level
    print()
    if ( verbose > 0):
        print('Atomic forces (eV/Ang)')
        for i in range(len(struct)):
            print('  {}{:12.6f}{:12.6f}{:12.6f}'.format(atom_names[i],atom_forces[i,0],atom_forces[i,1],atom_forces[i,2]))
        print()
        
    print('Cell stress tensor components (eV/Ang**3)')
    for i in range(5):
        print('  {} = {:11.6f}'.format(stress_tens_comp[i],stress_tens[i]))

###############################################################################
#
# Given Conquest output print data
#
###############################################################################
def cq_print_output_data(struct,verbose):
    # results are extracted from struct.calc.output !
    dft_energy   = struct.calc.results['energy']
    #fermi_energy = struct.calc.output.get_fermi_level()
    atom_names   = struct.get_chemical_symbols()
    atom_forces  = struct.calc.results['forces']
    stress_tens  = struct.calc.results['stress']
    stress_tens_comp = ['xx', 'yy', 'zz', 'yz', 'xz', 'xy']
    
    print('Total electronic energy = {:14.6f} eV'.format(dft_energy))
    print('Atoms kinetic energy    = {:14.6f} eV'.format(struct.get_kinetic_energy())) # of the atoms (relevant only with MD)
    print('Total energy = {:14.6f} eV'.format(dft_energy+struct.get_kinetic_energy())) # total electronic + atoms-kinetic energy
    #print('Fermi energy = {:14.6f} eV'.format(fermi_energy)) # Fermi level
    print()
    if ( verbose > 0):
        print('Atomic forces (eV/Ang)')
        for i in range(len(struct)):
            print('  {}{:12.6f}{:12.6f}{:12.6f}'.format(atom_names[i],atom_forces[i,0],atom_forces[i,1],atom_forces[i,2]))
        print()
        
    print('Cell stress tensor components (eV/Ang**3)')
    for i in range(5):
        print('  {} = {:11.6f}'.format(stress_tens_comp[i],stress_tens[i]))
        
###############################################################################
#
# Given ASE Atoms object print data
#
###############################################################################
def print_struct_data(struct,verbose):
    
    atom_names = struct.get_chemical_symbols()
    atom_positions_cart = struct.get_positions(wrap=False)
    atom_positions_cell = struct.get_scaled_positions(wrap=False)
    
    cell_param = struct.get_cell_lengths_and_angles()
    cell_tens  = struct.get_cell()

    if ( verbose > 0):

        print('Cartesian atomic positions (Ang.)')
        for i in range(len(struct)):
            # print formated cartesian
            print('  {:2}{:14.6f}{:14.6f}{:14.6f}'.format(atom_names[i],
                atom_positions_cart[i,0],atom_positions_cart[i,1],atom_positions_cart[i,2]))
        
        print()
        print('Fractional atomic positions (Adim.)')
        for i in range(len(struct)):
            # print formated fractional
            print('  {:2}{:14.6f}{:14.6f}{:14.6f}'.format(atom_names[i],
                atom_positions_cell[i,0],atom_positions_cell[i,1],atom_positions_cell[i,2]))
        
        # Check the Cartesian positions are recovered from the cell tensor
        print()
        print('Cell tensor (Ang.)')
        for i in range(len(cell_tens)):
            print('    {:14.6f}{:14.6f}{:14.6f}'.format(cell_tens[i,0],\
                                                        cell_tens[i,1],\
                                                        cell_tens[i,2]))
        print()
        print('Cartesian atomic positions (Ang.) from cell tensor')
        for i in range(len(struct)):
            # WARNING: xyz -> row
            tmp_positions = matmul(atom_positions_cell[i,:],cell_tens)
            print('  {:2}{:14.6f}{:14.6f}{:14.6f}'.format(atom_names[i],\
                    tmp_positions[0],tmp_positions[1],tmp_positions[2]))
    
    print()
    print('Cell parameters (Ang. and degree)')
    print('  a = {:14.6f}'.format(cell_param[0]))
    print('  b = {:14.6f}'.format(cell_param[1]))
    print('  c = {:14.6f}'.format(cell_param[2]))
    print('  alpha = {:8.4f}'.format(cell_param[3]))
    print('  beta  = {:8.4f}'.format(cell_param[4]))
    print('  gamma = {:8.4f}'.format(cell_param[5]))