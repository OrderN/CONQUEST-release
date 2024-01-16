# -*- coding: utf-8 -*-
"""
Created on Fri Dec 18 14:53:09 2020

@author: lioneltruflandier
"""

import os.path
import re, sys
sys.path.append('../../')

from pylab import plot, show, legend
import matplotlib.pylab as plt
from numpy import size, exp, array, zeros

from src.gto_fit          import gauss3, gto_fit_orb
from src.read_write_files import read_ion, orbital
from src.read_write_files import write_gto, read_gto
from src.plot_orbitals    import color_list_orbital, color_list_orbital_fit
from src.plot_orbitals    import plot_GTO_PAO, plot_GTO_PAO_dev, plot_GTO_prim, plot_PAO
from src.read_write_files import build_guess


#%% Edit file names ###########################################################

# file generated by  MakeIonFiles to extract obitals data 
filename_ion='NeCQ.ion'         
# file to write the GTO fit results formatted for Conquest
filename_gto_write='NeCQ.gto_TZTP_3G' 
# file to read in order to build an initial guess (not mandatory)
filename_gto_read ='NeCQ.gto_TZTP_3G_guess' 

#%% Read the ion file genrerated by MakeIonFiles
# return the chemical symbol, the number orbitals @(norb) 
# and the orbitals as the class orb 
# (cf. read_ion.py file in ../../src)  
symbol, norb, orb, kind = read_ion( filename_ion )

#%% Plot PAO only #############################################################
#plot_PAO( orb, filename_gto_write )

#%% Define the initial guess for each orb. ####################################
# (not mandatory)
#
orb_guess = read_gto(filename_gto_read)
build_guess( orb, orb_guess)
#
#for i in range(norb):
#   orb[i].guess = orb_guess[i].guess  
#
#orb[0].guess = array([0.340524815, 0.4463843004, 3.2450659895, 0.4537962082, 0.1019009042, 0.2253588404, 0.9568838539, 0.4357005884, -0.01, 0.02])
#orb[1].guess = array([ 0.1,-1,5,-1,1,1])

#%% Define the number of Gaussian primitives for each orb. ####################
# (not mandatory, default is 3)
#
#for i in range(norb):
#    orb[i].nG.append(3)
#orb[0].nG = 1
#orb[0].nG = 4
#orb[1].nG = 4
#
#%% Fit the radial part and plot ##############################################
for i in range(norb):
    x = array(orb[i].x)
    y = array(orb[i].y)
    # 
    # GTO fit ; plots of the GTO radial part is done in gto_fit_orb
    param, nG = gto_fit_orb( x, y, orb[i].nG, orb[i].guess, orb[i].n, orb[i].lname, orb[i].z, i )
    #
    # If the number of Gaussian has not been setup ; 
    # gto_fit_orb set it to 3 by default
    if ( orb[i].nG == 0 ):
        orb[i].nG = nG
    #
    # Store GTO exponent (a) and coefficient (d)
    for j in range(0,orb[i].nG*2,2):
        orb[i].gto_a.append( param[j]   )      
        orb[i].gto_d.append( param[j+1] )      
    # Generate and store GTO fit results on the radial grid
    orb[i].y_fit_gto()

#%% Write GTO file ############################################################
write_gto( filename_gto_write, symbol, orb, True )
    
#%% Plot radial part of the orbitals ##########################################
plot_GTO_PAO( orb, filename_gto_write )

#%% Plot PAO-GTO deviation ####################################################
plot_GTO_PAO_dev( orb, filename_gto_write )

#%% For each orb plot the GTO primitives ######################################
plot_GTO_prim( orb, filename_gto_write )

