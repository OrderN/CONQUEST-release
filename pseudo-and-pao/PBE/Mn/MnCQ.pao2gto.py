# -*- coding: utf-8 -*-
"""
Created on Fri Dec 18 14:53:09 2020

@author: lioneltruflandier
"""

#import os.path
import re, sys
sys.path.append('../../')

from pylab import plot, show, legend
import matplotlib.pylab as plt
from numpy import size, exp, array, zeros
from scipy.optimize import curve_fit

from src.gto_fit            import gauss3, gto_fit_orb
from src.read_write_files   import read_ion, orbital
from src.read_write_files   import write_gto, read_gto
from src.plot_orbitals      import color_list_orbital, color_list_orbital_fit
from src.plot_orbitals      import plot_GTO_PAO, plot_GTO_PAO_dev, plot_GTO_prim, plot_PAO
from src.read_write_files   import build_guess
from scipy.optimize         import minimize 

#%% Edit file names ###########################################################

# file generated by  MakeIonFiles to extract obitals data 
filename_ion='MnCQ.ion'         
# file to write the GTO fit results formatted for Conquest
filename_gto_write='MnCQ.gto_SZ_3G' 
# file to read in order to build an initial guess (not mandatory)
#filename_gto_read ='MnCQ.gto_SZP_3G_guess' 

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
#center = True
#orb_guess = read_gto(filename_gto_read)
#build_guess( orb, orb_guess, center =center)
#
#for i in range(norb):
#   orb[i].guess = orb_guess[i].guess  
orb[0].guess = array([0.159486,0.598452,2.0,-1.0,0.1,0.1])
#orb[1].guess = array([0.3341085034,0.2805108907,4.1100000000,1.0297932964,1.4588060470,2.0776528365])
#orb[2].guess = array([1.0,0.389613,0.389613,0.054788,0.1,0.1])
#orb[3].guess = array([3.155396,0.423250,0.501348,-0.663974,0.031197,0.198817])
#orb[4].guess = array([3.155396,0.423250,0.501348,-0.663974,0.031197,0.198817])
#orb[5].guess = array([3.155396,0.423250,0.501348,-0.663974,0.031197,0.198817])
#orb[6].guess = array([3.155396,0.423250,0.501348,-0.663974,0.031197,0.198817])
#orb[7].guess = array([3.155396,0.423250,0.501348,-0.663974,0.031197,0.198817])


#orb[9].guess = array([1.333977,-3.808716,0.120180,0.102885,1.256987,0.932507])
#orb[10].guess = array([1.333977,-3.808716,0.120180,0.102885,1.256987,0.932507])

#orb[0].bounds=([-1,-1,-6,-6,-1,-6],[1,4,6,3,6,6])
orb[1].bounds=([-1,-1,-6,-6,-1,-6],[1,2,6,3,6,6])
#%% Define the number of Gaussian primitives for each orb. ####################
# (not mandatory, default is 3)
#
#for i in range(norb):
#    orb[i].nG = 2
#orb[0].nG = 2
orb[1].nG = 2
orb[2].nG = 2
orb[3].nG = 2
#orb[4].nG = 3
#orb[5].nG = 3
#orb[6].nG = 3
# orb[7].nG = 4
#orb[8].nG = 3
#
#%% Fit the radial part and plot ##############################################
for i in range(norb):
    x = array(orb[i].x)
    y = array(orb[i].y)
    print(orb[i].nG)
    # 
    # GTO fit ; plots of the GTO radial part is done in gto_fit_orb
    center = True
    #bounds = ([-6,-6,-6,-6],[1,2,3,4]) #Bounds for a,d and c
    param, nG = gto_fit_orb( x, y, orb[i].nG, orb[i].guess, orb[i].n, orb[i].lname, orb[i].z, i, 
                            center=center, maxfev=90000, bounds=orb[i].bounds, method='trf')
    
    #param, nG = gto_fit_orb( x, y, orb[i].nG, orb[i].guess, orb[i].n, orb[i].lname, orb[i].z, i, center = center)

    #
    # If the number of Gaussian has not been setup ; 
    # gto_fit_orb set it to 3 by default
    if ( orb[i].nG == 0 ):
        orb[i].nG = nG
    #
    # Store GTO exponent (a) and coefficient (d)
    
    if (center == False):
        for j in range(0,orb[i].nG*3,3):
            orb[i].gto_a.append( param[j]   )      
            orb[i].gto_d.append( param[j+1] )
            orb[i].gto_c.append( param[j+2] )
    else:
        for j in range(0,orb[i].nG*2,2):
            orb[i].gto_a.append( param[j]   )      
            orb[i].gto_d.append( param[j+1] )
            orb[i].gto_c.append( 0.0 )
         
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

