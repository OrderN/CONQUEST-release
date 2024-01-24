# -*- coding: utf-8 -*-
"""
Created on Sat Jan 30 18:11:28 2021

@author: lioneltruflandier
"""
from __future__ import division

import os.path
import re, sys
from numpy import zeros, array, exp
import matplotlib.pylab as plt


#%% 
def color_list_orbital(n):
    color_list = ['black','darkred','darkgreen','darkblue','darkcyan','darkmagenta','gold',\
    'darkorange',]
    if ( n < 8):
        return color_list[n]
    else:
        return 'black'
    
def color_list_orbital_fit(n):
    color_list = ['gray','lightcoral','lightgreen','lightblue','lightcyan','magenta',\
    'lightyellow','moccasin']
    if ( n < 8 ):
        return color_list[n]
    else:
        return 'blue'

#%% Plot radial part of the orbitals ##########################################
def plot_PAO( orb, filename_gto_write, save=False ):
    #
    figure0 = plt.figure(0)
    plt.plot(orb[0].x,zeros(len(orb[0].x)),color='black',linestyle='dotted',linewidth=1.0)
    #
    for i in range(len(orb)):
        # Plot original PAO orbital
        tmp_label = str(orb[i].n)+orb[i].lname+'-'+str(orb[i].z)+'$\zeta$'+'-PAO'
        plt.plot(orb[i].x, orb[i].y,label=tmp_label, \
            color=color_list_orbital(i) )
    
    plt.xlabel('$r$ in au')
    plt.ylabel('radial part $\phi(r)$')
    plt.legend()
    plt.show()
    if ( save == True ): 
        figure0.savefig( filename_gto_write+'_radial.pdf', dpi=300, format='pdf', transparent=True,\
        bbox_inches='tight', frameon=None)    
    

#%% Plot radial part of the orbitals ##########################################
def plot_GTO_PAO( orb, filename_gto_write, save=False ):
    #
    figure1 = plt.figure(1)
    plt.plot(orb[0].x,zeros(len(orb[0].x)),color='black',linestyle='dotted',linewidth=1.0)
    #
    for i in range(len(orb)):
        # Plot original PAO orbital
        tmp_label = str(orb[i].n)+orb[i].lname+'-'+str(orb[i].z)+'$\zeta$'+'-PAO'
        plt.plot(orb[i].x, orb[i].y,label=tmp_label, \
            color=color_list_orbital(i) )
        # Plot the PAO-nG orbital
        tmp_label = str(orb[i].n)+orb[i].lname+'-'+str(orb[i].z)+'$\zeta$'+'-'+str(orb[i].nG)+'G'         
        plt.plot(orb[i].x, orb[i].y_fit,label=tmp_label,\
            color=color_list_orbital_fit(i),linestyle='--')
    
    plt.xlabel('$r$ in au')
    plt.ylabel('radial part $\phi(r)$')
    plt.legend()
    plt.show()
    if ( save == True ): 
        figure1.savefig( filename_gto_write+'_radial.pdf', dpi=300, format='pdf', transparent=True,\
        bbox_inches='tight')    
    

#%% Plot PAO-GTO deviation ####################################################
def plot_GTO_PAO_dev( orb, filename_gto_write, save=False ):
    #
    figure2 = plt.figure(2)
    plt.plot(orb[0].x,zeros(len(orb[0].x)),color='black',linestyle='dotted',linewidth=1.0)
    #
    for i in range(len(orb)):
        # Compute deviation from PAO
        y_dev = orb[i].y - orb[i].y_fit
        # Plot    
        tmp_label = str(orb[i].n)+orb[i].lname+'-'+str(orb[i].z)+'$\zeta$'
        plt.plot(orb[i].x, y_dev, label=tmp_label, color=color_list_orbital(i))

    plt.xlabel('$r$ in au')
    plt.ylabel('deviation $\phi^{PAO}(r)-\phi^{PAO}_{nG}(r)$')
    plt.legend()
    plt.show()
    if ( save == True ):        
        figure2.savefig( filename_gto_write+'_deviation.pdf', dpi=300, format='pdf', transparent=True,\
        bbox_inches='tight')    

#%% For each orb plot the GTO primitives ######################################
def plot_GTO_prim( orb, filename_gto_write, save=False):
    #
    for i in range(len(orb)):
        # One figure for each orbital
        figure = plt.figure(i+3)
        plt.plot(orb[0].x,zeros(len(orb[0].x)),color='black',linestyle='dotted',linewidth=1.0)
        # Plot each Gaussian function
        
        
        for j in range(orb[i].nG):
            y = orb[i].gto_d[j]*exp( -orb[i].gto_a[j]*(array(orb[i].x) - orb[i].gto_c[j])**2 )
            plt.plot(orb[i].x, y,label='$a=%12.6f$, $d=%12.6f$, $c=%12.6f$' %(orb[i].gto_a[j],orb[i].gto_d[j],orb[i].gto_c[j]))
    
        # Plot the PAO radial function
        tmp_label = str(orb[i].n)+orb[i].lname+'-'+str(orb[i].z)+'$\zeta$'+'-PAO'
        plt.plot(orb[i].x, orb[i].y,label=tmp_label, \
        linestyle='-', color=color_list_orbital(i))
        # Plot the GTO radial function as the sum over primitives
        tmp_label = str(orb[i].n)+orb[i].lname+'-'+str(orb[i].z)+'$\zeta$'+'-PAO-'+str(orb[i].nG)+'G'
        plt.plot(orb[i].x, orb[i].y_fit,label=tmp_label, \
        linestyle='dotted',color=color_list_orbital_fit(i))

        plt.xlabel('$r$ in au')
        plt.ylabel('radial part')
        plt.title('$\phi^{PAO}_{nG}(r)=\sum_i^n d_i \exp(-a_i r^2)$')
        plt.legend()
        plt.show()
        if ( save ==  True ):
            figure.savefig( filename_gto_write+'_'+str(orb[i].n)+orb[i].lname+'.pdf',\
                            dpi=300, format='pdf', transparent=True,\
                            bbox_inches='tight')           
        
