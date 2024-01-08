# -*- coding: utf-8 -*-
"""
Created on Fri Dec 18 14:53:09 2020

@author: lioneltruflandier
"""

from __future__ import division

import os.path
import re, sys
from pylab import plot, show, legend

from numpy import size, exp, array, append, delete
from scipy.optimize import curve_fit
from src.read_write_files import color_list_orbital_fit

def gauss1(r,a0,d0):
    return d0*exp(-a0*r**2)

def gauss2(r,a0,d0,a1,d1):
    return d0*exp(-a0*r**2) + d1*exp(-a1*r**2) 
    
def gauss3(r,a0,d0,a1,d1,a2,d2):
    return d0*exp(-a0*r**2) + d1*exp(-a1*r**2) + d2*exp(-a2*r**2)

def gauss4(r,a0,d0,a1,d1,a2,d2,a3,d3):
    return d0*exp(-a0*r**2) + d1*exp(-a1*r**2) + d2*exp(-a2*r**2) + d3*exp(-a3*r**2)

def gauss5(r,a0,d0,a1,d1,a2,d2,a3,d3,a4,d4):
    return d0*exp(-a0*r**2) + d1*exp(-a1*r**2) + d2*exp(-a2*r**2) + d3*exp(-a3*r**2) + d4*exp(-a4*r**2)

def gauss6(r,a0,d0,a1,d1,a2,d2,a3,d3,a4,d4,a5,d5):
    return d0*exp(-a0*r**2) + d1*exp(-a1*r**2) + d2*exp(-a2*r**2) + d3*exp(-a3*r**2) + d4*exp(-a4*r**2) + d5*exp(-a5*r**2)

def gauss7(r,a0,d0,a1,d1,a2,d2,a3,d3,a4,d4,a5,d5,a6,d6):
    return d0*exp(-a0*r**2) + d1*exp(-a1*r**2) + d2*exp(-a2*r**2) + d3*exp(-a3*r**2) + d4*exp(-a4*r**2) + d5*exp(-a5*r**2) \
    + d6*exp(-a6*r**2)

def gauss8(r,a0,d0,a1,d1,a2,d2,a3,d3,a4,d4,a5,d5,a6,d6,a7,d7)   :
    return d0*exp(-a0*r**2) + d1*exp(-a1*r**2) + d2*exp(-a2*r**2) + d3*exp(-a3*r**2) + d4*exp(-a4*r**2) + d5*exp(-a5*r**2) \
    + d6*exp(-a6*r**2) + d7*exp(-a7*r**2)


def gto_fit_orb( x, y, nG, guess, n, name, zeta, index ):

        print()        
        print('#####################################')
        print('orb[',index,'] => state '+str(n)+name+', zeta '+str(zeta))
        print('#####################################')

        if ( nG <= 0 ):
            nG = 3
            print('gto_fit_orb: default number for nG =',nG) 
        
        if ( size(guess) == 0 ):
            print('gto_fit_orb: proceed with no initial guess',size(guess)) 
                
            if  ( nG == 1 ):            
                popt, cov = curve_fit( gauss1, x, y )
                a0,b0 = popt
                y_fit = gauss1(x,a0,b0)
                
            elif( nG == 2 ):
                popt, cov = curve_fit( gauss2, x, y )
                a0,b0,a1,b1 = popt
                y_fit = gauss2(x,a0,b0,a1,b1)

            elif( nG == 3 ):
                popt, cov = curve_fit( gauss3, x, y )
                a0,b0,a1,b1,a2,b2 = popt
                y_fit = gauss3(x,a0,b0,a1,b1,a2,b2)
                
            elif( nG == 4 ):
                popt, cov = curve_fit( gauss4, x, y )
                a0,b0,a1,b1,a2,b2,a3,b3 = popt
                y_fit = gauss4(x,a0,b0,a1,b1,a2,b2,a3,b3)

            elif( nG == 5 ):
                popt, cov = curve_fit( gauss5, x, y )
                a0,b0,a1,b1,a2,b2,a3,b3,a4,b4 = popt
                y_fit = gauss5(x,a0,b0,a1,b1,a2,b2,a3,b3,a4,b4)

            elif( nG == 6 ):
                popt, cov = curve_fit( gauss6, x, y )
                a0,b0,a1,b1,a2,b2,a3,b3,a4,b4,a5,b5 = popt
                y_fit = gauss6(x,a0,b0,a1,b1,a2,b2,a3,b3,a4,b4,a5,b5)

            elif( nG == 7 ):
                popt, cov = curve_fit( gauss7, x, y )
                a0,b0,a1,b1,a2,b2,a3,b3,a4,b4,a5,b5,a6,b6 = popt
                y_fit = gauss7(x,a0,b0,a1,b1,a2,b2,a3,b3,a4,b4,a5,b5,a6,b6)

            elif( nG == 8 ):
                popt, cov = curve_fit( gauss8, x, y )
                a0,b0,a1,b1,a2,b2,a3,b3,a4,b4,a5,b5,a6,b6,a7,b7 = popt
                y_fit = gauss8(x,a0,b0,a1,b1,a2,b2,a3,b3,a4,b4,a5,b5,a6,b6,a7,b7)
                
            else:
                print('gto_fit_orb: no implemented for nG =',nG) 
                exit
                
                
        #print('size(guess)', (size(guess)), guess)
        if ( size(guess) != 0 ):            
            if ( (size(guess) % 2) == 0 ) :            
                if  ( size(guess) > nG*2 ):
                    print('gto_fit_orb: larger initial guess than expected ; shortened') 
                    for i in range( size(guess) - nG*2 ):
                        guess = delete(guess,-1)
    
                elif( size(guess) < nG*2 ):
                    print('gto_fit_orb: smaller initial guess than expected ; augmented')                         
                    #print guess
                    for i in range( nG*2 - size(guess) ):
                        guess = append(guess,1)

            else:
                print('gto_fit_orb: problem with the initial guess =',guess) 
                exit()
                
            print('gto_fit_orb: current initial guess =')
            print(guess)

            if  ( nG == 1 ):            
                popt, cov = curve_fit( gauss1, x, y, guess )
                a0,b0 = popt
                y_fit = gauss1(x,a0,b0)
                                                
            elif( nG == 2 ):
                popt, cov = curve_fit( gauss2, x, y, guess )
                a0,b0,a1,b1 = popt
                y_fit = gauss2(x,a0,b0,a1,b1)
                
            elif( nG == 3 ):
                popt, cov = curve_fit( gauss3, x, y, guess )
                a0,b0,a1,b1,a2,b2 = popt
                y_fit = gauss3(x,a0,b0,a1,b1,a2,b2)

            elif( nG == 4 ):
                popt, cov = curve_fit( gauss4, x, y, guess )
                a0,b0,a1,b1,a2,b2,a3,b3 = popt
                y_fit = gauss4(x,a0,b0,a1,b1,a2,b2,a3,b3)

            elif( nG == 5 ):
                popt, cov = curve_fit( gauss5, x, y, guess )
                a0,b0,a1,b1,a2,b2,a3,b3,a4,b4 = popt
                y_fit = gauss5(x,a0,b0,a1,b1,a2,b2,a3,b3,a4,b4)

            elif( nG == 6 ):
                popt, cov = curve_fit( gauss6, x, y, guess )
                a0,b0,a1,b1,a2,b2,a3,b3,a4,b4,a5,b5 = popt
                y_fit = gauss6(x,a0,b0,a1,b1,a2,b2,a3,b3,a4,b4,a5,b5)

            elif( nG == 7 ):
                popt, cov = curve_fit( gauss7, x, y, guess )
                a0,b0,a1,b1,a2,b2,a3,b3,a4,b4,a5,b5,a6,b6 = popt
                y_fit = gauss7(x,a0,b0,a1,b1,a2,b2,a3,b3,a4,b4,a5,b5,a6,b6)

            elif( nG == 8 ):
                popt, cov = curve_fit( gauss8, x, y, guess )
                a0,b0,a1,b1,a2,b2,a3,b3,a4,b4,a5,b5,a6,b6,a7,b7 = popt
                y_fit = gauss8(x,a0,b0,a1,b1,a2,b2,a3,b3,a4,b4,a5,b5,a6,b6,a7,b7)


            
        chi2 = sum((y - y_fit)**2)
        #print('state '+str(n)+name+' with '+str(nG)+' Gaussians: chi2 =', chi2)
        print('results')
        print(repr(popt))
        print('fit with '+str(nG)+' Gaussians: chi2 =', chi2)

        #plot(x, y_fit, linestyle='--', label=str(n)+name+'-'+str(nG)+'G', color=color_list_orbital_fit(index))

        return popt, nG
