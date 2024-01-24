 #-*- coding: utf-8 -*-
"""
Created on Fri Dec 18 14:53:09 2020

@author: lioneltruflandier
"""

from __future__ import division

import os.path
import re, sys
from pylab import plot, show, legend

from itertools import repeat
from numpy import size, exp, array, append, delete, ones, inf
from scipy.optimize import curve_fit
from src.read_write_files import color_list_orbital_fit
from src.read_write_files import orbital
from scipy import linalg

def gauss1(r,a0,d0):
    return d0*exp(-a0*r**2)

def gauss1_c(r,c0,a0,d0):
    return d0*exp(-a0*(r-c0)**2)

def gauss2(r,a0,d0,a1,d1):
    return d0*exp(-a0*r**2) + d1*exp(-a1*r**2)

def gauss2_c(r,a0,d0,c0,a1,d1,c1):
    return d0*exp(-a0*(r-c0)**2) + d1*exp(-a1*(r-c1)**2)
   
def gauss3(r,a0,d0,a1,d1,a2,d2):
    return d0*exp(-a0*r**2) + d1*exp(-a1*r**2) + d2*exp(-a2*r**2)

def gauss3_c(r,a0,d0,c0,a1,c1,d1,a2,d2,c2):
    return d0*exp(-a0*(r-c0)**2) + d1*exp(-a1*(r-c1)**2) + d2*exp(-a2*(r-c2)**2)

def gauss4(r,a0,d0,a1,d1,a2,d2,a3,d3):
    return d0*exp(-a0*r**2) + d1*exp(-a1*r**2) + d2*exp(-a2*r**2) + d3*exp(-a3*r**2)

def gauss4_c(r,a0,d0,c0,a1,d1,c1,a2,d2,c2,a3,d3,c3):
    return d0*exp(-a0*(r-c0)**2) + d1*exp(-a1*(r-c1)**2) + d2*exp(-a2*(r-c2)**2) + d3*exp(-a3*(r-c3)**2)

# def gauss5(r,a0,d0,c0,a1,d1,c1,a2,d2,c2,a3,d3,c3,a4,d4,c4):
#     return d0*exp(-a0*(r-c0)**2) + d1*exp(-a1*(r-c1)**2) + d2*exp(-a2*(r-c2)**2) + d3*exp(-a3*(r-c3)**2) + d4*exp(-a4*(r-c4)**2)

# def gauss6(r,a0,d0,c0,a1,d1,c1,a2,d2,c2,a3,d3,c3,a4,d4,c4,a5,d5,c5):
#     return d0*exp(-a0*(r-c0)**2) + d1*exp(-a1*(r-c1)**2) + d2*exp(-a2*(r-c2)**2) + d3*exp(-a3*(r-c3)**2) + d4*exp(-a4*(r-c4)**2) + d5*exp(-a5*(r-c5)**2)

# def gauss7(r,a0,d0,c0,a1,d1,c1,a2,d2,c2,a3,d3,c3,a4,d4,c4,a5,d5,c5,a6,d6,c6):
#     return d0*exp(-a0*(r-c0)**2) + d1*exp(-a1*(r-c1)**2) + d2*exp(-a2*(r-c2)**2) + d3*exp(-a3*(r-c3)**2) + d4*exp(-a4*(r-c4)**2) + d5*exp(-a5*(r-c5)**2) \
#     + d6*exp(-a6*(r-c6)**2)
    
# def gauss8(r,a0,d0,c0,a1,d1,c1,a2,d2,c2,a3,d3,c3,a4,d4,c4,a5,d5,c5,a6,d6,c6,a7,d7,c7)   :
#     return d0*exp(-a0*(r-c0)**2) + d1*exp(-a1*(r-c1)**2) + d2*exp(-a2*(r-c2)**2) + d3*exp(-a3*(r-c3)**2) + d4*exp(-a4*(r-c4)**2) + d5*exp(-a5*(r-c5)**2) \
#     + d6*exp(-a6*(r-c6)**2) + d7*exp(-a7*(r-c7)**2)

 
def gto_fit_orb( x, y, nG, guess, n, name, zeta, index, center, maxfev, bounds, method):
        print()        
        print('#####################################')
        print('orb[',index,'] => state '+str(n)+name+', zeta '+str(zeta))
        print('#####################################')

        #######################################################################
        # if nG is not defined exit
        #######################################################################
        if ( nG <= 0 ):
             nG = 3
             print('gto_fit_orb: default number for nG =', nG)
        
        #######################################################################
        # if nG guess not defined set it as 1 ; check the size of guess
        #######################################################################
        if ( size(guess) == 0 ):

            if (center == False):
                guess = ones(nG*3)
            else:
                guess = ones(nG*2)
        else:
            # check wether or not the provides guess has correct size ; exit otherwise
            if (center == False):
                if ( size(guess) != nG*3):
                    print('gto_fit_orb: problem with the initial guess =',size(guess),'vs. nparam expected =',nG*3)
                    sys.exit()
            else:
                if ( size(guess) != nG*2):
                    print('gto_fit_orb: problem with the initial guess =',size(guess),'vs. nparam expected =',nG*2)
                    sys.exit()
                
        #######################################################################
        # if bounds defined as empty tuple
        #######################################################################
        bounds_min = [] ; bounds_max = []
        if (size(bounds) == 0):
            for i in range(nG):
                if (center == False):
                    bounds_min.extend(3*[-inf])
                    bounds_max.extend(3*[+inf])
                else:
                    bounds_min.extend(2*[-inf])
                    bounds_max.extend(2*[+inf])
            bounds = (bounds_min,bounds_max)
        # otherwise try to use unconsistent guess (experimental)
        # else:            
        #     if ( center == False ):
        #         if ( (size(guess) % 3) != 0 ) :            
        #             if  ( size(guess) > nG*3 ):
        #                 print('gto_fit_orb: larger initial guess than expected ; shortened') 
        #                 for i in range( size(guess) - nG*2 ):
        #                     guess = delete(guess,-1)
    
        #             elif( size(guess) < nG*3 ):
        #                 print('gto_fit_orb: smaller initial guess than expected ; augmented')                         
        #                 #print guess
        #                 for i in range( nG*3 - size(guess) ):
        #                     guess = append(guess,1)

        #         else:
        #             print('gto_fit_orb: problem with the initial guess =',guess) 
        #             sys.exit()
                
        #     else:
        #         if ( (size(guess) % 2) != 0 ) :            
        #             if  ( size(guess) > nG*2 ):
        #                 print('gto_fit_orb: larger initial guess than expected ; shortened') 
        #                 for i in range( size(guess) - nG*2 ):
        #                     guess = delete(guess,-1)
    
        #             elif( size(guess) < nG*2 ):
        #                 print('gto_fit_orb: smaller initial guess than expected ; augmented')                         
        #                 #print guess
        #                 for i in range( nG*2 - size(guess) ):
        #                     guess = append(guess,1)

                #else:
                #    print('gto_fit_orb: problem with the initial guess =',guess) 
                #    sys.exit()
                
        # print guess and bounds
        print('gto_fit_orb: current initial guess =')
        print(guess)
            

        if  ( nG == 1 ):  
            
            if (center == False):
                popt, cov = curve_fit( gauss1, x, y, p0=guess, maxfev=maxfev, method=method, bounds=bounds)
                a0,d0,c0 = popt
                y_fit = gauss1_c(x,a0,d0,c0)
            else:
                popt, cov = curve_fit( gauss1, x, y, p0=guess, maxfev=maxfev, method=method, bounds=bounds)
                a0,d0 = popt
                y_fit = gauss1(x,a0,d0)
                                            
        elif( nG == 2 ):
            
            if (center == False) :                    
                popt, cov = curve_fit( gauss2_c, x, y, p0=guess, method=method, maxfev=maxfev, bounds=bounds)
                a0,d0,c0,a1,d1,c1 = popt
                y_fit = gauss2_c(x,a0,d0,c0,a1,d1,c1)
            else :
                popt, cov = curve_fit( gauss2,   x, y, p0=guess, method=method, maxfev=maxfev, bounds=bounds)
                a0,d0,a1,d1 = popt
                y_fit = gauss2(x,a0,d0,a1,d1)
            
        elif( nG == 3 ):
            
            if ( center == False ) :
                popt, cov = curve_fit( gauss3_c, x, y, p0=guess, method=method, maxfev=maxfev, bounds=bounds)
                a0,d0,c0,a1,d1,c1,a2,d2,c2 = popt
                y_fit = gauss3_c(x,a0,d0,c0,a1,d1,c1,a2,d2,c2)
            else:
                popt, cov = curve_fit( gauss3,   x, y, p0=guess, method=method, maxfev=maxfev, bounds=bounds)
                a0,d0,a1,d1,a2,d2 = popt
                y_fit = gauss3(x,a0,d0,a1,d1,a2,d2)

        elif( nG == 4 ):
            
             if(center == False) : 
                 popt, cov = curve_fit( gauss4_c, x, y, p0=guess, method=method, maxfev=maxfev, bounds=bounds) #, bounds=bounds )
                 a0,d0,c0,a1,d1,c1,a2,d2,c2,a3,d3,c3 = popt
                 y_fit = gauss4_c(x,a0,d0,c0,a1,d1,c1,a2,d2,c2,a3,d3,c3)
             else:
                 popt, cov = curve_fit( gauss4,   x, y, p0=guess, method=method, maxfev=maxfev, bounds=bounds) #bounds=bounds )
                 a0,d0,a1,d1,a2,d2,a3,d3 = popt
                 y_fit = gauss4(x,a0,d0,a1,d1,a2,d2,a3,d3)

        # elif( nG == 5 ):
        #     if(center == False):
        #         popt, cov = curve_fit( gauss5_r, x, y )
        #         a0,b0,a1,b1,a2,b2,a3,b3,a4,b4,r0,r1,r2,r3,r4 = popt
        #         y_fit = gauss5_r(x,a0,b0,a1,b1,a2,b2,a3,b3,a4,b4,r0,r1,r2,r3,r4)
        #     else:  
        #         popt, cov = curve_fit( gauss5, x, y, guess )
        #         a0,b0,a1,b1,a2,b2,a3,b3,a4,b4 = popt
        #         y_fit = gauss5(x,a0,b0,a1,b1,a2,b2,a3,b3,a4,b4)

        # elif( nG == 6 ):
        #     if(center ==False):
        #         popt, cov = curve_fit( gauss6_r, x, y )
        #         a0,b0,a1,b1,a2,b2,a3,b3,a4,b4,a5,b5,r0,r1,r2,r3,r4,r5 = popt
        #         y_fit = gauss6_r(x,a0,b0,a1,b1,a2,b2,a3,b3,a4,b4,a5,b5,r0,r1,r2,r3,r4,r5)
        #     else:
        #         popt, cov = curve_fit( gauss6, x, y, guess )
        #         a0,b0,a1,b1,a2,b2,a3,b3,a4,b4,a5,b5 = popt
        #         y_fit = gauss6(x,a0,b0,a1,b1,a2,b2,a3,b3,a4,b4,a5,b5)

        # elif( nG == 7 ):
        #     if(center == False):
        #         popt, cov = curve_fit( gauss7_r, x, y )
        #         a0,b0,a1,b1,a2,b2,a3,b3,a4,b4,a5,b5,a6,b6,r0,r1,r2,r3,r4,r5,r6 = popt
        #         y_fit = gauss7_r(x,a0,b0,a1,b1,a2,b2,a3,b3,a4,b4,a5,b5,a6,b6,r0,r1,r2,r3,r4,r5,r6)
        #     else:
        #         popt, cov = curve_fit( gauss7, x, y, guess )
        #         a0,b0,a1,b1,a2,b2,a3,b3,a4,b4,a5,b5,a6,b6 = popt
        #         y_fit = gauss7(x,a0,b0,a1,b1,a2,b2,a3,b3,a4,b4,a5,b5,a6,b6)

        # elif( nG == 8 ):
        #     if (center == False):
        #         popt, cov = curve_fit( gauss8_r, x, y )
        #         a0,b0,a1,b1,a2,b2,a3,b3,a4,b4,a5,b5,a6,b6,a7,b7,r0,r1,r2,r3,r4,r5,r6,r7 = popt
        #         y_fit = gauss8_r(x,a0,b0,a1,b1,a2,b2,a3,b3,a4,b4,a5,b5,a6,b6,a7,b7,r0,r1,r2,r3,r4,r5,r6,r7)
        #     else:
        #         popt, cov = curve_fit( gauss8, x, y, guess )
        #         a0,b0,a1,b1,a2,b2,a3,b3,a4,b4,a5,b5,a6,b6,a7,b7 = popt
        #         y_fit = gauss8(x,a0,b0,a1,b1,a2,b2,a3,b3,a4,b4,a5,b5,a6,b6,a7,b7)
      
        delta_r= x[1]-x[0]  
        chi2 = sum((y - y_fit)**2)
        
        print('state '+str(n)+name+' with '+str(nG)+' Gaussians: chi2 =', chi2,'error =',(chi2*delta_r), 'normalisation =',linalg.norm(sum(y_fit**2*delta_r)))
        
                                                                                    
        print('results')
        print(repr(popt))
        print('fit with '+str(nG)+' Gaussians: chi2 =', chi2)

        plot(x, y_fit, linestyle='--', label=str(n)+name+'-'+str(nG)+'G', color=color_list_orbital_fit(index))

        return popt, nG
