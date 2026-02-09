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
import numpy as np

#%% Plot radial part of the orbitals ##########################################
def plot_PAO( orb, normalize, filename ):
    #
    figure0 = plt.figure(0)
    plt.plot(orb[0].x,zeros(len(orb[0].x)),color='black',linestyle='dotted',linewidth=1.0)
    #
    for i in range(len(orb)):
        # Compute the norm
        if ( normalize ):
           norm_orb = np.trapz(orb[i].y,dx=orb[i].x[1]-orb[i].x[0])
        else:
           norm_orb = 1.0 
        # Plot original PAO orbital
        tmp_label = str(orb[i].n)+orb[i].lname+'-'+str(orb[i].z)+'$\zeta$'+'-PAO'+' rc = '+str(round(orb[i].rc,4))
        plt.plot(orb[i].x, np.array(orb[i].y)/norm_orb,label=tmp_label)
        #integ2= np.trapz(orb[i].y/norm_orb,dx=orb[i].x[1]-orb[i].x[0])   
        #print(integ,integ2)
 
    plt.xlabel('$r$ in au')
    plt.ylabel('radial part $\phi(r)$')
    plt.legend()
    plt.show()
    figure0.savefig( filename+'_radial.pdf', dpi=300, format='pdf', transparent=True,\
        bbox_inches='tight')    
    

#%% 
class orbital:
    def __init__(self,n,l,z,pol,pop,grid,rc,x,y):
        self.n = n
        self.l = l
        self.z = z
        self.pol  = pol
        self.pop  = pop        
        self.grid = grid
        self.rc = rc
        self.x = x
        self.y = y 
        self.y_fit = 0.0
        self.lname = ''
        self.gto_a = []
        self.gto_d = []
        self.guess = []
        self.nzeta = 0
        self.nG    = 0
        
    def ang_mom(self):
        if   ( self.l == 0 ):
            self.lname = 'S'
        elif ( self.l == 1 ):
            self.lname = 'P'
        elif ( self.l == 2 ):
            self.lname = 'D'
        elif ( self.l == 3 ):
            self.lname = 'F'
        else:            
            self.lname = 'X'

    def ang_mom_reverse(self):
        if   ( self.lname == 'S' ):
            self.l = 0
        elif ( self.lname == 'P' ):
            self.l = 1
        elif ( self.lname == 'D' ):
            self.l = 2
        elif ( self.lname == 'F' ):
            self.l = 3
        else:            
            self.l =-1

    def y_fit_gto(self):
        #self.y_fit = zeros(int(self.grid))
        for i in range(self.nG):
            self.y_fit = self.y_fit + \
            self.gto_d[i]*exp( -self.gto_a[i]*array(self.x)**2 )
        
              
#%% 
def read_ion( filename_ion ):             
    
    if os.path.isfile(filename_ion):
        print ('read_ion: file'+' '+filename_ion+' '+'exists')
    else:
        print ('read_ion: file'+' '+filename_ion+' '+'does not exist')
        exit()
    # search for data: symbol, valence electron and lmax ##########################
    file = open(filename_ion, 'r')
            
    pattern1 = 'Element'
    pattern2 = 'Valence'
    pattern3 = 'Lmax for basis'
    
    for line in file:    
        if re.search(pattern1, line):
            res1   = line.split()
            symbol = res1[0]
        if re.search(pattern2, line):
            res2 = line.split()
            nele = float(res2[0])
        if re.search(pattern3, line):
            res3 = line.split()
            norb = int(res3[1])

    # search for pao data and extract radial part #############################
    n    = [ ] # principal quantum number
    l    = [ ] # angular momentum
    z    = [ ] # zeta
    pol  = [ ] # polarization function ?
    pop  = [ ] # population
    grid = [ ] # grid definition of each orbital
    rc   = [ ] # cutoff radius
    #
    x_orb= []
    y_orb= []
    #
    file.seek(0) # rewind position at file
    #
    pattern4='orbital l, n, z' ; numl = [ ]
    #
    i = 0
    file = open(filename_ion, 'r')
    for num, line in enumerate(file, 1):
        if (pattern4 in line):
            print('read_ion: found orbital at line:', num)
            numl.append(num)
            i += 1                
            
    if (i != norb):
        print('read_ion: error number of orbitals =', i, norb)
        sys.exit('exit')
    
    file.seek(0) # rewind position at file
    #
    orb = [ ]
    #
    for i in range(norb): # loop over the norb orbitals
        #
        #print(numl[i])
        for j in range(numl[i]):
            line = file.readline()
            res  = line.split()
            #
        n.append( int(res[1]) )
        l.append( int(res[0]) )
        z.append( int(res[2]) )    
        if ( int(res[3]) == 0):
            pol.append( False )
        else:
            pol.append( True )            
        pop.append( float(res[4]) )    
        #
        line2 = file.readline()
        res2  = line2.split()
        grid.append( res2 )
        rc.append( float(res2[2]) )
        #     
        x = [ ] ; y = [ ]
        for j in range( int(grid[i][0]) ):
            line = file.readline()
            res  = line.split()
            x.append( float(res[0]) )
            y.append( float(res[1]) )
        #print 'i orb', i, x
        orb.append( orbital(n=n[i],l=l[i],z=z[i],pol=pol[i],pop=pop[i],grid=grid[i],rc=rc[i],x=x,y=y) )
        
        orb[i].ang_mom()

        print(n[i],orb[i].lname,' zeta = ',z[i],'polarized =',pol[i],'occ =',pop[i], 'rc = %.3f'%rc[i])    

        file.seek(0)

    file.close()
                       
    basis_kind = orbital_detection( orb )            
            
    return symbol, norb, orb, basis_kind
    
    
#%%    
def kind_zeta(n):
    if ( n == 1):
        return 'SZ'
    elif ( n == 2):
        return 'DZ'
    elif ( n == 3):
        return 'TZ'
    elif ( n == 4):
        return 'QZ'
    elif ( n == 5):
        return '5Z'
    elif ( n == 6):
        return '6Z'
    else:
        return 'XZ'
        
#%%            
def kind_pol(n):
    if ( n == 1):
        return 'P'
    elif ( n == 2):
        return 'DP'
    elif ( n == 3):
        return 'TP'
    elif ( n == 4):
        return 'QP'
    elif ( n == 5):
        return '5P'
    elif ( n == 6):
        return '6P'
    elif ( n == 0):
        return ''
    else:
        return 'XP'       
        
#%%
def orbital_detection( orb ):    
    #
    norb = len(orb)
    # Find out basis kind
    max_z = 1
    for i in range(norb-1): # loop over the norb orbitals        
        if ( orb[i].z > orb[i+1].z ):
            max_z = orb[i].z       
    #print max_z    
    
    tmp_kind_pol = [ ] ; max_kind_zeta = 0 ; max_kind_pol = 0

    if ( max_z > 1 ):                
        
        for i in range(norb): # loop over the norb orbitals
            #
            if ( orb[i].z == max_z and orb[i].pop == 0.0 and orb[i].pol == False ):    
                #   
                tmp_kind = [ ] ; orb[i].nzeta = max_z
                for j in range(1,max_z):
                    #
                    if (orb[i].n == orb[i-j].n and orb[i].lname == orb[i-j].lname ):
                        tmp_kind.append(True)
                        orb[i-j].nzeta = max_z
                        
                if ( all(tmp_kind) ):
                    #print 'tmp_kind', tmp_kind
                    #tmp_kind_zeta = kind_zeta( len(tmp_kind)+1 )
                    max_kind_zeta = max(max_kind_zeta,len(tmp_kind)+1)
                    #print('state',n[i],orb[i].lname,'is',tmp_kind_zeta)

                #elif ( len(tmp_kind) == 0 ):
                #    tmp_kind_zeta = kind_zeta( 1 )
                #    print('state',n[i],orb[i].lname,'is',tmp_kind_zeta)
                
            if ( orb[i].z == max_z and orb[i].pop == 0.0 and orb[i].pol == True ):    
                #
                tmp_kind = [ ] ; orb[i].nzeta = max_z
                for j in range(1,max_z):
                    #
                    if (orb[i].n == orb[i-j].n and orb[i].lname == orb[i-j].lname ):
                        tmp_kind.append(True)
                        #tmp_kind_pol.append(1)
                        orb[i-j].nzeta = max_z
                    
                if ( all(tmp_kind) ):
                    #print 'tmp_kind', tmp_kind
                    #tmp_kind_zeta = kind_zeta( len(tmp_kind)+1 )   
                    #tmp_kind_pol  = kind_pol ( len(tmp_kind)+1 ) 
                    max_kind_pol  = max(max_kind_pol,len(tmp_kind)+1)  
                    #print('state',n[i],orb[i].lname,'is',tmp_kind_zeta,' and is polarisation')
             
    tmp_kind_pol = [ ] #; max_kind_zeta = 0 ; max_kind_pol = 0
    for i in range(norb): # loop over the norb orbitals
        #        
        #print  i, orb[i].nzeta, orb[i].pop, orb[i].pol
        if ( orb[i].nzeta == 0 and orb[i].pop > 0.0 and orb[i].pol == False ):    
            orb[i].nzeta = 1
            #tmp_kind_zeta = kind_zeta( orb[i].nzeta )             
            max_kind_zeta = max(max_kind_zeta,orb[i].nzeta)            
            #print('state',n[i],orb[i].lname,'is',tmp_kind_zeta)

        if ( orb[i].nzeta == 0 and orb[i].pop == 0.0 and orb[i].pol == True ):
            #print 'here'
            orb[i].nzeta = 1
            #tmp_kind_zeta = kind_zeta( orb[i].nzeta )                         
            #tmp_kind_pol.append(True)
            max_kind_zeta = max(max_kind_zeta,orb[i].nzeta)
            max_kind_pol  = max(max_kind_pol,orb[i].nzeta) #; print max_kind_pol
            #print('state',n[i],orb[i].lname,'is',tmp_kind_zeta,' and is polarisation')
            
    #print max_kind_zeta, max_kind_pol
    final_kind_zeta = kind_zeta( max_kind_zeta )
    final_kind_pol  = kind_pol ( max_kind_pol )    
    #
    orbital_kind = final_kind_zeta+final_kind_pol
    print(orbital_kind+' detected')   
    
    #     
    return orbital_kind

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
    
    
    
