# -*- coding: utf-8 -*-
"""
Created on Sat Jan 30 18:11:28 2021

@author: lioneltruflandier
"""
from __future__ import division

import os.path
import re, sys

class orbital:
    def __init__(self,n,l,z,grid,x,y):
        self.n = n
        self.l = l
        self.z = z
        self.grid = grid
        self.x = x
        self.y = y 
        self.lname = ''
        self.gto_a = []
        self.gto_d = []
        self.guess = []
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

def read_data_ion( filename_ion ):
    
    if os.path.isfile(filename_ion):
        print ('File'+' '+filename_ion+' '+'exists')
    else:
        print ('File'+' '+filename_ion+' '+'does not exist')
        exit
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

    # search for pao data and extract radial part #################################
    n    = [ ] # principal quantum number
    l    = [ ] # angular momentum
    z    = [ ] # number of zeta
    grid = [ ] # grid definition of each orbital
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
            print('found orbital at line:', num)
            numl.append(num)
            i += 1                
            
    if (i != norb):
        print('error number of orbitals =', i, norb)
        sys.exit('exit')
    
    file.seek(0) # rewind position at file
    #
    orb = [ ]
    #
    for i in range(norb): # loop over the norb orbitals
        #
        for j in xrange(numl[i]):
            line = file.readline()
            res  = line.split()
            #
        n.append( int(res[1]) )
        l.append( int(res[0]) )
        z.append( int(res[2]) )    
        #
        line = file.readline()
        res  = line.split()
        grid.append( res )
        #     
        x = [ ] ; y = [ ]
        for j in range( int(grid[i][0]) ):
            line = file.readline()
            res  = line.split()
            x.append( float(res[0]) )
            y.append( float(res[1]) )
    
        orb.append( orbital(n=n[i],l=l[i],z=z[i],grid=grid[i],x=x,y=y) )
        
        orb[i].ang_mom()

        print(n[i],orb[i].lname,' zeta = ',z[i])    
        orb[i].ang_mom()
        file.seek(0)

    file.close()
    
    return symbol, norb, orb
