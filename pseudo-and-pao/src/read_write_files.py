# -*- coding: utf-8 -*-
"""
Created on Sat Jan 30 18:11:28 2021

@author: lioneltruflandier
"""
from __future__ import division

import os.path
import re, sys
from numpy import zeros, array, exp

#%% 
class orbital:
    def __init__(self,n,l,z,pol,pop,grid,x,y):
        self.n = n
        self.l = l
        self.z = z
        self.pol  = pol
        self.pop  = pop        
        self.grid = grid
        self.x = x
        self.y = y 
        self.y_fit = 0.0
        self.lname = ''
        self.gto_a = []
        self.gto_d = []
        self.gto_c = []
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
        #print 'i orb', i, x
        orb.append( orbital(n=n[i],l=l[i],z=z[i],pol=pol[i],pop=pop[i],grid=grid[i],x=x,y=y) )
        
        orb[i].ang_mom()

        print(n[i],orb[i].lname,' zeta = ',z[i],'polarized =',pol[i],'occ =',pop[i])    

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
def write_gto( filename_gto, symbol, orb, overwrite ):
        
    if ( overwrite == False ):
        #
        if os.path.isfile(filename_gto):
            print('write_gto:',filename_gto,' file exists, do you want to overwrite it? (Y/N)')
            res = input()
            #
            if ( 'y' in res or 'Y' in res ):
                file = open(filename_gto, 'w')
            else:
                exit()
        else:
            file = open(filename_gto, 'w')                
    else:
        print('write_gto:',filename_gto,' will be overwritten!')
        file = open(filename_gto, 'w')                
        
    norb = len(orb)
    #file = open(filename_gto, 'w')    
    file.write(symbol+'\n')
    file.write(str(norb)+'\n')
    for i in range(norb):
        #
        file.write('{:3.2}{:2d}{:6.4}{:4d}\n'.format(str(orb[i].n)+orb[i].lname, orb[i].nG, orb[i].pop, orb[i].pol))
        #file.write(str(orb[i].lname)+' '+str(nG)+'\n')    
        #file.write(str(orb[i].lname)+str(orb[i].n)+str(orb[i].l)+str(orb[i].z)+'\n')
        for j in range(orb[i].nG):
            #file.write(str(orb[i].gto_a[j])+'_'+str(orb[i].gto_d[j])+'\n')
            file.write('{:18.10f}{:18.10f}{:18.10f}\n'.format(orb[i].gto_a[j],orb[i].gto_d[j],orb[i].gto_c[j]))
            
    print('write_gto:',filename_gto,' written')      

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
def read_gto( filename_gto ):

    if os.path.isfile(filename_gto):
        print ('read_gto: file'+' '+filename_gto+' '+'exists')
    else:
        print ('read_gto: file'+' '+filename_gto+' '+'does not exist')
        exit()
        
    file = open(filename_gto, 'r')
    #
    symbol = str( file.readline() ).strip() ; print('read_gto: symbol',symbol)
    norb   = int( file.readline() ) ; print('read_gto: norb',norb)

    orb = [ ] ; zeta = 1
    for i in range(norb):
        lname, nG, pop, pol = file.readline().split()        
        #print (lname, nG, pop)
        #orb.append( orbital(n=0,l=0,z=0,pol=False,pop=0.0,grid=0.0,x=0.0,y=0.0) )
        tmp_split    = list(lname)
        orb.append( orbital(n=int(tmp_split[0]),l=0,z=zeta,pol=False,pop=float(pop),\
        grid=0.0,x=0.0,y=0.0) )
        #orb[i].n     = int(tmp_split[0])
        orb[i].lname = tmp_split[1]           
        orb[i].nG    = int(nG)
        orb[i].pol   = bool(int(pol))
        orb[i].ang_mom_reverse()
        #
        if ( i > 0 and orb[i].n == orb[i-1].n and orb[i].lname == orb[i-1].lname):
            orb[i].z = orb[i-1].z+1     
        
        print ('read_gto :', orb[i].n, orb[i].l, orb[i].lname, 'nG', orb[i].nG,\
        'zeta', orb[i].z, 'pop', orb[i].pop, 'pol', orb[i].pol)
        for j in range( orb[i].nG ):
            a, d, c = file.readline().split()   
            orb[i].gto_a.append( float(a) )
            orb[i].gto_d.append( float(d) )
            orb[i].gto_c.append( float(c) )             
            orb[i].guess.append( float(a) )
            orb[i].guess.append( float(d) )
            orb[i].guess.append( float(c) )

    orbital_kind = orbital_detection( orb )
        
    return orb
    
#%% Check consistency between curent orb and orb_guess and build a guess
def build_guess(orb, orb_guess, center):
    
    if ( len(orb) != len(orb_guess) ):
        #
        print('build_guess: orbs size are inconsistent, trying to resolve',len(orb), len(orb_guess))
        #        
    elif ( len(orb) == len(orb) ):
        #
        print('build_guess: orbs size are consistent', len(orb))    
    
    finish = False ; i = 0 ; j = 0
    #for i in range(norb):        
    while ( finish == False ):
        #
        norb = max(len(orb),len(orb_guess))
        #        
        #print('i',i,'/',len(orb_guess),'j',j,'/',len(orb))
        if ( j < len(orb) and i < len(orb_guess) ):
            #if ( orb[j].n == orb_guess[i].n and orb[j].z == orb_guess[i].z and  \
            #    orb[j].lname == orb_guess[i].lname ):
            if ( orb[j].z == orb_guess[i].z and  \
                orb[j].lname == orb_guess[i].lname ):

                #print('build_guess:', str(orb[j].n)+orb[j].lname+'-'+'zeta'+str(orb[j].z)+' recognized, guess =',orb[k].guess)
                orb[j].nG    = orb_guess[i].nG
                #orb[j].gto_a = orb_guess[i].gto_a             
                #orb[j].gto_d = orb_guess[i].gto_d
                for k in range(orb_guess[i].nG):
                    orb[j].guess.append(orb_guess[i].gto_a[k])
                    orb[j].guess.append(orb_guess[i].gto_d[k])
                #print 'j', j, orb[j].guess
                print('build_guess:', str(orb[j].n)+orb[j].lname+'-'+'zeta'+str(orb[j].z)+' recognized, guess =',orb[j].guess)

                if (orb[j].nzeta > orb_guess[i].nzeta):                
                    #print('build_guess: duplicate orb[',j,'].zeta')
                    for k in range(j+1,j+orb[j].nzeta):                                    
#                        print('build_guess: duplicate orb[',j,'].zeta for orb[',k,'].zeta')
                        #print('build_guess:', str(orb[k].n)+orb[k].lname+'-'+'zeta'+str(orb[k].z)+' duplicated, guess =',orb[k].guess)
                        orb[k].nG    = orb_guess[i].nG
                        #orb[k].gto_a = orb_guess[i].gto_a             
                        #orb[k].gto_d = orb_guess[i].gto_d
                        if (center == False):
                            for l in range(orb_guess[i].nG):
                                orb[k].guess.append(orb_guess[i].gto_a[l])
                                orb[k].guess.append(orb_guess[i].gto_d[l])
                                orb[k].guess.append(orb_guess[i].gto_c[l])
                        else:
                            for l in range(orb_guess[i].nG):
                                orb[k].guess.append(orb_guess[i].gto_a[l])
                                orb[k].guess.append(orb_guess[i].gto_d[l])
                                
                        #print 'k',k, orb[k].guess
                        print('build_guess:', str(orb[k].n)+orb[k].lname+'-'+'zeta'+str(orb[k].z)+' duplicated, guess =',orb[k].guess)

                        
                    i += 1
                    j += orb[j].nzeta

                elif( orb[j].nzeta < orb_guess[i].nzeta ):
                    
                    i += orb_guess[i].nzeta
                    j += 1

                else:
                    j += 1                                    
                    i += 1

        elif ( j < len(orb) and i == len(orb_guess) ): 
            orb[j].nG = 3
            for k in range(orb[j].nG):
                orb[j].guess.append(1.0)
                orb[j].guess.append(1.0)
            print('build_guess: create orb[',j,'].zeta')         
            #print 'j',j, orb[j].guess, 'i', i,len(orb_guess), len(orb)
            j += 1
        if ( j == len(orb) ):
            finish = True
    
    
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
    
    
    