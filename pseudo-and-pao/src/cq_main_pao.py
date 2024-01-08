# -*- coding: utf-8 -*-
"""
Created on Sat Feb  6 19:29:54 2021

@author: lioneltruflandier
"""
  
import getopt, sys
import os.path
import subprocess

from src.cq_periodic_table import element_base

# List of basis size
basis_size_list = ['minimal', 'small', 'medium', 'large']
basis_name_list = ['SZ', 'SZP', 'DZP', 'TZTP']

# List of available XC functional
xc_functional_list = ['LDA', 'PBE', 'PBEsol', 'BLYP']

def check_path(path):
    if (path[-1] != '/'):
        path = path+'/'
    return path
 


#%%
def count_line_file(filename):
    lines = 0
    for line in open(filename):
        lines += 1
    return lines
    
#%%
def check_atom( atom ):

    found = False
    #print('check_atom',atom)
    for Z, (name, symbol, ions, uncommon_ions) in element_base.items(): 
        if ( atom == symbol ):
            print(name,'found ( Z =',Z,')')
            found = True
    if ( not found ):
        print(atom,'not found, exiting!') ; exit()
        Z = 0
    return Z
    
#%%    
def main_args( ):
    
    # tmp atom list
    tmp_atom_list_name = 'tmp_atom.list'    
         
    # Remove 1st argument from the
    arguments_list = sys.argv[1:]
     
    # Options
    options = 'hx:b:l:a:g:rp:'
    
    # Defaults values
    xc = [ ] ; atom_list = [ ] ; basis_size = [ ] ; atom = [ ] 
    gto = False ; lib = False ; path='./'
    try:
        # Parsing argument
        args, vals = getopt.getopt(arguments_list, options)
         
        # checking each argument
        for arg, val in args:
     
            if arg in ('-h', '--help'):
                print ('main_args: usage')
                print ('$ python cq_test.py -x xc_functional -b basis_size -l atom_list -a atom -g')
                print ('')            
                print ('examples:')
                print ('for a list of atom contained in cq_pao.list')            
                print ('$ python cq_test.py -x PBE -b small -l cq_pao.list')
                print ('')            
                print ('for a single atom')            
                print ('$ python cq_test.py -x PBE -b small -a H')
                print ('')            
                print ('for all the atoms contained in the xc_functional directory')           
                print ('$ python cq_test.py -x PBE -b small -a all')
                print ('')
                print ('to generate GTO fit of the PAO radial part')           
                print ('$ python cq_test.py -x PBE -b small -a H -g')
                print ('')
                print ('  xc_functional available:')
                for xc in xc_functional_list:
                    print('  ',xc)
                print(' ') 
                print ('  basis_size:')
                for basis in basis_size_list:
                    print('  ',basis)            
                 
            elif arg in ('-x'):
                if (val not in xc_functional_list):
                    print('main_args:',val,'xc not available, exiting!') 
                    exit()
                else:
                    xc.append(val)
                                
            elif arg in ('-l'):
                #print (('list of atoms %s') %val)
                atom_list.append(val)
    
            elif arg in ('-a'):
                #print (('atom: %s') %val)
                atom.append(val)
                Z = check_atom( atom[0] )
    
            elif arg in ('-b'):
                if (val not in basis_size_list and val not in basis_name_list):
                    print('main_args:',val,'basis not available, exiting!') 
                    exit()
                #elif (val not in basis_name_list):
                #    print(val,'basis not available, exiting!') 
                #    exit()                    
                else:
                    if (val in basis_name_list):
                        i = 0
                        for test in basis_name_list:
                            if ( val == test):
                                val = basis_size_list[i]
                            i += 1
                            
                    basis_size.append(val)

            elif arg in ('-p'):
                print('main_args: path = ',val) 
                path = val
                    
            elif arg in ('-g'):
                gto = True         

            elif arg in ('-r'):
                lib = True         
                
    except getopt.error as err: # except try
        # output error, and return with an error code
        print (str(err)) ; print('exiting') ; exit()    

    if   ( len(xc) == 0 ):
            print('main_args: xc_functional is mandatory, exiting!') 
            exit()          
            
    elif ( (len(atom_list) or len(atom)) == 0 ):
            print('main_args: either atom_list or atom must be given, exiting!') 
            exit()            

    if( len(basis_size) == 0 ):
        print('no basis provided')
        basis_size.append('not provided')        
    
    path = check_path(path)   
    if not os.path.isdir(xc[0]):
        print ('main_args:', xc[0], 'directory is missing, exiting')
        exit()
        
    else:

        print (('main_args: xc functional= %s') % xc[0])
        print (('main_args: basis set    = %s') % basis_size[0])
        
        if ( len(atom) != 0 ):   
            atom_list = [ ] 
            file = open(tmp_atom_list_name, 'w') 
            file.write(atom[0])            
            file.close()
            print (('main_args: atom = %s') % atom[0])
            atom_list.append( atom[0] )                        
        else:
            print (('main_args: list of atoms = %s') % atom_list[0])
            if os.path.isfile(atom_list[0]):
                print('main_args:',atom_list[0],'file exists')
                n = count_line_file(atom_list[0]) 
                if ( n==0 ):
                    print ('main_args:',atom_list[0],'file empty, exiting!')
                    exit()
                else:
                    print('main_args: found',n,'elements:')
                    file = open(atom_list[0], 'r')                     
                    atom_list = [ ]                    
                    for line in file:
                        if ( line != "\n" ):                                              
                            Z = check_atom( line.strip() )
                            #if not os.path.isdir('./'+xc[0]+'/'+str(line.strip())):
                            #    print ('./'+xc[0]+'/'+str(line.strip())+'/', 'directory is missing, skip this element')
                            #    #exit()
                            #else:
                            atom_list.append( line.strip() )

                    file.close()
                    file = open(tmp_atom_list_name, 'w') 
                    for i in range(len(atom_list)):
                        file.write(str(atom_list[i])+'\n')
                    file.close()
                        
                        
#                file = open(atom_list[0], 'r') 
#                data = file.readline() ; print(data)
                
            else:
                print('main_args:',atom_list[0],'file does not exist, exiting!')
                exit()

    return path, atom_list, xc[0], basis_size[0], gto, lib
                
#%%
def main_process( path, alist, xc, bs, gto, lib ):
    #
    check_files = ['Conquest_ion_input', '.in', '.pot' ]        
    #
    cq_pao_command = 'MakeIonFiles'
    #
    CQ_ION_FILE  = 'CQ.ion'
    CQ_SPEC_FILE = 'CQ.spec'
    CQ_BLOCK_FILE= 'CQ.block'
    #
    if (path != './'):
        cq_pao_bindir = path
        if not os.path.isdir(cq_pao_bindir+'bin/') : print('main_process:',\
            cq_pao_bindir+'bin/','not found, exiting') ; exit() 
    else:
        cq_pao_bindir = '../'

    cq_pao_command = cq_pao_bindir+'bin/'+cq_pao_command
    print(cq_pao_command)

    if not os.path.isfile(cq_pao_bindir+'bin/'+cq_pao_command) : print('main_process:',\
        cq_pao_bindir+'bin/'+cq_pao_command,'not found, exiting') ; exit()        
    #
    if ( len(alist) == 0): print('main_process: no atoms to process, exiting');exit()
    #
    atom_list = [ ]
    for atom in alist:
        #
        if not os.path.isdir(xc+'/'+str(atom)):
            print('main_process: ',xc+'/'+str(atom)+'/', 'directory is missing, skip this element')
        #
        else:    
            if ( len(os.listdir(xc+'/'+str(atom)) ) == 0 ):
                print('main_process: ',xc+'/'+str(atom)+'/', 'directory is empty, skip this element')
            else:
                print('main_process: ',xc+'/'+str(atom)+'/')
                #
                atom_list.append(atom)
                #
                check_files[1]=atom+check_files[1]
                check_files[2]=atom+check_files[2]
                #print(check_files)
                #
                i = 0 ; check_files_res = [False,False,False]
                for file in os.listdir(xc+'/'+str(atom)):
                    if ( file == check_files[0] ):
                        print('main_process: ',file,'file found') ; i+=1
                        check_files_res[0] = True
                    if ( file == check_files[1] ):

                        print('main_process: ',file,'file found') ; i+=1
                        check_files_res[1] = True
                    if ( file == check_files[2] ):
                        print('main_process: ',file,'file found') ; i+=1
                        check_files_res[2] = True

                j = 0
                for res in check_files_res:                    
                    if ( res == False):
                        print('main_process: ',check_files[j],'not found')
                    j += 1

                if ( i < 3):
                    print('main_process: input file(s) are missing, exiting')
                    exit()
                        
    for atom in atom_list:
        print('')
        print('############## processing',atom,'##############')
        os.chdir(xc+'/'+str(atom)+'/')
        cq_pao_command = '../../'+cq_pao_command  
        f = open("blah.txt", "w")
        process = subprocess.Popen(cq_pao_command.split(), stdout=f)        
        output, error = process.communicate()
        if ( error): print('main_process: error in executing MakeIonFiles')
        os.chdir('../../')

     
#%%        
#if __name__ == "__main__":
#        
#   alist, xc, bs, gto, lib = main_args(  )
#   #subprocess.call(['python', 'helloworld.py', 'arg'])
#   main_process( toto, xc, bs, gto, lib )





