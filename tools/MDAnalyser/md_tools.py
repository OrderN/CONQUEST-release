#!/usr/local/bin/python3

import scipy as sp
import numpy as np
import os
import sys
import glob
import matplotlib.pyplot as plt
from scipy.linalg import norm
from scipy.integrate import cumtrapz
from scipy.signal import correlate
from scipy.fft import fft, fftshift
from scipy.constants import physical_constants
from scipy import histogram
from math import ceil, pi
from md_frame import Frame
from pdb import set_trace


######################################################
# Conversions and global threshold #
bohr2ang = 0.529177249
small = 1.0e-3
###################################
# Function used in multiple class #
###################################
def autocorr(x, y=None):
  """Autocorrelation function"""
  if y.any():
    result = correlate(x, y, mode='full')
  else:
    result = correlate(x, x, mode='full')
  return result[result.size // 2:]

def diff_mic(pos1, pos2, cell):
  """Minimum image convention relative vector (orthorhombic cell only)"""
  diff = pos2 - pos1
  for i in range(3):
    diff[i] -= round(diff[i]/cell[i])*cell[i]
  return diff

def center_cell(pos1,pos2,cell):
    """Minimum image convention relative vector (orthorhombic cell only)"""
    dx = pos2 - pos1
    for i in range(3):
      if (dx[i] > cell[i]*0.5):
        dx[i] = dx[i]-cell[i]
      if (dx[i] <= -cell[i]*0.5):
        dx[i] = dx[i] + cell[i]
    return dx

def disp_mic_npt(pos1, pos2, cell1, cell2):
  """MIC displacement when cell dimensions change"""
  disp = pos2/cell2 - pos1/cell1
  for i in range(3):
    disp[i] -= round(disp[i]/cell2[i])
  return disp

def linear_regression_diffusion(x,a,b):
  """Function of a linear function"""
  return a*x + b

def Logical_to_string(logic):
  """"Function converting logical to readable strings for dynpro."""
  if logic:
    return "T"
  if not logic:
    return "F"

def Atom_selector(nspec,species):
    """ Function to choose the central atom to plot"""
    # Display all the atoms
    print("###########################")
    print("Species you can choose:")
    for ispec in range(nspec):
        print(species[ispec+1])
    selection = 0
    print("###########################")
    # Choosing the central atom
    while selection == 0:
      choosen_atom1 = input("Which species do you wan to take? Type 'all' to show every datas. ")
      if(choosen_atom1 == "All" or choosen_atom1 == "all" or choosen_atom1 == "ALL" or choosen_atom1 == "0"):
        Central_atom_to_print = "all"
        selection = 1
        print("###########################")
        return selection, Central_atom_to_print
      for ispec in range(nspec):
        if(str(species[ispec+1]) == choosen_atom1):
          Central_atom_to_print = choosen_atom1
          selection = 2
          print("###########################")
          return selection, Central_atom_to_print
      if(selection == 0):
        print("The atom entered isn't found. Please, enter the data again.")

def logic_selector(opt):
  """Function that will return a logic operator in function of an input."""
  print("Do you want to see all the "+str(opt)+"(s) graphs during the process? ")
  print("y or n")
  show_opt = input()
  if (show_opt == "y" or show_opt == "Y" or show_opt == "Yes" 
  or show_opt == "YES" or show_opt == "1" or show_opt == "T"
  or show_opt == "True" or show_opt == "TRUE"):
    print("-----------------------------------------------")
    print("The graphic(s) will appear between each steps.")
    print("-----------------------------------------------")
    return True
  else:
    return False

def Unknown_atom(Atomic_name):
  """Verify the input of the selected atom."""
  working_directory          = "./datas_dynpro/"
  if ( not os.path.isdir(working_directory) ):
    os.makedirs(working_directory)
    print("Create the working directory for the input files.")
  Mendel_table = [
            "H", "He", "Li", "Be",  "B",  "C",  "N",  "O",  "F", "Ne", 
            "Na", "Mg", "Al", "Si",  "P",  "S", "Cl", "Ar",  "K", "Ca", 
            "Sc", "Ti",  "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", 
            "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr",  "Y", "Zr", 
            "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "St", 
            "Sb", "Te",  "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", 
            "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", 
            "Lu", "Hf", "Ta",  "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", 
            "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", 
            "Pa",  "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", 
            "Md", "No", "Lr"
      ]
  try:
    Atomic_number = Mendel_table.index(Atomic_name)+1
  except:
    while True:
      print("The name of the central element isn't valid")
      print("Enter the name of the central atom as a placeholder.")
      print("Valid atoms:")
      print(Mendel_table)
      Placeholder = input("-> ")
      old_Atomic_name = Atomic_name
      os.system("cp ./Conquest_input ./Conquest_input_dynpro")
      f = open("./Conquest_input_dynpro","r")
      line_copy       = np.array([])
      we_are_in_block = False
      effect = False
      for line in f:
          words   = line.split()
          temp_line = []
          if len(words) >= 1:
            if words[0] == "%block":
              we_are_in_block = True
            if words[0] == "%endblock":
              we_are_in_block = False
          if len(words) >= 3:
            if words[2] == Atomic_name and we_are_in_block:
              temp_line.insert(0,words[:2])
              temp_line[0].insert(2,Placeholder)
              line_copy = np.append(line_copy,temp_line[0])
              line_copy = np.append(line_copy,"\n")
              effect = True
          if len(words) >= 2 and len(words) <= 3:
            if words[1] == Atomic_name and we_are_in_block:
              temp_line.insert(0,words[:1])
              temp_line[0].insert(2,Placeholder)
              line_copy = np.append(line_copy,temp_line[0])
              line_copy = np.append(line_copy,"\n")
              effect = True
          if not effect:
            line_copy = np.append(line_copy,words)
            line_copy = np.append(line_copy,"\n")
          else:
            effect = False
      f.close()

      f = open("./Conquest_input_dynpro","w")
      index = 0
      for line in line_copy:
        index +=1
        if line == "\n" or line_copy[index] == "\n":
          f.write(line)
        elif line == "%block":
          f.write(line+" ")
        else:
          f.write(line+"   ")
      Atomic_name = Placeholder
      f.close()
      try:
        Atomic_number = Mendel_table.index(Atomic_name)+1
      except:
        continue
      else:
        break
    os.system("mv ./Conquest_input_dynpro "+working_directory+"Conquest_input")
    print("Modification of the document with success.")
    return Atomic_name, old_Atomic_name, Atomic_number 
    # The function will output the name of the atom that need to change and it's new name.
  else:
    print("The atomic number of your chosen atom is: ", Atomic_number)
    return Atomic_name, Atomic_name, Atomic_number
  

######################################################

class Dynpro_gestion:
  """Class that will prepare and launch the dynpro program."""

  def __init__(self, nspec, rcut, binwidth, species, 
              species_count,Starting_point,Ending_point,Number_of_sample,jump,Atomic_name,r_min,Temperature,
              Prefix,Input_dynpro_name,md_frame_name,md_stats_name,lcenter,ct_atom,
              RDF_data,VACF_data,Sample_data, substi,Atom_number,Old_atom_name):
    # Global datas
    self.nspec = nspec                  # Number of different species.
    self.species = species              # Name of all the species present in the cell.
    self.rcut = rcut                    # Maximum distance taked into account for the RDF.
    self.binwidth = binwidth            # Size of the histogram step.
    self.nbins = ceil(rcut/binwidth)+1  # Number of histogram steps needed, it will be truncated to the superior.
    self.spec_count = species_count     # Number of atom per species.
    # For sampling
    self.Starting_point   = float(Starting_point)        # The starting point for the sampling, is 0 by default.
    self.Ending_point     = float(Ending_point)          # The ending point for the sampling, is 0 by default.
    self.Number_of_sample = float(Number_of_sample)      # Number of necessary samples is 0 by default.
    # For RDF and VACF
    self.jump             = int(jump)                   # Number of steps skiped during the printing of the result with conquest.
    self.Atomic_name      = str(Atomic_name)            # Central atom chosen (initialize it into a function?)
    self.Atomic_name_old  = str(Old_atom_name)          # Original name if the Conquest input was edited for Dynpro.
    self.Atomic_number    = int(Atom_number)            # Atomic number of the selected atom
    self.r_min            = float(r_min)                # Minimal radius for the RDF
    self.nr               = ceil(rcut/binwidth)+1       # Number of steps to form the histogram
    constant_for_the_time = (1e-15)/(float(physical_constants["atomic unit of time"][0])) # The time selected is in (fs), it will be converted in atomic unit for dynpro.
    self.time_step        = self.jump*constant_for_the_time # MD time step
    self.temperature      = float(Temperature)           # Temperature of the simulation in the input.log
    # For naming
    self.Prefix           = str(Prefix)                # Prefix of the MD data files.
    self.Input_dynpro_name= str(Input_dynpro_name)     # Name of the dynpro file.
    self.md_frame_name    = str(md_frame_name)         # Name of the md.frames file.
    self.md_stats_name    = str(md_stats_name)         # Name of the md.stats file.
    # Centering
    self.center           = lcenter               # Do we center the cell around an atom?
    self.ct_atom          = int(ct_atom)               # ID of the centering atom. If its 0, it will be taken automatically.
    # Selections
    self.RDF              = RDF_data              # Do we want to do the RDF?
    self.VACF             = VACF_data             # Do we want to do the VACF (PSD method is included)?
    self.Sample           = Sample_data           # Do we sample the datas?
    #
    self.RDF_data         = Logical_to_string(self.RDF)
    self.VACF_data        = Logical_to_string(self.VACF)
    self.Sample_data      = Logical_to_string(self.Sample)
    # Substitution
    self.subti            = substi                # Substitution option
    self.subtitution      = Logical_to_string(self.subti)
    self.Real_atom_name   = np.array([])
    # For debugging
    self.debug = False
    # Database
    power2table = [
      1,2,4,8,16,32,64,128,256,512,1024
    ]
    
    # Evaluation of sp_n, the number of sample to be taken
    if self.Ending_point > self.Starting_point and self.Sample:
      i_last = 0
      for i in power2table:
        if (i >= self.Ending_point):
          self.sp_n = i_last
          break
        if (i == 1024):
          self.sp_n = i
          break
        i_last = i

    # Directories selection
    self.working_directory          = "./datas_dynpro/"       # Directory for the inputs for the dynpro
    self.working_directory_scratch  = "./md_files_dynpro/"    # Directory for the inputs with md_files for the dynpro
    self.working_directory_output_  = "./Outputs_dynpro_"     # Directory for the output
    self.working_directory_output   = "./Outputs_dynpro_/"     # Directory for the output
    if ( not os.path.isdir(self.working_directory) ):
      os.makedirs(self.working_directory)
      print("Create the working directory for the input files.")
    if ( not os.path.isdir(self.working_directory_scratch) ):
      os.makedirs(self.working_directory_scratch)
      print("Create the working directory for the MD files.")
    if ( not os.path.isdir(self.working_directory_output) ):
      os.makedirs(self.working_directory_output)
      print("Create the Output directory.")

  def substitution(self,Input_Atom,ID_atom_replace):
    """Function that do the substitution of a selected atom"""
    if self.subti:
      # If we have the substitution we initialise the new atomic name
      self.substi_Atom_name = Input_Atom
      #
      if not os.path.isfile(self.working_directory+"sub.struct_perm.lk"):
        print("The file 'sub.struct_perm.lk' is needed in the working directory.")
        print("Please, run at least once dynpro to generate the file before doing a substitution.")
        print("*-----------Job cancelled-----------*")
        sys.exit()
      
      Mendel_table = [
                "H", "He", "Li", "Be",  "B",  "C",  "N",  "O",  "F", "Ne", 
                "Na", "Mg", "Al", "Si",  "P",  "S", "Cl", "Ar",  "K", "Ca", 
                "Sc", "Ti",  "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", 
                "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr",  "Y", "Zr", 
                "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "St", 
                "Sb", "Te",  "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", 
                "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", 
                "Lu", "Hf", "Ta",  "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", 
                "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", 
                "Pa",  "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", 
                "Md", "No", "Lr"
          ]
      try:
        self.sub_Atomic_number = Mendel_table.index(Input_Atom)+1
      except:
        while True:
          print("The name of the element of substitute isn't valid")
          print("Enter the name of the central atom as a placeholder.")
          print("Valid atoms:")
          print(Mendel_table)
          Input_Atom = input("-> ")
          try:
            self.sub_Atomic_number = Mendel_table.index(Input_Atom)+1
          except:
            continue
          else:
            self.substi_Atom_name = Input_Atom
            break
      # The ID is verified and the Atom name is verified. We can replace one atom.
      line_copy = np.array([])
      f = open(self.working_directory+"sub.struct_perm.lk","r")
      for line in f:
        words   = line.split()
        temp_line = []
        if len(words) < 6:
          line_copy = np.append(line_copy, words)
          line_copy = np.append(line_copy, "\n")
          continue
        if words[5] == str(ID_atom_replace):
          temp_line.insert(0,words[1:4])
          temp_line[0].insert(0,Input_Atom)
          temp_line[0].insert(4,str(self.sub_Atomic_number))
          temp_line[0].insert(5, words[5])
          line_copy = np.append(line_copy,temp_line[0])
          line_copy = np.append(line_copy,"\n")
        else:
          line_copy = np.append(line_copy, words)
          line_copy = np.append(line_copy, "\n")
      f.close()

      f = open(self.working_directory+"sub.struct_perm.in","w")
      index = 0
      reset = 0
      space = " "
      for line in line_copy:
        index +=1
        # The first lines
        if index == 1:
          f.write(space*10+line)
        if index == 3:
          f.write(space*11+line)
        if line == "\n":
          f.write(line)
          reset = 2
          
        # Document layout
        if index > 3:
          if reset == 2 and line != "\n":
              f.write(space+line)
              f.write(space*(6-len(line)))
              reset += 1
              continue
          if reset >= 3 and reset < 6:
            if (len(line_copy[index]) == 19 or len(line_copy[index]) == 24):
              f.write(line)
              f.write(space*(25-len(line)))
            elif(len(line_copy[index]) == 23 or len(line_copy[index]) == 24):
              f.write(line)
              f.write(space*(26-len(line)))
            elif (len(line_copy[index]) >= 1 and  len(line_copy[index]) <= 3 ):
              f.write(line)
              f.write( space*(17-len(line_copy[index])) )
            else:
              f.write(line)
              f.write(space*(26-len(line)))
            reset +=1
            continue
          if reset == 6:
            f.write(line)
            f.write(space*(12-len(line_copy[index])))
            reset +=1
            continue
          if reset == 7:
            f.write(line)
            reset +=1
            continue
      f.close()

  def File_generation(self,substitute):
    """Generation of the dynpro.in and all the necessities"""
    #############################
    # Conquest coord generation #
    #############################
    # If the file is already created the program will not generate it again.
    Conquest_coord_is_present = False
    if (os.path.isfile(self.working_directory+"Conquest_coord") ):
      Conquest_coord_is_present = True
    if ( not os.path.isfile("Conquest_coord") and not Conquest_coord_is_present):
      # Will look at all the .in files that can't be a dynpro.in file.
      coord_list = []
      files_list = glob.glob("*.in")
      files_list_dir = glob.glob(self.working_directory+"*.in")
      for i in files_list_dir:
          files_list.append(i)
      for file in files_list:
          select_file = open(file,mode="r")
          line = select_file.readline()   # We read the first line only
          if (line != "&inputdynpro\n"):  # We only take the files that can't be in a dynpro format.
              coord_list.append(file)
          select_file.close()

      # Selection of the file to convert.
      if (len(coord_list) == 0):
          print("No possible .in files translatable found.")
      if (len(coord_list) == 1):
          Choosen_coord = coord_list[0]
      if (len(coord_list) > 1):
          More100 = True
          if (len(coord_list) >= 100): # Selection if we have a lot of files to show.
              print("More than 100 files are detected, do you wish to select one?")
              while True:
                  Shall_we = str(input("y or n"))
                  if (Shall_we == "y" or Shall_we == "Y" or Shall_we == "Yes" or Shall_we == "YES" or Shall_we == "1"):
                      More100 = True
                      break
                  if (Shall_we == "n" or Shall_we == "N" or Shall_we == "No" or Shall_we == "NO" or Shall_we == "0"):
                      More100 = False
                      print("The following file is then selected:")
                      print(coord_list[0])
                      Choosen_coord = coord_list[0]
                      break
          if More100 == True:
              print("Multiple possible .in files can be translated, please chose one.")
              print("--------------")
              for i in range(len(coord_list)):
                  print(i," ",coord_list[i])
              print("--------------")
              try:
                  Answer = input("Please, enter only a integer. ")
                  Answer = int(Answer)
              except:
                  while (type(Answer) is not int):
                      try:
                          Answer = input("Please, enter only a integer. ")
                          Answer = int(Answer)
                      except:
                          continue
              Choosen_coord = coord_list[Answer]

      # Generating the Conquest_coord from the file choosen.
      # If the fortran code is edited and doesn't need the "0 0.0" at the begining of the the file
      # the program only have to copy the .in file and rename it into Conquest_coord
      cq_coord = open("Conquest_coord","w")
      cq_coord.write("0 0.0\n")

      input_coord = open(Choosen_coord,"r")
      for lines in input_coord:
          cq_coord.write(lines)

      cq_coord.close()
      input_coord.close()
          
      if ( os.path.isdir(self.working_directory) and len(coord_list) != 0):
          os.system("cp ./Conquest_coord "+self.working_directory)
    elif ( os.path.isfile("Conquest_coord") and not Conquest_coord_is_present):
      os.system("cp ./Conquest_coord "+self.working_directory)
  
    # ct_atom selection, it will take the first atom with the same atomic nature.
    if self.ct_atom <= 0:
      f = open(self.working_directory+"Conquest_coord","r")
      nature_of_atom = np.array([])
      for line in f:
          words   = line.replace(",",".").split()
          if len(words) < 4:
            continue
          if len(words) == 0:
              continue
          if (words[0] == "#"):
              continue
          nature_of_atom       = np.append(nature_of_atom,float(words[3]))
      f.close()

      for ispec in range(self.nspec):
        if (self.species[ispec+1] == self.Atomic_name):
          atom_ID = ispec+1
      index = 0
      for nature in nature_of_atom:
        index += 1
        if str(int(atom_ID)) == str(int(nature)):
          self.ct_atom = index

          print("The automatic central atom as been selected with the ID: "+str(self.ct_atom))
          break
    
    if substitute: #The substitution is only the last atom.
      self.Atomic_number = self.sub_Atomic_number
      self.Atomic_name   = self.substi_Atom_name


    # Folder organisation
    if (not os.path.isfile("Conquest_input") or not os.path.isfile("Conquest_coord")):
        if (not os.path.isfile(self.working_directory+"Conquest_input") or not os.path.isfile(self.working_directory+"Conquest_coord")):
            print("You must have a 'Conquest_input' and a 'Conquest_coord' to effectuate the calculation.")
            print("It should be present in the terminal position or in the directory: "+self.working_directory)
            print("The program will not run due to this error.")
            sys.exit()
    if ( not os.path.isdir(self.working_directory) ):
        print("The files must be in the directory the bash is in.")
        os.makedirs(self.working_directory)
        print("Create the working directory.")
        os.system("cp ./Conquest_input ./Conquest_coord "+self.working_directory)
    
    if ( not os.path.isfile(self.working_directory+"Conquest_input") ):
      os.system("cp ./Conquest_input "+self.working_directory)
    if ( not os.path.isdir(self.working_directory_scratch) ):
        print("The md files must be in the directory the bash is in.")
        os.makedirs(self.working_directory_scratch)
        os.system("cp ./"+self.md_frame_name+" ./"+self.md_stats_name+" "+self.working_directory_scratch)
        print("Create working directory for the MD files.")
    if (not os.path.isfile(self.working_directory_scratch+self.md_frame_name) or not os.path.isfile(self.working_directory_scratch+self.md_stats_name)):
      os.system("cp ./"+self.md_frame_name+" ./"+self.md_stats_name+" "+self.working_directory_scratch)
    # dynpro.in generation
    dynpro_in = open(self.working_directory+self.Input_dynpro_name,"w")
    dynpro_in.write("&inputdynpro\n")
    dynpro_in.write("       prefix    = '"+str(self.Prefix)+"',  		! Prefix of MD data files.\n")
    dynpro_in.write("       outdir    = '."+str(self.working_directory_scratch)+"',     ! Output directory containing MD data files.\n")
    # Force data and it's dependencies will stay commented, but can be added later on.
    dynpro_in.write("       !prefix_cqdata = 'CoordForce'	! Output for the force coordinate file.\n")
    dynpro_in.write("       prefix_cqdata  = '"+str(self.md_frame_name)+"'	! The name of the frame file.\n")
    dynpro_in.write("       lcq_coordforce = F		! Do we use the force coordinate file?\n")
    dynpro_in.write("       lcq_mdframes   = T		! Do we use the frame file?\n")
    dynpro_in.write("       time_step =  "+str(self.time_step)+"d0,   	! MD time step in atomic unit\n")
    if self.center:
      dynpro_in.write("       lcenter   =  T			! Do we center the atom?\n")
      dynpro_in.write("       ct_atom   =  "+str(self.ct_atom)+"			! ID of the centered atom\n")
    else:
      dynpro_in.write("       lcenter   =  F			! Do we center the atom?\n")
      dynpro_in.write("       ct_atom   =  0			! ID of the centered atom\n")
    dynpro_in.write("       lgrace    =  T\n")
    dynpro_in.write("       lsub      =  "+str(self.subtitution)+"    ! Do we do a substitution?\n")
    dynpro_in.write("       lau_unit  =  T        ! We are in atomic units\n")
    dynpro_in.write("       lsample   =  "+str(self.Sample_data)+"			! Do we sample our data?\n")
    dynpro_in.write("       lsphere   =  F\n")
    dynpro_in.write("       lrdf      =  "+str(self.RDF_data)+"			! Do we do the RDF analysis?\n")
    dynpro_in.write("       lvacf     =  "+str(self.VACF_data)+"			! Do we do the VACF analysis (FFT)?\n")
    dynpro_in.write("       lpsd      =  "+str(self.VACF_data)+"			! Another method for the fourier transform.\n")
    dynpro_in.write("       ldebug    =  F      ! Debugging option.\n")
    dynpro_in.write("/\n")
    dynpro_in.write("&inputsample\n")
    if (self.Starting_point < 0 or not self.Sample):
        dynpro_in.write("       !sp_start = 0 		! When does the sample start to be taken?	(if not precised, i.e. commented: 0)\n")
    else:
        dynpro_in.write("       sp_start = "+str(int(self.Starting_point))+" 		! When does the sample start to be taken?	(if not precised, i.e. commented: 0)\n")

    if (self.Ending_point < self.Starting_point or not self.Sample):
        dynpro_in.write("       !sp_end   = 0		! When does the sample finish his selection? 	(if not precised, i.e. commented: the end of the document)\n")
    else:
        dynpro_in.write("       sp_end   = "+str(int(self.Ending_point))+"		! When does the sample finish his selection? 	(if not precised, i.e. commented: the end of the document)\n")

    if ( not self.Sample ):
        dynpro_in.write("       !sp_n = 0 		! Number of sample to be taken\n")
    else:
        dynpro_in.write("       sp_n = "+str(int(self.sp_n))+" 		! Number of sample to be taken\n")

    dynpro_in.write("       iprint    = "+str(1)+"			! Number of steps skiped between each frames.\n")
    dynpro_in.write("/\n")
    dynpro_in.write("&inputrdf\n")
    
    dynpro_in.write("       rd_atom =  "+str(self.Atomic_number)+",			! Atomic number of the atom ("+str(self.Atomic_name)+")\n")
    dynpro_in.write("       rmin    =  "+str(self.r_min)+",			! Minimal radius (cannot be 0)\n")
    dynpro_in.write("       nr      =  "+str(self.nr)+",			! Number of steps to form the histogram\n")
    dynpro_in.write("/\n")
    # Inputpsd
    
    dynpro_in.write("&inputpsd\n")
    dynpro_in.write("       lspe       = "+str(self.VACF_data)+" 			! Species \n")
    dynpro_in.write("       m_cof      =  			! Low sp_n ex. sp_n/2\n")
    dynpro_in.write("       n_seg      =  1\n")
    dynpro_in.write("       loverlap   =  F\n")
    dynpro_in.write("       win        =  0\n")
    if self.temperature < 0:
      dynpro_in.write("       temp       = 298.15d0\n")
      print("An unpexpected error as been identified in the temperature.")
    else:
      dynpro_in.write("       temp       = "+str(self.temperature)+"d0\n")
    dynpro_in.write("       qcfactor   =  0\n")
    dynpro_in.write("/\n")
    # Inputsphere (non-automatised)
    dynpro_in.write("&inputsphere\n")
    dynpro_in.write("       sph_atom   = 197\n")
    dynpro_in.write("       sph_rcut   = 1.2d0\n")
    dynpro_in.write("       sph_n      = 0\n")
    dynpro_in.write("       sph_atn(1) = 8\n")
    dynpro_in.write("       sph_atn(2) = 1\n")
    dynpro_in.write("       sph_atn(4) = 17\n")
    dynpro_in.write("       sph_atn(3) = 78\n")
    dynpro_in.write("       sph_atn(8) = 6\n")
    dynpro_in.write("       sph_atn(5) = 7\n")
    dynpro_in.write("       sph_atn(6) = 2\n")
    dynpro_in.write("       sph_atn(7) = 16\n")
    dynpro_in.write("       langle_spef = T\n")
    dynpro_in.write("       ldist_spef  = T\n")
    dynpro_in.write("/\n")

    dynpro_in.close()

    os.system("mv "+self.working_directory+self.Input_dynpro_name+" "+self.working_directory)

  def dynpro_launch(self):
    """Function to launch dynpro"""
    os.chdir(self.working_directory)  # Go inside the testing file.
    python_loc = sys.argv[0]          # Identify where the python is located
    new_loc, fname = os.path.split(python_loc)  # Will find the dynpro.x executable
    new_loc = new_loc + "/cq_dynpro1.0/src/dynpro.x"
    print("Launching dynpro with the following center atom: "+self.Atomic_name)
    print(">>==========================================================")
    os.system(new_loc+" < "+self.Input_dynpro_name) # We run the program with one possible center.
    os.chdir("../")                   # We go back in the original folder.

    os.system("cp "+self.working_directory+"out.* "+self.working_directory+"UpdatedAtoms* "+self.working_directory_output)
    # Removing the temporary files generated by dynpro
    if (os.path.isfile("./"+self.working_directory+"tmp.txt") and not self.debug ):
      os.system("rm ./"+self.working_directory+"tmp.txt")
    if (os.path.isfile("./"+self.working_directory+"tmp2.txt") and not self.debug ):
      os.system("rm ./"+self.working_directory+"tmp2.txt")

  def custom_name(self):
    # To have a possible customizable name for the marked atom
    print("What is the name of your atom that will be printed on the graphic?")
    print("Write # if you want to take the name of your substitution.")
    Custom_name = str(input("-> "))
    if Custom_name == "#" :
      Custom_name = str(self.substi_Atom_name)
    self.Real_atom_name = np.append(self.Real_atom_name,Custom_name)

  def plot_dynpro(self,graph_show_VACF,graph_show_RDF,central_atom_selected,species,advanced):
    """Will plot the graphic of each datas for the VACF and RDF"""
    #graph_show_VACF => Do we show the VACF graph?
    #graph_show_RDF => Do we show the RDF graph?
    #central_atom_selected => The central atom wished to be printed str().
    # The function need that all the centered atoms are printed.
    #
    # For the graphic display
    VACF_color = "109727"
    grid_color = "d7d7d7"
    line_width_choice = 1.3
    ###################
    # Advanced option #
    ###################
    if advanced:
      print_results = True # To return the maximum coordinance.
      Calculate_max = True # To return the maximums of the RDF plot.
    else:
      print_results = False # To return the maximum coordinance.
      Calculate_max = False # To return the maximums of the RDF plot.

    # The outputs will be in a .txt file.
    #
    atom_list = np.array([])
    #All_done_RDF  = False
    All_done_VACF = False
    ##############################
    # In case of a substitution. #
    ##############################
    if self.subti:
      for ispec in range(self.nspec):
          atom_list = np.append(atom_list,self.species[ispec+1])
      for isubsti in range(len(self.substi_Atom_name)-1):
        atom_list = np.append(atom_list,self.substi_Atom_name)
    else:
      for ispec in range(self.nspec):
          atom_list = np.append(atom_list,self.species[ispec+1])
    
    for atom_c in atom_list:
      if atom_c != central_atom_selected and central_atom_selected != "all":
        continue
      os.chdir(self.working_directory_output) # We go to the good directory.
      if self.VACF:
          central = atom_c
          method = "psd"
          path = "./"
          workplace = path+"out."+method+"."+central+".dat"
          if ( not os.path.isdir("./Plots/VACF/") ):
              os.makedirs("./Plots/VACF/")

          filename_vacf = "./Plots/VACF/VACF_"+central
          print("Treatment of the data in: ",workplace)
          f = open(workplace,"r")
          f.readline()

          Wave       = np.array([]) #   WaveNumber [cm-1]
          Intens     = np.array([]) #   Intensity [arb. unit]

          for line in f:
              words   = line.replace(",",".").replace("NaN","-1").split()
              # We change the "," into "." to avoid localisation problem, and "NaN" to "-1" to be able to plot and we split in function of spaces.
              # For now if the value is negative, its NOT a possible result.
              if (words[0] == "#"):
                  continue
              Wave        = np.append(Wave,float(words[0]))
              Intens      = np.append(Intens,float(words[1]))
          f.close()

          if not All_done_VACF and central_atom_selected == "all":
              workplace_all = path+"out."+method+"_spc"+".dat"
              filename_vacf_all = "./Plots/VACF/"+method+"_all"

              print("Treatment of the data in: ",workplace_all)
              ##################################
              # Extract the datas from a file  #
              ##################################
              f_ = open(workplace_all,"r")
              f_.readline()

              Wave_all        = np.array([]) # WaveNumber [cm-1]
              Intens_all      = np.array([]) #  Intensity [arb. unit]

              for line in f_:
                  words   = line.replace(",",".").replace("NaN","-1").split()
                  # We change the "," into "." to avoid localisation problem, and "NaN" to "-1" to be able to plot and we split in function of spaces.
                  # For now if the value is negative, its NOT a possible result.
                  if (words[0] == "#"):
                      continue
                  Wave_all       = np.append(Wave_all,float(words[0]))
                  Intens_all    = np.append(Intens_all,float(words[1]))
          if self.subti: # Adapt the data name for the visualisation of the graph in case of a substitution.
            for isubsti in range(len(self.substi_Atom_name)):
              if central == self.substi_Atom_name:
                central = self.Real_atom_name[isubsti]
          for ispec in range(self.nspec):
            if self.species[ispec+1] == central and self.species[ispec+1] != species[ispec+1]:
              central = species[ispec+1]
              break

          figVACF,axl = plt.subplots()

          axl.minorticks_on()
          axl.grid( which='major', axis='x', color='#'+grid_color, linestyle='-')
          axl.grid( which='minor', axis='x', color='#'+grid_color, linestyle='--')
          axl.grid( which='major', axis='y', color='#'+grid_color, linestyle='-')

          axl.set_ylabel("Intensity [arb. unit]", color="#"+VACF_color)
          axl.set_xlabel("WaveNumber [cm-1]")

          axl.plot(Wave, Intens, color="#"+VACF_color,linestyle="solid", label="VACF", linewidth=line_width_choice)
          axl.set_ylim(bottom=0)
          axl.set_xlim(left=Wave[0],right=Wave[-1])

          plt.title('VACF with PSD method of '+central)

          if graph_show_VACF:
              plt.show()
          figVACF.savefig(filename_vacf, bbox_inches='tight')
          plt.close()

          if not All_done_VACF and central_atom_selected == "all":
              figVACF_,axl_ = plt.subplots()
              All_done_VACF = True
              axl_.minorticks_on()
              axl_.grid( which='major', axis='x', color='#'+grid_color, linestyle='-')
              axl_.grid( which='minor', axis='x', color='#'+grid_color, linestyle='--')
              axl_.grid( which='major', axis='y', color='#'+grid_color, linestyle='-')


              axl_.set_ylabel("Intensity [arb. unit]", color="#"+VACF_color)
              axl_.set_xlabel("WaveNumber [cm-1]")

              axl_.plot(Wave_all, Intens_all, color="#"+VACF_color,linestyle="solid", label="VACF", linewidth=line_width_choice)
              axl_.set_ylim(bottom=0)
              axl_.set_xlim(left=Wave_all[0],right=Wave_all[-1])

              plt.title('VACF with PSD method of all the atoms')
              figVACF_.savefig(filename_vacf_all, bbox_inches='tight')
              if graph_show_VACF:
                  plt.show()
              plt.close()
      if self.RDF:
          for atom_t in atom_list:
              # The terminal needs to be in the folder directory to work.
              central = atom_c
              targuet = atom_t
              bond = central+"_"+targuet
              method = "rdf"
              path = "./"
              workplace = path+"out."+method+"."+bond+".dat"
              ##################################
              # Extract the datas from a file  #
              ##################################
              print("Treatment of the data in: ",workplace)
              try:
                f = open(workplace,"r")
              except:
                print("Error 404, file not found.")
                continue
              f.readline()

              r       = np.array([]) # In Ang
              g_ab    = np.array([]) # g_ab (r) Normalized
              g_int   = np.array([]) # g_int (r) Normalized

              for line in f:
                  words   = line.replace(",",".").replace("NaN","-1").split()
                  # We change the "," into "." to avoid localisation problem, and "NaN" to "-1" to be able to plot and we split in function of spaces.
                  # For now if the value is negative, its NOT a possible result.
                  if (words[0] == "#"):
                      continue
                  r       = np.append(r,float(words[0]))
                  g_ab    = np.append(g_ab,float(words[1]))
                  g_int   = np.append(g_int,float(words[2]))

              f.close()

              if self.subti: # Adapt the data name for the visualisation of the graph in case of a substitution.
                for isubsti in range(len(self.substi_Atom_name)):
                  if central == self.substi_Atom_name:
                    central = self.Real_atom_name[isubsti]
                    bond = central+"_"+targuet
                  if targuet == self.substi_Atom_name:
                    targuet = self.Real_atom_name[isubsti]
                    bond = central+"_"+targuet
              
              for ispec in range(self.nspec):
                if self.species[ispec+1] == central and self.species[ispec+1] != species[ispec+1]:
                  central = species[ispec+1]
                  bond = central+"_"+targuet
                if self.species[ispec+1] == targuet and self.species[ispec+1] != species[ispec+1]:
                  targuet = species[ispec+1]
                  bond = central+"_"+targuet
              ################################################
              # Calculation of the maximums of the g(r) pics #
              ################################################
              if Calculate_max:
                g_last = 0.0
                g_last_min = 0.0
                g = 0.0
                found = False
                rank_maxim = 0
                Verify = False
                if ( not os.path.isdir("./data_table_dynpro/") ):
                    os.makedirs("./data_table_dynpro/")
                #print("******************")
                #print("Maximums detected:")
                if (self.r_min >= r[0] and Calculate_max):
                    Verify = True
                for i in r:
                    rank_maxim += 1
                    if self.r_min <= float(i) and Verify:
                        # What is the maximum?
                        g = g_ab[rank_maxim-1]
                        Coord = g_int[rank_maxim-1]
                        if g != 0 and g <= g_last and not found:
                            Maximum_g = g
                            #print(i,"A ->",Maximum_g," =>", Coord)
                            if print_results:
                                if ( not os.path.isfile("./data_table_dynpro/log_max_length"+central+"-"+targuet+".txt")):
                                    deposit_file = open("./data_table_dynpro/log_max_length"+central+"-"+targuet+".txt","w")
                                    deposit_file.write("        r       max_g       Coord\n")
                                else:
                                    deposit_file = open("./data_table_dynpro/log_max_length"+central+"-"+targuet+".txt","a")
                                deposit_file.write("    "+str(i)+"   "+str(Maximum_g)+"       "+str(Coord)+"\n")
                                deposit_file.close()
                                found = True
                        else:
                            g_last = g
                        if g != 0 and g >= g_last and found:
                            Minimum_g = g
                            found = False
                        else:
                            g_last_min = g
                #print("******************")
              ###############################
              # Give the final coordinance. #
              ###############################
              if print_results:
                if ( not os.path.isdir("./data_table_dynpro/") ):
                  os.makedirs("./data_table_dynpro/")
                if ( not os.path.isfile("./data_table_dynpro/log_coord"+targuet+".txt")):
                    deposit_file = open("./data_table_dynpro/log_coord"+targuet+".txt","w")
                    deposit_file.write("        "+"X"+"-"+targuet+"       Coord\n")
                    deposit_file.write("        "+central+"-"+targuet+"     "+str(g_int[-1])+"\n")
                else:
                    deposit_file = open("./data_table_dynpro/log_coord"+targuet+".txt","a")
                    deposit_file.write("        "+central+"-"+targuet+"     "+str(g_int[-1])+"\n")
                    deposit_file.close()
              ####################
              # Graphical output #
              ####################
              figRDF,axl = plt.subplots()
              if ( not os.path.isdir("./Plots/RDF/") ):
                  os.makedirs("./Plots/RDF/")
              filename = "./Plots/RDF/RDF_"+bond

              axl.minorticks_on()
              axl.grid(which='major', axis='x', color='#'+grid_color, linestyle='-')
              axl.grid(which='minor', axis='x', color='#'+grid_color, linestyle='--')
              axl.grid(which='major', axis='y', color='#'+grid_color, linestyle='-')

              axr = axl.twinx()

              axl.set_ylabel("g(r)", color="#314ad2")
              axr.set_ylabel("Coordination", color='#d23131')
              axl.set_xlabel("r (A)")


              axl.plot(r, g_ab, color="#314ad2",linestyle="solid", label="g(r)", linewidth=line_width_choice)
              axr.plot(r, g_int, color='#d23131',linestyle="solid", label="Coord", linewidth=line_width_choice)
              plt.xlim((0,r[-1]))
              axl.set_ylim(bottom=0)
              axr.set_ylim(bottom=axl.get_ylim()[0], top=axr.get_ylim()[1])

              plt.title('g(r) and coordinance of '+central+'-'+targuet)

              figRDF.savefig(filename, bbox_inches='tight')
              if graph_show_RDF:
                  plt.show()
              plt.close()
      os.chdir("../")

class MSER:
  """Marginal Standard Error Rule heuristic 
  --- K P White, Simulation 69, 323 (1997)"""

  def __init__(self, nframes, varname, var_traj):
    self.n_j = nframes-1
    self.propname = varname
    self.traj = var_traj
    self.mser = sp.zeros(self.n_j, dtype='float')
    # stop before the end otherwise the MSER becomes very noisy
    self.mser_cut = 200

  def get_point(self, d_j):
    prefac = 1.0/(self.n_j-d_j)**2
    ybar_ij = sp.mean(self.traj[d_j:])
    variance = 0.0
    for i in range(d_j+1,self.n_j):
      variance += (self.traj[i] - ybar_ij)**2
    return prefac*variance

  def get_mser(self):
    for i in range(self.n_j):
      self.mser[i] = self.get_point(i)

  def mser_min(self):
    return sp.argmin(self.mser[:-self.mser_cut])

  def plot_mser(self, steps):
    plt.figure("{} MSER".format(self.propname))
    plt.xlabel("step")
    plt.ylabel("MSER ({})".format(self.propname))
    plt.plot(steps[:-200], self.mser[:-200], 'k-')
    mser_min = self.mser_min()
    lab = "Minimum at step {}".format(mser_min)
    plt.axvline(x=mser_min, label=lab)
    plt.legend(loc="upper left")
    plt.savefig("mser.pdf", bbox_inches='tight')

  def dump_mser(self, steps):
    mser_fmt = "{0:>8d}{1:>16.6f}\n"
    filename = "mser.dat"
    with open(filename, 'w') as outfile:
      for i in range(self.n_j):
        outfile.write(mser_fmt.format(steps[i], self.mser[i]))

class MSD:
  """Mean Squared Deviation"""

  def __init__(self, nat, dt, init_frame,ID_select,msdskip,jump):
    self.nframes = 0
    self.nat = nat
    self.dt = dt
    self.init_r = init_frame.r
    self.r_prev = np.copy(self.init_r)
    self.msd = []
    self.steps = []
    self.init_cell = np.zeros(3, dtype='float')
    self.r_diff = np.zeros((self.nat,3), dtype='float')
    ##For the Diffusion coefficient##
    self.msdskip  = msdskip   # Input from the user, time in (fs)
    self.jump   = jump        # How much data is printed by conquest each time
    self.min_index = int(round(self.msdskip/self.jump))
    for i in range(3):
      self.init_cell[i] = init_frame.lat[i,i]

    ##For multiple species.##
    self.ID_select = ID_select
    self.msd_spec = [[] for _ in range(self.nat)] # Array to store the msd for one specific species.


  def update_msd(self, step, frame):
    self.steps.append(step)
    self.nframes += 1
    cell = np.zeros(3, dtype='float')
    for i in range(3):
      cell[i] = frame.lat[i,i]
    self.msd.append(0.0)

    for i in range(self.nat):
      diff = diff_mic(frame.r[i,:], self.r_prev[i,:], cell) # To verify?
      self.r_diff[i,:] += diff
      self.msd[-1] += np.sum(self.r_diff[i,:]**2)

      self.msd_spec[i].append((np.sum(self.r_diff[i,:]**2))) # For each individual atoms

    self.r_prev = np.copy(frame.r)

  def norm_msd(self):
    self.msd = np.array(self.msd)/self.nat
    self.msd_spec = np.array(self.msd_spec)
    self.time = np.array(self.steps, dtype='float')*self.dt

  def Diffusion_msd(self):
    x,y = self.time[self.min_index:], self.msd[self.min_index:]
    popt,_ = sp.optimize.curve_fit(linear_regression_diffusion,x,y)
    self.a,self.b = popt
    self.diffusion = (self.a*((bohr2ang*1e-8)**2)*1e+15)/6.0 # We have the diffusion in cm^2/s
    print("Diffusion coefficient: ",self.diffusion)

  def Diffusion_msd_multiple(self):
    #######################################################################
    # The function will be taken in function of the ID of the atom chosen.
    #######################################################################
    x,y = self.time[self.min_index:], self.msd_spec[self.ID_select-1][self.min_index:]
    popt,_ = sp.optimize.curve_fit(linear_regression_diffusion,x,y)
    self.a_multi,self.b_multi = popt
    self.diffusion_spec = (self.a_multi*((bohr2ang*1e-8)**2)*1e+15)/6.0 # We have the diffusion in cm^2/s
    print("Diffusion coefficient: ",self.diffusion_spec)

  def plot_msd(self):
    filename = "msd.pdf"
    plt.figure("MSD")
    plt.xlabel("t (fs)")
    plt.ylabel("MSD ($Bhor^2$)") # the data in input is in Bhor.
    plt.xlim((self.time[0],self.time[-1]))
    plt.plot(self.time, self.msd)
    plt.plot(self.time,self.a*self.time+self.b,linestyle="--", color="#e37362",label="Diffusion coefficient = "+"{0:>10.4e}".format(self.diffusion)+" $cm^2.s^{-1}$")
    plt.plot((0,self.time[-1]), (0, 0), 'k-')
    plt.ylim(ymin=0)
    plt.legend(prop={"size":8})
    plt.savefig(filename, bbox_inches='tight')

  def plot_msd_multiple(self):
    filename2 = "msd_atm"+str(self.ID_select)+".pdf"
    plt.figure("MSD")
    plt.xlabel("t (fs)")
    plt.ylabel("MSD ($Bhor^2$)") # the data in input is in Bhor.
    plt.xlim((self.time[0],self.time[-1]))
    plt.plot(self.time, self.msd_spec[self.ID_select-1][:])
    plt.plot(self.time,self.a_multi*self.time+self.b_multi,linestyle="--", color="#62b0e3",label="Diffusion coefficient of the selected atom = "+"{0:>10.4e}".format(self.diffusion_spec)+" $cm^2.s^{-1}$")
    plt.plot((0,self.time[-1]), (0, 0), 'k-')
    plt.ylim(ymin=0)
    plt.legend(prop={"size":8})
    plt.savefig(filename2, bbox_inches='tight')

  def dump_msd(self):
    msd_fmt = "{0:>12.4f}{1:>16.6f}\n"
    filename = "msd.dat"
    with open(filename, 'w') as outfile:
      for i in range(self.nframes):
        outfile.write(msd_fmt.format(self.time[i], self.msd[i]))
