#!/usr/local/bin/python3

import argparse
import sys
import re
import os.path
import scipy as sp
import matplotlib.pyplot as plt
from md_tools import Pairdist
from frame import Frame
from pdb import set_trace

# Parsing functions

specblock_re = re.compile(r'%block ChemicalSpeciesLabel\n(.*?)\n%endblock', re.M | re.S | re.I)
qe_data_re = re.compile(r'Begin final coordinates(.*?)End final coordinates', re.M | re.S)
qe_mass_re = re.compile(r'mass(.*?)\n\n', re.M | re.S)

bohr2ang = 0.529177249
ang2cm = 1.0e-8
avogadro = 6.022140857e23

def parse_qe_output(qe_output_file):
  frac = False
  angstrom = False
  with open(qe_output_file, 'r') as infile:
    text = infile.read()

  cell_data = {}
  m = re.search(qe_data_re, text)
  data = m.group(1).strip().splitlines()
  for i, line in enumerate(data):
    if "CELL_PARAMETERS" in line:
      i_cell = i
    if "ATOMIC_POSITIONS" in line:
      i_atom = i
      if 'crystal' in line:
        frac = True
      elif 'angstrom' in line:
        angstrom = True
  m = re.search(qe_mass_re, text)
  mdata = m.group(1).splitlines()
  mass = {}
  for line in mdata[1:]:
    bits = line.split()
    spec = bits[0]
    m = bits[2]
    mass[spec] = float(m)

  cell = data[i_cell+1:i_cell+4]
  pos = data[i_atom+1:]
  spec_list = []
  spec_ind = {}
  scount = {}
  species = []
  coords = []
  a = [float(bit) for bit in cell[0].split()]
  b = [float(bit) for bit in cell[1].split()]
  c = [float(bit) for bit in cell[2].split()]
  for line in pos:
    s, x, y, z = line.split()
    x = float(x)
    y = float(y)
    z = float(z)
    if not s in spec_list:
      spec_list.append(s)
    si = spec_list.index(s)+1
    if not si in scount.keys():
      scount[si] = 1
    else:
      scount[si] += 1
    coords.append([x,y,z])
    species.append(si)

  natoms = len(species)
  for i, spec in enumerate(spec_list):
    spec_ind[i+1] = spec

  coords = sp.array(coords)
  lat = sp.array([a,b,c])
  if angstrom:
    lat = lat/bohr2ang
  if frac:
    for i in range(natoms):
      coords[i,:] = sp.matmul(lat, coords[i,:])

  data = {}
  data['latvec'] = lat
  data['coords'] = coords
  data['species'] = sp.array(species)
  data['natoms'] = natoms
  data['nspecies'] = len(spec_list)
  data['spec_list'] = spec_ind
  data['species_count'] = scount
  data['mass'] = mass
  return data

def strip_comments(line, separator):
  for s in separator:
    i = line.find(s)
    if i >= 0:
      line = line[:i]
  return line.strip()

def parse_cq_input(cq_input_file):
  cq_params = {}
  with open(cq_input_file, 'r') as cqip:
    for line in cqip:
      stripped = strip_comments(line, "#%!")
      if stripped:
        bits = stripped.split()
        if len(bits[1:]) == 1:
          cq_params[bits[0]] = bits[1:][0]
        else:
          cq_params[bits[0]] = bits[1:]

  # get the species labels
  cq_params['species'] = {}
  cq_params['mass'] = {}
  with open(cq_input_file, 'r') as cqip:
    text = cqip.read()
    m = re.search(specblock_re, text)
    specinfo = m.group(1).splitlines()
    for line in specinfo:
      bits = line.split()
      index = int(bits[0])
      mass = float(bits[1])
      spec = bits[2]
      cq_params['species'][int(index)] = spec
      cq_params['mass'][spec] = mass

  return cq_params

def parse_init_config(conf_filename, frac):
  data = {}
  with open(conf_filename, 'r') as infile:
    a = [float(bit) for bit in infile.readline().strip().split()]
    b = [float(bit) for bit in infile.readline().strip().split()]
    c = [float(bit) for bit in infile.readline().strip().split()]
    data['latvec'] = sp.array([a,b,c])
    cell = sp.array([a[0], b[1], c[2]])
    natoms = int(infile.readline().strip())
    data['natoms'] = natoms
    coords = []
    species = []
    for i in range(natoms):
      bits = infile.readline().strip().split()
      x = float(bits[0])
      y = float(bits[1])
      z = float(bits[2])
      spec = int(bits[3])
      pos = sp.array([x,y,z])
      if frac:
        pos *= cell
      coords.append(pos)
      species.append(int(spec))
    data['coords'] = sp.array(coords)
    data['species'] = sp.array(species)
    scount = {}
    for i in range(natoms):
      if data['species'][i] in scount.keys():
        scount[data['species'][i]] += 1
      else:
        scount[data['species'][i]] = 1
    data['species_count'] = scount
    data['nspecies'] = len(scount.keys())
  return data

# Command line arguments
parser = argparse.ArgumentParser(description='Compute bond lengths for a  single \
          configuration', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-i', '--infile', action='store', dest='infile',
                    default='coord_next.dat', help='Conquest format structure file')
parser.add_argument('--nbins', action='store', dest='nbins', default=100,
                    help='Number of histogram bins')
parser.add_argument('-c', '--cutoff', nargs='+', action='store', default=None,
                    dest='cutoff', 
                    help='Bond length cutoff matrix (upper triangular part, in rows')
parser.add_argument('--printall', action='store_true', default=False,
                    dest='printall', help='Printa all bond lengths')
parser.add_argument('--qe', action='store_true', default=False, dest='flag_qe',
                    help='Parse a Quantum Espresso output')
parser.add_argument('--cq', action='store_true', default=False, dest='flag_cq',
                    help='Parse a Conquest file')


opts = parser.parse_args()

cutoff = [float(bit) for bit in opts.cutoff]

if opts.flag_cq:
  cq_input_file = "Conquest_input"
  cq_params = parse_cq_input(cq_input_file)
  if cq_params['IO.FractionalAtomicCoords'][0] == 'T':
    frac = True
  else:
    frac = False
  cell_data = parse_init_config(opts.infile, frac)
  cell_data['spec_list'] = cq_params['species']
  cell_data['mass'] = cq_params['mass']
elif opts.flag_qe:
  cell_data = parse_qe_output(opts.infile)
else:
  sys.exit("No calculation type flag selected! (--qe, --cq)")
nspec = cell_data['nspecies']

bondcut = sp.zeros((nspec,nspec))
k = 0
for i in range(nspec):
  for j in range(i, nspec):
    bondcut[i,j] = cutoff[k]
    bondcut[j,i] = cutoff[k]
    k+=1

f = Frame(cell_data['natoms'], 1)
f.r = cell_data['coords']
f.lat = cell_data['latvec']
f.species = cell_data['species']
small = 0.1
rdfwidth = 0.1
rdfcut = min(f.lat[0,0], f.lat[1,1], f.lat[2,2])/2.0 + small

pairdist = Pairdist(cell_data['natoms'], cell_data['nspecies'], rdfcut,
                    rdfwidth, cell_data['spec_list'], cell_data['species_count'])
pairdist.update_rdf(f)
pairdist.get_bondlength(bondcut, f, opts.printall)

v = sp.dot(f.lat[:,0], sp.cross(f.lat[:,1],f.lat[:,2]))
m = 0.0
for s in cell_data['species_count'].keys():
  spec = cell_data['spec_list'][s]
  m += cell_data['species_count'][s]*cell_data['mass'][spec]
v = v * bohr2ang**3 * ang2cm**3
m = m / avogadro
density = m/v
print(f"Density: {density:>10.4f} g/cm^3")
