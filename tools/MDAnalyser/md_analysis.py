#!/usr/local/bin/python3

import argparse
import sys
import re
import os.path
import scipy as sp
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.offsetbox as osb
from md_frame import Frame
from md_tools import MSER, MSD, Dynpro_gestion, autocorr, Atom_selector, logic_selector, Unknown_atom
from pdb import set_trace

# Constants.
ha2ev = 27.211399
ha2k = 3.15737513e5

# Regular expressions
frame_re = re.compile('frame')
endframe_re = re.compile('end frame')
specblock_re = re.compile(r'%block ChemicalSpeciesLabel\n(.*?)\n%endblock',
                          re.M | re.S | re.I)

# The input.log and the Conquest input file is needed for the analysis.
cq_input_file = 'input.log'
cq_input_species = 'Conquest_input'

####################
# Parsing functions
####################
def strip_comments(line, separator): # Function to remove all the comments.
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
  with open(cq_input_species, 'r') as cqip:
    m = re.search(specblock_re, cqip.read())
    specinfo = m.group(1).splitlines()
    for line in specinfo:
      try:
        index, mass, spec  = line.split()
      except ValueError:
        index, mass, spec, ionfile = line.split()
      cq_params['species'][int(index)] = spec

  return cq_params

def parse_cq_input_edit(cq_input_file):
  cq_params_ = {}
  with open(cq_input_file, 'r') as cqip:
    for line in cqip:
      stripped = strip_comments(line, "#%!")
      if stripped:
        bits = stripped.split()
        if len(bits[1:]) == 1:
          cq_params_[bits[0]] = bits[1:][0]
        else:
          cq_params_[bits[0]] = bits[1:]

  # get the species labels
  cq_params_['species'] = {}
  with open("./datas_dynpro/"+cq_input_species, 'r') as cqip:
    m = re.search(specblock_re, cqip.read())
    print(m)
    specinfo = m.group(1).splitlines()
    for line in specinfo:
      try:
        index, mass, spec  = line.split()
      except ValueError:
        index, mass, spec, ionfile = line.split()
      cq_params_['species'][int(index)] = spec

  return cq_params_


def parse_init_config(conf_filename):
  data = {}
  with open(conf_filename, 'r') as infile:
    a = [float(bit) for bit in infile.readline().strip().split()]
    b = [float(bit) for bit in infile.readline().strip().split()]
    c = [float(bit) for bit in infile.readline().strip().split()]
    data['latvec'] = np.array([a,b,c])
    natoms = int(infile.readline().strip())
    data['natoms'] = natoms
    coords = []
    species = []
    for i in range(natoms):
      x, y, z, spec, cx, cy, cz = infile.readline().strip().split()
      coords.append([float(x), float(y), float(z)])
      species.append(int(spec))
    data['coords']  = np.array(coords)
    data['species'] = np.array(species)
    scount = {}
    for i in range(natoms):
      if data['species'][i] in scount.keys():
        scount[data['species'][i]] += 1
      else:
        scount[data['species'][i]] = 1
    data['species_count'] = scount
    data['nspecies'] = len(scount.keys())
    data['volume'] = data['latvec'][0,0]*data['latvec'][1,1]*data['latvec'][2,2]
  return data

def read_stats(stats_file, nstop):
  nstep = 0
  data = {}
  header = True
  with open(stats_file, 'r') as statfile:
    for line in statfile:
      if nstop != -1:
        if nstep > nstop:
          break
      if header:
        col_id = line.strip().split()
        for col in col_id:
          data[col] = []
        header = False
      else:
        bits = line.strip().split()
        for i, bit in enumerate(bits):
          if i==0:
            info = int(bit)
          else:
            info = float(bit)
          data[col_id[i]].append(info)
      nstep += 1
    for key in data:
      data[key] = np.array(data[key])
  return nstep, data

###########################
# Command line arguments
###########################
parser = argparse.ArgumentParser(description='Analyse a Conquest MD \
        trajectory', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-c', '--compare', action='store_true', default=False,
                    dest='compare', help='Compare statistics of trajectories \
                    in directories specified by -d')
parser.add_argument('-d', nargs='+', default=['.',], dest='dirs',
                    action='store', help='Directories to compare')
parser.add_argument('--description', nargs='+', default='', dest='desc',
                    action='store', help='Description of graph for legend \
                    (only if using --compare)')
                    # Command to store the frame file.
parser.add_argument('-f', action='store', dest='framesfile',
                    default='Frames', help='MD frames file')
                    # Command to store the stats file.
parser.add_argument('-s', action='store', dest='statfile',
                    default='Stats', help='MD statistics file')
                    # Command to skip a number of step to not take into account the equilibration steps.
parser.add_argument('--skip', action='store', dest='nskip', default=0,
                    type=int, help='Number of equilibration steps to skip')
                    # Command to analyse every nth steps from the frame file.
parser.add_argument('--stride', action='store', dest='stride', default=1,
                    type=int, help='Only analyse every nth step of frames file')
                    #
parser.add_argument('--snap', action='store', dest='snap', default=-1, 
                    type=int, help='Analyse Frame of a single snapshot')
                    #
parser.add_argument('--stop', action='store', dest='nstop', default=-1, 
                    type=int, help='Number of last frame in analysis')
                    #
parser.add_argument('--equil', action='store', dest='nequil', default=0, 
                    type=int, help='Number of equilibration steps')
                    # Command to plot the statistics.
parser.add_argument('--stats', action='store_true', dest='stats',
                    help='Plot statistics')
                    # Command to compute and plot the VACF.
parser.add_argument('--vacf', action='store_true', dest='vacf', 
                    help='Plot velocity autocorrelation function')
                    # Command to compute and plot the HFACF.
parser.add_argument('--hfacf', action='store_true', dest='hfacf', 
                    help='Plot heat flux autocorrelation function')

parser.add_argument('--acfwindow', action='store', dest='acfwindow', default=0.0, 
                    type=float, help='window for autocorrelation')
                    # Command to compute and plot the mean squared deviation.
parser.add_argument('--msd', action='store_true', dest='msd', 
                    help='Plot mean squared deviation')
                    # Command to skip the begining time to evaluate the diffusion.
parser.add_argument('--msd_skip', action='store', dest='msd_skip', default=0,
                    type=int, help='Time in femto-seconds to skip for the evaluation \
                    of the diffusion')
                    # Command to compute and plot the mean squared deviation of a specific atom.
parser.add_argument('--msd_multi', action='store', dest='msd_multi_ID', default=0, type=int, 
                    help='ID of the atom to plot the mean squared deviation.')
                    # Command to compute and plot the RDF.
parser.add_argument('--rdf', action='store_true', dest='rdf', 
                    help='Plot radial distribution function')
                    # Command to plot the stress.
parser.add_argument('--stress', action='store_true', dest='stress', 
                    help='Plot stress')
                    # Modify the orientation of the stat plot.
parser.add_argument('--landscape', action='store_true', dest='landscape', 
                    help='Generate plot with landscape orientation')
                    # Option to modify the color of the stat plot to take in account colors visible for color blind.
parser.add_argument('--colorblind', action='store_true', dest='colorblind', 
                    help='The plot of the statistics generated will take colors visible \
                    for color blind viewers.')                   
                    #
parser.add_argument('--pub', action='store_true', dest='pub',
                    help='Publication text size')
                    # RDF command: Number of histograms needed.
parser.add_argument('--nbins', action='store', dest='nbins', default=100,
                    help='Number of histogram bins')
                    # RDF command: width of the RDF histograms.
parser.add_argument('--rdfwidth', action='store', dest='rdfwidth',
                    default=0.05, help='RDF histogram bin width (A)')
                    # RDF command: cutoff distance.
parser.add_argument('--rdfcut', action='store', dest='rdfcut', default=10.0,
                    help='Distance cutoff for RDF')
                    # RDF command: minimal distance necessary for the RDF
parser.add_argument('--r_min', action='store', dest='r_min', default=0.1,
                    help='Minimal distance for the RDF')
                    # VACF command: Does the sampling will be done?
parser.add_argument('--vacf_sample', action='store_true', dest='vacf_sample', 
                    help='Sampling of the data for the VACF.')
                    # VACF command: Starting frame of the VACF
parser.add_argument('--vacf_start', action='store', dest='vacf_start', default=0,
                    help='Starting frame of the VACF with dynpro (need --vacf_end)')
                    # VACF command: Last frame calculated of the VACF
parser.add_argument('--vacf_end', action='store', dest='vacf_end', default=0,
                    help='Ending frame of the VACF with dynpro.')
                    # VACF command: Number of sample wished by the user.
parser.add_argument('--vacf_nbr', action='store', dest='vacf_nbr_sample', default=0,
                    help='Number of sample with dynpro, must be a power of 2.')
                    # VACF command: Allow the centering of the atom with dynpro
parser.add_argument('--vacf_center', action='store_true', dest='lcenter', 
                    help='Centering the atom with dynpro')
                    # VACF command: Input of the ID of the central atom we want to take.
parser.add_argument('--central_atom', action='store', dest='ct_atom', default=0,
                    help='ID of the central atom taken (need --vacf_center).')
                    # Substitution option (to be fully implemented)
parser.add_argument('--substi', action='store_true', dest='dynpro_substi', 
                    help='A substitution of atoms need to be done.')
                    # Name of the dynpro.in
parser.add_argument('-dynpro_in', action='store', dest='dynpro_in', default="dynpro.in",
                    help='Name of the .in file taken by dynpro. Must be written with\
                    .in at the end.')
                    # Dynpro command: Only use the graphical display for dynpro output.
parser.add_argument('--dynpro_graph', action='store_true', dest='dynpro_graph_only', 
                    help='Display only the graphical display without launching dynpro')
                    # Prefix for the md data files treated by dynpro
parser.add_argument('-dynpro_md', action='store', dest='dynpro_md', default="DynPro",
                    help='Name of the MD treated files by dynpro.')
                    # Command to allow the output of the maximums of g(r) of the RDF and the coordinance of each atoms.
parser.add_argument('--advanced_rdf', action='store_true', dest='advanced_rdf',
                    help='Print the coordination and the localisation of the g(r) maximums (RDF).') # The coordination printed is the one below the maximum only.
                    #                  
parser.add_argument('--dump', action='store_true', dest='dump', 
                    help='Dump secondary data used to generate plots')
                    #
parser.add_argument('--mser', action='store', dest='mser_var', default=None,
                    type=str, help='Compute MSER for the given property')


opts = parser.parse_args()
if (opts.vacf or opts.msd or opts.stress or opts.rdf):
  read_frames = True
else:
  read_frames = False

if opts.msd_multi_ID > 0 or opts.msd or opts.rdf or opts.vacf:
  cq_params = parse_cq_input(cq_input_file)
  cq_jump = float(cq_params['AtomMove.OutputFreq']) # Take the frequency of output of conquest
  cq_Temperature = float(cq_params['AtomMove.IonTemperature'])

# if opts.nskip > 0:
#   opts.nequil = opts.nskip

if opts.nequil == 0:
  opts.nequil = opts.nskip

if not opts.compare:
  # Parse the input structure and Conquest_input files
  cq_params = parse_cq_input(cq_input_file)
  init_config = parse_init_config(cq_params['IO.Coordinates'])
  natoms = init_config['natoms']
  dt = float(cq_params['AtomMove.Timestep'])
  species = cq_params['species']
  extended_system = False
  if 'MD.Thermostat' in cq_params.keys():
    if cq_params['MD.Thermostat'] == 'nhc':
      extended_system = True
    if cq_params['MD.Thermostat'] == 'ssm':
      extended_system = True
    if cq_params['MD.Thermostat'] == 'svr':
      extended_system = False
  if 'MD.Barostat' in cq_params.keys():
    if cq_params['MD.Barostat'] == 'iso-ssm':
      extended_system = True
    if cq_params['MD.Barostat'] == 'ortho-ssm':
      extended_system = True
    if cq_params['MD.Barostat'] == 'iso-mttk':
      extended_system = True

  # Parse the statistics file
  nsteps, data = read_stats(opts.statfile,opts.nstop)
  avg = {}
  std = {}
  for key in data:
    data[key] = np.array(data[key])
    avg[key]  = np.mean(data[key][opts.nequil:-1])
    std[key]  = np.std(data[key][opts.nequil:-1])
    
  time = [float(s)*dt for s in data['step']]
  data['time'] = np.array(time)
  
  #######################
  # Plot the statistics #
  #######################
  if opts.stats:
    # The color chosen can be distinguished by colorblinds.
    if opts.colorblind == True:
      color_red = "F05039"
      color_blue= "1F449C"
      color_green= "7CA1CC"
      color_pink = "EEBAB4"
      color_cyan = "A8B6CC"
    else:
      color_red = "a63a30"
      color_blue= "183bad"
      color_green= "1c721f"
      color_pink = "9325ac"
      color_cyan = "0da8c2"

    if opts.landscape:
      # Parameters to obtain a good landscape graphic.
      figsizex = 13
      figsizey = 7
      numb_row = 2
      numb_col = 2
      fontsizelegend = 10
      axis_size = fontsizelegend*(7)/10
      left_location = [-0.18,0.5]
      right_location = [1.21,0.5]
      linewidth_stat = 0.5
      legend_correction = 0.06
      fig1, ((ax1, ax2), (ax3, ax4)) = plt.subplots(nrows=numb_row, ncols=numb_col, sharex=True, figsize=(figsizex,figsizey))
    else:
      # Parameters to obtain a good graphic with only one colomn.
      figsizex = 8
      figsizey = 15
      numb_row = 4
      numb_col = 1
      fontsizelegend= 10
      axis_size = fontsizelegend*(7)/10
      left_location = [-0.10,0.5]
      right_location = [1.12,0.5]
      linewidth_stat = 0.5
      legend_correction = 0.03
      fig1, (ax1, ax2, ax3, ax4) = plt.subplots(nrows=numb_row, ncols=numb_col, sharex=True, figsize=(figsizex,figsizey), gridspec_kw={"height_ratios": [2,1,1,1]})

    # Legend system.
    ybox_neut = osb.TextArea("E (Ha)", textprops=dict(color="k",size=fontsizelegend*1.05,rotation=90, ha='left', va='bottom'))
    ybox_neut2 = osb.TextArea("E (Ha)", textprops=dict(color="k",size=fontsizelegend*1.05,rotation=270, ha='right', va='bottom'))
    ybox_pe = osb.TextArea("Potential energy", textprops=dict(color="#"+color_red,size=fontsizelegend,rotation=90, ha='left', va='bottom'))
    
    ybox_ke = osb.TextArea("Kinetic energy", textprops=dict(color="#"+color_blue,size=fontsizelegend,rotation=270, ha='right', va='bottom'))
    Ensemble_box = np.array([])

    Ensemble_box = np.append(Ensemble_box, ybox_ke)

    ax1a = ax1.twinx()
    ax1.plot(data['time'][opts.nskip:], data['pe'][opts.nskip:], color="#"+color_red, linestyle="-", linewidth=linewidth_stat, label='Potential energy')
    ax1a.plot(data['time'][opts.nskip:], data['ke'][opts.nskip:], color="#"+color_blue, linestyle="-", linewidth=linewidth_stat, label='Kinetic energy')

    if cq_params['MD.Ensemble'][2] == 't':
      if cq_params['MD.Thermostat'] == 'nhc':
        ax1a.plot(data['time'][opts.nskip:], data['thermostat'][opts.nskip:], color="#"+color_green, linestyle="-", linewidth=linewidth_stat, label='Thermostat energy')
      if cq_params['MD.Thermostat'] == 'svr':
        ax1a.plot(data['time'][opts.nskip:], data['thermostat'][opts.nskip:],color="#"+color_green, linestyle="-", linewidth=linewidth_stat, label='Thermostat energy')
      # For the legend
      ybox_EnsembleT = osb.TextArea(", Thermostat energy", textprops=dict(color="#"+color_green,size=fontsizelegend,rotation=270, ha='right', va='bottom'))
      Ensemble_box = np.append(Ensemble_box,ybox_EnsembleT)

    if cq_params['MD.Ensemble'][1] == 'p':
      if extended_system:
        ax1a.plot(data['time'][opts.nskip:], data['box'][opts.nskip:], color="#"+color_cyan, linestyle="-", linewidth=linewidth_stat, label='Barostat energy')
        # For the legend
        ybox_EnsembleBaro = osb.TextArea(", Barostat energy", textprops=dict(color="#"+color_cyan,size=fontsizelegend,rotation=270, ha='right', va='bottom'))
        Ensemble_box = np.append(Ensemble_box,ybox_EnsembleBaro)
      ax1a.plot(data['time'][opts.nskip:], data['pV'][opts.nskip:],color="#"+color_pink, linestyle="-", linewidth=linewidth_stat, label='pV')
      # For the legend
      ybox_EnsemblePV = osb.TextArea(", pV", textprops=dict(color="#"+color_pink,size=fontsizelegend,rotation=270, ha='right', va='bottom'))
      Ensemble_box = np.append(Ensemble_box,ybox_EnsemblePV)

    ax2.plot(data['time'][opts.nskip:], data['H\''][opts.nskip:])
    ax2.plot((opts.nskip,data['time'][-1]), (avg['H\''],avg['H\'']), '-',
          label=r'$\langle H\' \rangle$ = {0:>12.4f} $\pm$ {1:<12.4f}'.format(avg['H\''], std['H\'']))

    ax3.plot(data['time'][opts.nskip:], data['T'][opts.nskip:])
    ax3.plot((opts.nskip,data['time'][-1]), (avg['T'],avg['T']), '-',
          label=r'$\langle T \rangle$ = {0:>12.4f} $\pm$ {1:<12.4f}'.format(avg['T'], std['T']))

    ax4.plot(data['time'][opts.nskip:], data['P'][opts.nskip:], color="#"+color_blue, linestyle="-", linewidth=linewidth_stat)
    ax4.plot((opts.nskip,data['time'][-1]), (avg['P'],avg['P']),color="#"+color_cyan, linestyle="--",
          label=r'$\langle P \rangle$ = {0:>12.4f} $\pm$ {1:<12.4f}'.format(avg['P'], std['P']))

    if cq_params['MD.Ensemble'][1] == 'p':
      ax4a = ax4.twinx()
      ax4a.plot(data['time'][opts.nskip:], data['V'][opts.nskip:], color="#"+color_red, linestyle="-", linewidth=linewidth_stat)
      ax4a.plot((opts.nskip,data['time'][-1]), (avg['V'],avg['V']), color="#"+color_pink, linestyle="--",
                label=r'$\langle V \rangle$ = {0:>12.4f} $\pm$ {1:<12.4f}'.format(avg['V'], std['V']))


    # Adapt ifself in function of the ensemble used

    ybox = osb.VPacker(children=[ybox_neut], align="bottom", pad=0,sep=5)
    anchored_ybox = osb.AnchoredOffsetbox(loc=5, child=ybox, pad=0.,frameon=False,bbox_to_anchor=(float(left_location[0]),float(left_location[1])),
                                          bbox_transform=ax1.transAxes,borderpad=0.)
    ax1.add_artist(anchored_ybox)
    ybox = osb.VPacker(children=[ybox_pe], align="bottom", pad=0,sep=5)
    anchored_ybox = osb.AnchoredOffsetbox(loc=5, child=ybox, pad=0.,frameon=False,bbox_to_anchor=(float(left_location[0])+legend_correction,float(left_location[1])),
                                          bbox_transform=ax1.transAxes,borderpad=0.)
    ax1.add_artist(anchored_ybox)

    ybox = osb.VPacker(children=[ybox_neut2], align="bottom", pad=0,sep=5)
    anchored_ybox = osb.AnchoredOffsetbox(loc=5, child=ybox, pad=0.,frameon=False,bbox_to_anchor=(float(right_location[0]),float(right_location[1])),
                                          bbox_transform=ax1.transAxes,borderpad=0.)
    ax1a.add_artist(anchored_ybox)

    ybox = osb.VPacker(children=Ensemble_box, align="bottom", pad=0,sep=5)
    anchored_ybox = osb.AnchoredOffsetbox(loc=5, child=ybox, pad=0.,frameon=False,bbox_to_anchor=(float(right_location[0])-legend_correction,float(right_location[1])),
                                          bbox_transform=ax1.transAxes,borderpad=0.)
    ax1a.add_artist(anchored_ybox)

    ###########
    ##Legend.##
    ###########
    # Set the size of the y ticks of the axis
    for label in (ax1.get_yticklabels()):
      label.set_fontsize(axis_size)
    for label in (ax1a.get_yticklabels()):
      label.set_fontsize(axis_size)
    for label in (ax2.get_yticklabels()):
      label.set_fontsize(axis_size)
    for label in (ax3.get_yticklabels()):
      label.set_fontsize(axis_size)
    for label in (ax4.get_yticklabels()):
      label.set_fontsize(axis_size)
    if cq_params['MD.Ensemble'][1] == 'p':
      for label in (ax4a.get_yticklabels()):
        label.set_fontsize(axis_size)
    # Special conditions in case of a landscape graphic
    if opts.landscape:
      plt.tight_layout(pad=7)
      fig1.subplots_adjust(hspace=0,top=0.97,bottom=0.08,left=0.105,right=0.9)
      ax3.set_xlabel("time (fs)")
    ax2.set_ylabel("H$'$ (Ha)")
    ax3.set_ylabel("T (K)")
    ax4.set_ylabel("P (GPa)", color="#"+color_blue)
    if cq_params['MD.Ensemble'][1] == 'p':
      ax4a.set_ylabel("V ($a_0^3$)", color="#"+color_red)
    ax4.set_xlabel("time (fs)")
    ax2.legend()
    ax3.legend()
    ax4.legend(loc="upper left")
    if cq_params['MD.Ensemble'][1] == 'p':
      ax4a.legend(loc="lower right")
    plt.xlim((opts.nskip,data['time'][-1]))
    if not opts.landscape:
      fig1.subplots_adjust(hspace=0,top=0.97,bottom=0.06,left=0.105,right=0.9)
    fig1.savefig("stats.pdf", bbox_inches='tight')
else:
  # If we're comparing statistics in several directories, use a simplified plot
  fig1, (ax1, ax2, ax3) = plt.subplots(nrows=3, ncols=1, sharex=True, figsize=(7,7))
  ax1a = ax1.twinx()
  for ind, d in enumerate(opts.dirs):
    path = os.path.join(d, cq_input_file)
    cq_params = parse_cq_input(path)
    path = os.path.join(d, cq_params['IO.Coordinates'])
    init_config = parse_init_config(path)
    natoms = init_config['natoms']
    dt = float(cq_params['AtomMove.Timestep'])
    species = cq_params['species']
  
    path = os.path.join(d, opts.statfile)
    nsteps, data = read_stats(path,opts.nstop)
    time = [float(s)*dt for s in data['step']]
    data['time'] = np.array(time)
    print(opts.dirs)
    print(opts.desc)
    
    ax1.plot(data['time'][opts.nskip:], data['H\''][opts.nskip:],
            linewidth=0.5, label=opts.desc[ind])
    y1,y2 = ax1.get_ylim()
    ax1a.set_ylim(y1*ha2k,y2*ha2k)
    ax2.plot(data['time'][opts.nskip:], data['T'][opts.nskip:],
            linewidth=0.5, label=opts.desc[ind])
    ax3.plot(data['time'][opts.nskip:], data['P'][opts.nskip:],
            linewidth=0.5, label=opts.desc[ind])

    ax1.set_ylabel("H$'$ (Ha)")
    ax1a.set_ylabel("H$'$ (K)")
    ax2.set_ylabel("T (K)")
    ax3.set_ylabel("P (GPa)")
    ax3.set_xlabel("time (fs)")
    ax1.legend()
    plt.xlim((opts.nskip,data['time'][-1]))
    fig1.subplots_adjust(hspace=0)
    fig1.savefig("stats.pdf", bbox_inches='tight')

# Plot MSER
if opts.mser_var:
  traj = MSER(nsteps, opts.mser_var, data[opts.mser_var])
  traj.get_mser()
  traj.plot_mser(data['step'])
  if opts.dump:
    traj.dump_mser(data['step'])

# Plot heat flux autocorrelation function
if opts.hfacf:
  window = int(opts.acfwindow // dt)
  G = sp.zeros((3,3,window))
  time = np.array([float(i)*dt for i in range(window)])
  nruns = 0
  for ind, d in enumerate(opts.dirs):
    path = os.path.join(d, heatfluxfile) # Undefined name.

    J = []
    t = []
    nsteps = 0
    with open(path, 'r') as infile:
      for line in infile:
        step, Jx, Jy, Jz = line.split()
        step = int(step)
        Jx = float(Jx)
        Jy = float(Jy)
        Jz = float(Jz)
        J.append([Jx, Jy, Jz])
        t.append(step*dt)
        nsteps += 1
    J = np.array(J)
    t = np.array(t)

    nwindows = int((nsteps - opts.nskip) // window)
    for i in range(3):
      for j in range(3):
        for k in range(nwindows):
          nruns += 1
          start = opts.nskip + k*window
          finish = opts.nskip + k*window + window
          G += autocorr(J[start:finish,i],J[start:finish,j])
  G = G / float(nruns)
  plt.figure("HFACF")
  plt.xlabel("t (fs)")
  plt.ylabel("HFACF")
  plt.xlim((0, time[-1]))
  plt.plot(time[:], G[0,0,:], 'r-', label='G_{xx}', linewidth=1.0)
  plt.plot(time[:], G[1,1,:], 'g-', label='G_{yy}', linewidth=1.0)
  plt.plot(time[:], G[2,2,:], 'b-', label='G_{zz}', linewidth=1.0)
  plt.plot(time[:], G[0,1,:], 'r--', label='G_{xy}', linewidth=1.0)
  plt.plot(time[:], G[0,2,:], 'b--', label='G_{xz}', linewidth=1.0)
  plt.plot(time[:], G[1,2,:], 'g--', label='G_{yz}', linewidth=1.0)
  plt.legend(loc='upper right')
  plt.savefig("hfacf.pdf", bbox_inches='tight')

# Parse the frames file

# Automatisation of the RDF datas
if opts.rdf and not opts.dynpro_graph_only:
  print("---Dynpro calculations---")
  print()
  print("#---------------------------------------------------------------------------#")
  print("Do you want the input for the RDF to be in function of your results?")
  print("It will make in sort that the cutoff will have the size of the smallest cell at the first frame divided by two.")
  print("If the smallest step is higher than half of the cutoff (entered or by default), the size of the step will be setup in function of the number of histogram bin you selected in the command line.")
  print("By default the number of histogram is 100 and the width of the histogram is 0.05 A.")
  while True:
    Shall_we = str(input("y, n or ?. "))
    if (Shall_we == "y" or Shall_we == "Y" or Shall_we == "Yes" or Shall_we == "YES" or Shall_we == "1"):
      Automatic_rdf = True
      break
    if (Shall_we == "n" or Shall_we == "N" or Shall_we == "No" or Shall_we == "NO" or Shall_we == "0"):
      Automatic_rdf = False
      print("#---------------------------------------------------------------------------#")
      break
    if (Shall_we == "?"):
      print("The size of the cutoff: "+str(opts.rdfcut)+" the command to input the value you want is: --rdfcut")
      print("The width of the histogram: "+str(opts.rdfwidth)+" the command to input the value you want is: --rdfwidth")
      print("The number of histogram: "+str(opts.nbins)+" the command to input the value you want is: --nbins")
if opts.vacf and not opts.rdf and not opts.dynpro_graph_only: # We make in sort that if we only want the VACF we will take the automatic values.
  Automatic_rdf = True


###################################
###      Dynpro_gestion.py      ###
##-------------------------------##
#dynpro_gestion.py is a script for#
# RDFs and VACF.                  #
# It will lauch a fortran script ##
#to evaluate the datas accurately.#
##-------------------------------##
###################################

if opts.rdf or opts.vacf: # Dynpro_gestion is only taken if we want an RDF or a VACF.
  if not opts.dynpro_graph_only:
    if read_frames:
      nframes = 0
      newframe = True
      buf = ""
      time = []
      stress = []
      lat = []
      first_frame = True
      done = False
      with open(opts.framesfile, 'r') as framesfile:
        while not done:
          line = framesfile.readline()
          if not line:
            break
          if re.match(frame_re, line):
            n = int(line.split()[1])

          if re.match(endframe_re, line):
            newframe = True
            if opts.snap != -1:
              if n != opts.snap:
                continue
              else:
                done = True
            if n <= opts.nskip:
              continue
            elif n%opts.stride != 0:
              continue
            else:
              nframes += 1
            if opts.nstop != -1:
              if n > opts.nstop:
                done = True
            sys.stdout.write("Processing frame {}\r".format(n))
            if first_frame:
              first_frame = False
              f1 = Frame(natoms,n)
              f1.parse_frame(buf)
              
              if Automatic_rdf and nframes <= 1: # We only take the first cell for consistancy.
                # The cutoff
                cell = np.zeros(3)
                for i in range(3):
                  cell[i] = f1.lat[i,i]
                rdfcutoff = min(cell)/2
                # The width of the step
                # Because the number of histogram doesn't need to be necessary put. We will first verify if the width is valid.
                # If it isn't, the number of histogram by default OR choosen will be taken.
                if (2*rdfcutoff < float(opts.rdfwidth)):
                  print("The width, "+str(opts.rdfwidth)+" is not valid.")
                  print("The number of histogram: "+str(opts.nbins)+" will be taken.")
                  rdfwidth = rdfcutoff/float(opts.nbins)
                else:
                  rdfwidth = float(opts.rdfwidth)
                # To inform the user of the variation of the values.
                print("As a result, we have:")
                print("A cutoff equal to: "+str(rdfcutoff))
                print("An histogram width equal to: "+str(rdfwidth))
                print("#---------------------------------------------------------------------------#")
              # The pair dist taken will be different if the automatisation is asked.
            
            f = Frame(natoms, n)
            f.parse_frame(buf)
            time.append(n*dt)
            continue
          if newframe:
            buf = ""
            newframe = False
          else:
            buf += line

      time = data['time']
      time = time - time[0]
      print()
      print("Analysing {} frames...".format(nframes))
      last_frame = nframes

    if float(opts.vacf_end) <= 0:
      ending_point = last_frame
    else:
      ending_point = opts.vacf_end
    All_atom_old        = []
    All_number_atom_old = []
    for ispec in range(init_config['nspecies']):
      Atomic_name_center_atom = cq_params['species'][ispec+1]

      Atom_name_old, Atom_name_new, Atom_number = Unknown_atom(cq_params['species'][ispec+1])
      
      if Atom_name_new != Atom_name_old:
        cq_params_ = parse_cq_input_edit(cq_input_file) # In case one atom isn't the good one.
      if Atom_name_new == Atom_name_old:
        cq_params_ = parse_cq_input(cq_input_file)      # In case all atoms are good.
      All_atom_old.append(Atom_name_old)
      All_number_atom_old.append(Atom_number)

    if opts.dynpro_substi:
      # ID_atom_replace
      print("What is the ID of the atom you want to modify?")
      try:
        while True:
          ID_atom_replace = int(input("-> "))
          if ID_atom_replace > init_config['natoms'] or ID_atom_replace <= 0:
            print("ID out of bonds.")
            print("Please, enter the ID again.")
            continue
          else:
            break
      except:
        while True:
          print("Invalid ID, please enter an integer.")
          try:
            ID_atom_replace = int(input("-> "))
          except:
            continue
          else:
            break
      # Input_Atom
      print("What is the placeholder name for this atom?")
      Input_Atom = str(input("-> "))

    for ispec in range(init_config['nspecies']):
      Atomic_name_center_atom = cq_params_['species'][ispec+1]
      Atom_name_old_useless, Atom_name_new, Atom_number = Unknown_atom(cq_params_['species'][ispec+1])
      for j in range(len(All_number_atom_old)):
        if (All_number_atom_old[j] == Atom_number):
          Atom_name_old = All_atom_old[j]
      dynpro = Dynpro_gestion(init_config['nspecies'],
                    rdfcutoff, rdfwidth,
                    cq_params_['species'], init_config['species_count'],
                    opts.vacf_start,ending_point,
                    opts.vacf_nbr_sample,cq_jump,Atomic_name_center_atom,
                    opts.r_min,cq_Temperature,opts.dynpro_md,opts.dynpro_in,
                    opts.framesfile,opts.statfile,opts.lcenter,opts.ct_atom,
                    opts.rdf,opts.vacf,opts.vacf_sample,opts.dynpro_substi,Atom_number,Atom_name_old)
      dynpro.File_generation(False)
      dynpro.dynpro_launch()
      if opts.dynpro_substi:
        dynpro.substitution(Input_Atom,ID_atom_replace)
        print(ispec+1)
        print(init_config['nspecies'])
        if ispec+1 == init_config['nspecies']:
          dynpro.File_generation(True)
          dynpro.dynpro_launch()
  else: # If we only want to plot the datas
      ending_point = opts.vacf_end
      rdfcutoff = opts.rdfcut
      rdfwidth  = opts.rdfwidth
      All_atom_old        = []
      All_number_atom_old = []
      for ispec in range(init_config['nspecies']):
        Atomic_name_center_atom = cq_params['species'][ispec+1]

        Atom_name_old, Atom_name_new, Atom_number = Unknown_atom(cq_params['species'][ispec+1])
        
        if Atom_name_new != Atom_name_old:
          cq_params_ = parse_cq_input_edit(cq_input_file) # In case one atom isn't the good one.
        if Atom_name_new == Atom_name_old:
          cq_params_ = parse_cq_input(cq_input_file)      # In case all atoms are good.
        All_atom_old.append(Atom_name_old)
        All_number_atom_old.append(Atom_number)

      for ispec in range(init_config['nspecies']):
        Atomic_name_center_atom = cq_params_['species'][ispec+1]
        Atom_name_old_useless, Atom_name_new, Atom_number = Unknown_atom(cq_params_['species'][ispec+1])
        for j in range(len(All_number_atom_old)):
          if (All_number_atom_old[j] == Atom_number):
            Atom_name_old = All_atom_old[j]
        dynpro = Dynpro_gestion(init_config['nspecies'],
                      rdfcutoff, rdfwidth,
                      cq_params_['species'], init_config['species_count'],
                      opts.vacf_start,ending_point,
                      opts.vacf_nbr_sample,cq_jump,Atomic_name_center_atom,
                      opts.r_min,cq_Temperature,opts.dynpro_md,opts.dynpro_in,
                      opts.framesfile,opts.statfile,opts.lcenter,opts.ct_atom,
                      opts.rdf,opts.vacf,opts.vacf_sample,opts.dynpro_substi,Atom_number,Atom_name_old)
  ######################
  # Plotting the datas #
  ######################
  if opts.dynpro_substi:
    dynpro.custom_name()

  if ( init_config['nspecies'] > 1):

    # Loop to choose the number of graphs to plot.
    print("*-----------------------------*")
    print("How many graphics do you want?") 
    print("Write 0 or ? if you don't know.")
    print("Write r if you want to print everything without seing the graphics.")
    try:
      Number_of_graph = input("-> ")
      if (str(Number_of_graph) == "?"):
        Number_of_graph = 0
      if (str(Number_of_graph) == "r" or str(Number_of_graph) == "R"):
        Number_of_graph = -1
      Number_of_graph = int(Number_of_graph)
      if Number_of_graph < -1:
            Number_of_graph = abs(Number_of_graph)
    except:
      while( type(Number_of_graph) is not int):
        try:
          Number_of_graph = input("Unsupported number, please try again. ")
          if (str(Number_of_graph) == "?"):
            Number_of_graph = 0
          if (str(Number_of_graph) == "r" or str(Number_of_graph) == "R"):
            Number_of_graph = -1
          Number_of_graph = int(Number_of_graph)
          if Number_of_graph < -1:
            Number_of_graph = abs(Number_of_graph)
        except:
          continue
    if Number_of_graph == -1:
      Answer_nbr_graph = "All the"
    elif Number_of_graph == 0:
      Answer_nbr_graph = "Unknown number of"
    else:
      Answer_nbr_graph = Number_of_graph
    print(str(Answer_nbr_graph)+" graphic(s) will be generated.")
    loop = 0
    while((loop < Number_of_graph or Number_of_graph == 0) and Number_of_graph != -1):
      selection, Central_atom_to_print = Atom_selector(init_config['nspecies'],cq_params_['species'])
      if selection == 2:
        if opts.vacf:
          show_VACF = True
        else:
          show_VACF = False
        if opts.rdf:
          show_RDF  = True
        else:
          show_RDF  = False
      if selection == 1:
        print("All the graphics will be generated.")
        if opts.vacf:
          show_VACF = logic_selector("VACF")
        else:
          show_VACF = False
        if opts.vacf:
          show_RDF = logic_selector("RDF")
        else:
          show_RDF = False
      dynpro.plot_dynpro(show_VACF,show_RDF,Central_atom_to_print,cq_params['species'],opts.advanced_rdf)

      loop = loop+1
      if (Number_of_graph == 0):
        quitt_plot = str(input("Do you want to to quitt? y or n. "))
        if (quitt_plot == "y" or quitt_plot == "Y" or quitt_plot == "Yes" or quitt_plot == "YES" or quitt_plot == "1"):
          break
    if (Number_of_graph == -1):
      print("All the graphics will be generated.")
      dynpro.plot_dynpro(False,False,"all",cq_params['species'],opts.advanced_rdf)
  else:
    Central_atom_to_print = cq_params['species'][1]
    if opts.vacf:
      show_VACF = True
    else:
      show_VACF = False
    if opts.rdf:
      show_RDF  = True
    else:
      show_RDF  = False
    dynpro.plot_dynpro(show_VACF,show_RDF,Central_atom_to_print,cq_params['species'],opts.advanced_rdf)

if opts.advanced_rdf:
  print("****")
  print("The advanced options for the RDF is used.")
  print("If you want to redo the calculation, then, delete the files generated previously.")
  print("****")

###########################
# MSD and stress analysis #
###########################

if read_frames and (opts.msd or not opts.stress):
  if opts.msd and not opts.stress :
    print("---MSD calculations---")
  if opts.msd and opts.stress :
    print("---MSD and stress calculations---")
  if not opts.msd and opts.stress :
    print("---Stress calculations---") 
  nframes = 0
  newframe = True
  buf = ""
  time = []
  stress = []
  lat = []
  first_frame = True
  done = False
  with open(opts.framesfile, 'r') as framesfile:
    while not done:
      line = framesfile.readline()
      if not line:
        break
      if re.match(frame_re, line):
        n = int(line.split()[1])

      if re.match(endframe_re, line):
        newframe = True
        if opts.snap != -1:
          if n != opts.snap:
            continue
          else:
            done = True
        if n <= opts.nskip:
          continue
        elif n%opts.stride != 0:
          continue
        else:
          nframes += 1
        if opts.nstop != -1:
          if n > opts.nstop:
            done = True
        sys.stdout.write("Processing frame {}\r".format(n))
        if first_frame:
          first_frame = False
          f1 = Frame(natoms,n)
          f1.parse_frame(buf)
          
          if opts.msd:
            if opts.msd_multi_ID == 0 :
              if opts.msd_skip >= 1: # If the skip as been set up, then we can put it inside the variable.
                time_skip = opts.msd_skip
              else: # Else we put a skip of 0 (fs).
                time_skip = 0
              m = MSD(natoms, dt, f1,
                      opts.msd_multi_ID,time_skip,cq_jump)
            if opts.msd_multi_ID > 0 :
              if opts.msd_skip >= 1: # If the skip as been set up, then we can put it inside the variable.
                time_skip = opts.msd_skip
              else: # Else we put a skip of 0 (fs).
                time_skip = 0
              m = MSD(natoms, dt, f1,
                      opts.msd_multi_ID,time_skip,cq_jump)

        f = Frame(natoms, n)
        f.parse_frame(buf)

        time.append(n*dt)
        if opts.stress:
          stress.append(f.stress)
          lat.append(f.lat)
        if opts.msd:
          m.update_msd(n, f)
        continue
      if newframe:
        buf = ""
        newframe = False
      else:
        buf += line

  time = data['time']
  time = time - time[0]
  print()
  print("Analysing {} frames...".format(nframes))

# Plot the stress
  if opts.stress:
    stress = np.array(stress)
    lat    = np.array(lat)
    mean_stress = np.zeros((3,3))
    mean_lat    = np.zeros((3,3))
    for i in range(3):
      for j in range(3):
        mean_stress[i,j] = sp.mean(stress[:,i,j])
        mean_lat[i,j] = sp.mean(lat[:,i,j])
    plt.figure("Stress")

    if cq_params['MD.Ensemble'][1] == "p":
      fig2, (ax1, ax2) = plt.subplots(nrows=2, ncols=1, sharex=True)
    else:
      fig2, (ax1,) = plt.subplots(nrows=1, ncols=1)

    plt.xlabel("t (fs)")
    ax1.set_ylabel("Stress (GPa)")
    ax2.set_ylabel("Cell dimension ($a_0$)")
    plt.xlim((time[opts.nskip], time[-1]))
    ax1.plot(time[opts.nskip:], stress[:,0,0], 'r-', label='xx', linewidth=1.0)
    ax1.plot(time[opts.nskip:], stress[:,1,1], 'g-', label='yy', linewidth=1.0)
    ax1.plot(time[opts.nskip:], stress[:,2,2], 'b-', label='zz', linewidth=1.0)
    ax1.plot((time[0],time[-1]), (mean_stress[0,0], mean_stress[0,0]), 'r-',
            label=r'$\langle S_{{xx}} \rangle$ = {0:<10.4f}'.format(mean_stress[0,0]))
    ax1.plot((time[0],time[-1]), (mean_stress[1,1], mean_stress[1,1]), 'g-',
            label=r'$\langle S_{{yy}} \rangle$ = {0:<10.4f}'.format(mean_stress[1,1]))
    ax1.plot((time[0],time[-1]), (mean_stress[2,2], mean_stress[2,2]), 'b-',
            label=r'$\langle S_{{zz}} \rangle$ = {0:<10.4f}'.format(mean_stress[2,2]))

    if cq_params['MD.Ensemble'][1] == "p":
      ax2.plot(time[opts.nskip:], lat[:,0,0], 'r-', label='a', linewidth=1.0)
      ax2.plot(time[opts.nskip:], lat[:,1,1], 'g-', label='b', linewidth=1.0)
      ax2.plot(time[opts.nskip:], lat[:,2,2], 'b-', label='c', linewidth=1.0)
      ax2.plot((time[0],time[-1]), (mean_lat[0,0], mean_lat[0,0]), 'r-',
              label=r'$\langle a \rangle$ = {0:<10.4f}'.format(mean_lat[0,0]))
      ax2.plot((time[0],time[-1]), (mean_lat[1,1], mean_lat[1,1]), 'g-',
              label=r'$\langle b \rangle$ = {0:<10.4f}'.format(mean_lat[1,1]))
      ax2.plot((time[0],time[-1]), (mean_lat[2,2], mean_lat[2,2]), 'b-',
              label=r'$\langle c \rangle$ = {0:<10.4f}'.format(mean_lat[2,2]))
      ax1.legend(bbox_to_anchor=(1.05,1), loc=2, borderaxespad=0.)
      ax2.legend(bbox_to_anchor=(1.05,1), loc=2, borderaxespad=0.)
      fig2.subplots_adjust(hspace=0)
      plt.setp([a.get_xticklabels() for a in fig1.axes[:-1]], visible=False)
      fig2.savefig("stress.pdf", bbox_inches='tight')

# Plot the MSD
  if opts.msd:
    m.norm_msd()
    m.Diffusion_msd()
    if opts.msd_multi_ID > 0:
      m.Diffusion_msd_multiple()
      m.plot_msd_multiple()
    m.plot_msd()
    if opts.dump:
      m.dump_msd()

plt.show()
