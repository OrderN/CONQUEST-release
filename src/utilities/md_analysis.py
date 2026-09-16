#!/usr/local/bin/python3

import argparse
import sys
import re
import os.path
import numpy as np
import matplotlib.pyplot as plt
from frame import Frame
from md_tools import Pairdist, MSER, VACF, MSD, autocorr

ha2ev = 27.211399
ha2k = 3.15737513e5
ha_bohr3_to_gpa = 29421.01549104606

# Regular expressions
frame_re = re.compile('frame')
endframe_re = re.compile('end frame')
specblock_re = re.compile(r'%block ChemicalSpeciesLabel\n(.*?)\n%endblock',
                          re.M | re.S | re.I)

cq_input_file = 'Conquest_input'

# Parsing functions
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
  with open(cq_input_file, 'r') as cqip:
    m = re.search(specblock_re, cqip.read())
    specinfo = m.group(1).splitlines()
    for line in specinfo:
      try:
        index, mass, spec  = line.split()
      except ValueError:
        index, mass, spec, ionfile = line.split()
      cq_params['species'][int(index)] = spec

  return cq_params

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
    data['coords'] = np.array(coords)
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
  data = {}
  header = True
  with open(stats_file, 'r') as statfile:
    for line in statfile:
      if header:
        col_id = line.strip().split()
        for col in col_id:
          data[col] = []
        header = False
      else:
        bits = line.strip().split()
        step = int(bits[0])
        if nstop != -1 and step > nstop:
          break
        for i, bit in enumerate(bits):
          if i==0:
            info = step
          else:
            info = float(bit)
          data[col_id[i]].append(info)
    for key in data:
      data[key] = np.array(data[key])
  return len(data['step']), data

# Command line arguments
parser = argparse.ArgumentParser(description='Analyse a Conquest MD \
        trajectory', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-c', '--compare', action='store_true', default=False,
                    dest='compare', help='Compare statistics of trajectories \
                    in directories specified by -d')
parser.add_argument('-d', '--dirs', nargs='+', default=['.',], dest='dirs',
                    action='store', help='Directories to compare')
parser.add_argument('--description', nargs='+', default=[], dest='desc',
                    action='store', help='Description of graph for legend \
                    (only if using --compare)')
parser.add_argument('-f', '--frames', action='store', dest='framesfile',
                    default='md.frames', help='MD frames file')
parser.add_argument('-s', '--stats-file', action='store', dest='statfile',
                    default='md.stats', help='MD statistics file')
parser.add_argument('--heatflux-file', action='store', dest='heatfluxfile',
                    default='md.heatflux', help='MD heat-flux file')
parser.add_argument('--skip', action='store', dest='nskip', default=0,
                    type=int, help='Number of equilibration steps to skip')
parser.add_argument('--stride', action='store', dest='stride', default=1,
                    type=int, help='Only analyse every nth step of frames file')
parser.add_argument('--snap', action='store', dest='snap', default=-1, 
                    type=int, help='Analyse Frame of a single snapshot')
parser.add_argument('--stop', action='store', dest='nstop', default=-1, 
                    type=int, help='Number of last frame in analysis')
parser.add_argument('--equil', action='store', dest='nequil', default=0, 
                    type=int, help='Number of equilibration steps')
parser.add_argument('--stats', action='store_true', dest='stats',
                    help='Plot statistics')
parser.add_argument('--vacf', action='store_true', dest='vacf', 
                    help='Plot velocity autocorrelation function')
parser.add_argument('--hfacf', action='store_true', dest='hfacf', 
                    help='Plot heat flux autocorrelation function')
parser.add_argument('--acfwindow', action='store', dest='acfwindow', default=0.0, 
                    type=float, help='window for autocorrelation')
parser.add_argument('--msd', action='store_true', dest='msd', 
                    help='Plot mean squared deviation')
parser.add_argument('--rdf', action='store_true', dest='rdf', 
                    help='Plot radial distribution function')
parser.add_argument('--stress', action='store_true', dest='stress', 
                    help='Plot stress')
parser.add_argument('--landscape', action='store_true', dest='landscape', 
                    help='Generate plot with landscape orientation')
parser.add_argument('--pub', action='store_true', dest='pub', 
                    help='Publication text size')
rdf_bins = parser.add_mutually_exclusive_group()
rdf_bins.add_argument('--nbins', action='store', dest='nbins', default=100,
                      type=int, help='Number of RDF histogram bins')
rdf_bins.add_argument('--rdfwidth', action='store', dest='rdfwidth',
                      default=None, type=float,
                      help='RDF histogram bin width (A)')
parser.add_argument('--rdfcut', action='store', dest='rdfcut', default=10.0,
                    type=float,
                    help='Distance cutoff for RDF')
parser.add_argument('--dump', action='store_true', dest='dump', 
                    help='Dump secondary data used to generate plots')
parser.add_argument('--mser', action='store', dest='mser_var', default=None,
                    type=str, help='Compute MSER for the given property')

opts = parser.parse_args()
if (opts.vacf or opts.msd or opts.stress or opts.rdf):
  read_frames = True
else:
  read_frames = False

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
  # Parse the statistics file
  nsteps, data = read_stats(opts.statfile,opts.nstop)
  avg = {}
  std = {}
  for key in data:
    data[key] = np.array(data[key])
  time = [float(s)*dt for s in data['step']]
  data['time'] = np.array(time)
  plot_mask = data['step'] >= opts.nskip
  equil_mask = data['step'] >= opts.nequil
  for key in data:
    avg[key] = np.mean(data[key][equil_mask])
    std[key] = np.std(data[key][equil_mask])
  plot_start = data['time'][plot_mask][0]

  # Plot the statistics
  if opts.landscape:
    fig1, ((ax1, ax2), (ax3, ax4)) = plt.subplots(nrows=2, ncols=2, sharex=True, figsize=(11,7))
    plt.tight_layout(pad=6.5)
  else:
    fig1, (ax1, ax2, ax3, ax4) = plt.subplots(nrows=4, ncols=1, sharex=True, figsize=(7,10))

  ax1.plot(data['time'][plot_mask], data['pe'][plot_mask], 'r-', label='Potential energy')
  ax1a = ax1.twinx()
  ax1a.plot(data['time'][plot_mask], data['ke'][plot_mask], 'b-', label='Kinetic energy')
  if cq_params['MD.Ensemble'][2] == 't':
    if cq_params['MD.Thermostat'] == 'nhc':
      ax1a.plot(data['time'][plot_mask], data['thermostat'][plot_mask], 'g-', label='Thermostat energy')
    if cq_params['MD.Thermostat'] == 'svr':
      ax1a.plot(data['time'][plot_mask], data['thermostat'][plot_mask], 'g-', label='Thermostat energy')
  if cq_params['MD.Ensemble'][1] == 'p':
    if 'barostat' in data:
      ax1a.plot(data['time'][plot_mask], data['barostat'][plot_mask], 'c-', label='Barostat energy')
    ax1a.plot(data['time'][plot_mask], data['pV'][plot_mask], 'm-', label='pV')
  ax2.plot(data['time'][plot_mask], data['H\''][plot_mask])
  ax2.plot((plot_start,data['time'][-1]), (avg['H\''],avg['H\'']), '-',
        label=r'$\langle H\' \rangle$ = {0:>12.4f} $\pm$ {1:<12.4f}'.format(avg['H\''], std['H\'']))
  ax3.plot(data['time'][plot_mask], data['T'][plot_mask])
  ax3.plot((plot_start,data['time'][-1]), (avg['T'],avg['T']), '-',
        label=r'$\langle T \rangle$ = {0:>12.4f} $\pm$ {1:<12.4f}'.format(avg['T'], std['T']))
  ax4.plot(data['time'][plot_mask], data['P'][plot_mask], 'b-')
  ax4.plot((plot_start,data['time'][-1]), (avg['P'],avg['P']), 'b--',
        label=r'$\langle P \rangle$ = {0:>12.4f} $\pm$ {1:<12.4f}'.format(avg['P'], std['P']))
  if cq_params['MD.Ensemble'][1] == 'p':
    ax4a = ax4.twinx()
    ax4a.plot(data['time'][plot_mask], data['V'][plot_mask], 'r-')
    ax4a.plot((plot_start,data['time'][-1]), (avg['V'],avg['V']), 'r--',
              label=r'$\langle V \rangle$ = {0:>12.4f} $\pm$ {1:<12.4f}'.format(avg['V'], std['V']))
  ax1.set_ylabel("E (Ha)")
  ax2.set_ylabel("H$'$ (Ha)")
  ax3.set_ylabel("T (K)")
  ax4.set_ylabel("P (GPa)", color='b')
  if cq_params['MD.Ensemble'][1] == 'p':
    ax4a.set_ylabel("V ($a_0^3$)", color='r')
  ax4.set_xlabel("time (fs)")
  ax1.legend(loc="upper left")
  ax1a.legend(loc="lower right")
  ax2.legend()
  ax3.legend()
  ax4.legend(loc="upper left")
  if cq_params['MD.Ensemble'][1] == 'p':
    ax4a.legend(loc="lower right")
  plt.xlim((plot_start,data['time'][-1]))
  fig1.subplots_adjust(hspace=0)
  fig1.savefig("stats.pdf", bbox_inches='tight')
else:
  # If we're comparing statistics in several directories, use a simplified plot
  fig1, (ax1, ax2, ax3) = plt.subplots(nrows=3, ncols=1, sharex=True, figsize=(7,7))
  ax1a = ax1.twinx()
  labels = opts.desc if opts.desc else opts.dirs
  if len(labels) != len(opts.dirs):
    parser.error('--description must provide one label per comparison directory')
  time_limits = []
  for ind, d in enumerate(opts.dirs):
    path = os.path.join(d, cq_input_file)
    cq_params = parse_cq_input(path)
    dt = float(cq_params['AtomMove.Timestep'])
  
    path = os.path.join(d, opts.statfile)
    nsteps, data = read_stats(path,opts.nstop)
    time = [float(s)*dt for s in data['step']]
    data['time'] = np.array(time)
    plot_mask = data['step'] >= opts.nskip
    time_limits.append((data['time'][plot_mask][0], data['time'][-1]))

    ax1.plot(data['time'][plot_mask], data['H\''][plot_mask],
             linewidth=0.5, label=labels[ind])
    ax2.plot(data['time'][plot_mask], data['T'][plot_mask],
             linewidth=0.5, label=labels[ind])
    ax3.plot(data['time'][plot_mask], data['P'][plot_mask],
             linewidth=0.5, label=labels[ind])

  y1,y2 = ax1.get_ylim()
  ax1a.set_ylim(y1*ha2k,y2*ha2k)
  ax1.set_ylabel("H$'$ (Ha)")
  ax1a.set_ylabel("H$'$ (K)")
  ax2.set_ylabel("T (K)")
  ax3.set_ylabel("P (GPa)")
  ax3.set_xlabel("time (fs)")
  ax1.legend()
  plt.xlim((min(limit[0] for limit in time_limits),
            max(limit[1] for limit in time_limits)))
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
  flux_runs = []
  flux_dt = None
  for d in opts.dirs:
    run_params = parse_cq_input(os.path.join(d, cq_input_file))
    run_dt = float(run_params['AtomMove.Timestep'])
    if flux_dt is None:
      flux_dt = run_dt
    elif not np.isclose(run_dt, flux_dt):
      parser.error('all heat-flux trajectories must use the same timestep')
    path = os.path.join(d, opts.heatfluxfile)
    J = []
    with open(path, 'r') as infile:
      for line in infile:
        step, Jx, Jy, Jz = line.split()
        step = int(step)
        if step < opts.nskip:
          continue
        if opts.nstop != -1 and step > opts.nstop:
          break
        J.append([float(Jx), float(Jy), float(Jz)])
    if J:
      flux_runs.append(np.array(J))

  if not flux_runs:
    parser.error('no heat-flux samples remain after step selection')
  if opts.acfwindow > 0.0:
    window = int(opts.acfwindow // flux_dt)
  else:
    window = min(len(run) for run in flux_runs)
  if window < 1 or any(len(run) < window for run in flux_runs):
    parser.error('autocorrelation window exceeds the available heat-flux data')

  G = np.zeros((3,3,window))
  nruns = 0
  for J in flux_runs:
    for start in range(0, len(J)-window+1, window):
      sample = J[start:start+window]
      for i in range(3):
        for j in range(3):
          G[i,j,:] += autocorr(sample[:,i], sample[:,j])
      nruns += 1
  if nruns == 0:
    parser.error('heat-flux data do not contain a complete correlation window')
  G /= float(nruns)
  time = np.arange(window, dtype=float)*flux_dt
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
        if n < opts.nskip:
          continue
        elif n%opts.stride != 0:
          continue
        if opts.nstop != -1 and n > opts.nstop:
          done = True
          continue
        nframes += 1
        sys.stdout.write("Processing frame {}\r".format(n))
        if first_frame:
          first_frame = False
          f1 = Frame(natoms,n)
          f1.parse_frame(buf)
          if opts.rdf:
            if opts.rdfwidth is None:
              rdfwidth = opts.rdfcut/opts.nbins
            else:
              rdfwidth = opts.rdfwidth
            pairdist = Pairdist(natoms, init_config['nspecies'],
                                opts.rdfcut, rdfwidth,
                                cq_params['species'],
                                init_config['species_count'])
          if opts.vacf:
            c = VACF(natoms, dt, f1)
          if opts.msd:
            m = MSD(natoms, dt, f1)

        f = Frame(natoms, n)
        f.parse_frame(buf)

        time.append(n*dt)
        if opts.stress:
          stress.append(f.stress)
          lat.append(f.lat)
        if opts.rdf:
          pairdist.update_rdf(f)
        if opts.vacf:
          c.update_vacf(n, f)
        if opts.msd:
          m.update_msd(n, f)
        continue
      if newframe:
        buf = ""
        newframe = False
      else:
        buf += line

  time = np.array(time)
  time = time - time[0]
  print()
  print("Analysing {} frames...".format(nframes))

# Plot the stress
  if opts.stress:
    stress = np.array(stress)*ha_bohr3_to_gpa
    lat = np.array(lat)
    mean_stress = np.mean(stress, axis=0)
    mean_lat = np.mean(lat, axis=0)
    variable_cell = cq_params['MD.Ensemble'][1] == "p"

    if variable_cell:
      fig2, (ax1, ax2) = plt.subplots(nrows=2, ncols=1, sharex=True)
    else:
      fig2, ax1 = plt.subplots(nrows=1, ncols=1)

    ax1.set_ylabel("Stress (GPa)")
    ax1.set_xlim((time[0], time[-1]))
    ax1.plot(time, stress[:,0,0], 'r-', label='xx', linewidth=1.0)
    ax1.plot(time, stress[:,1,1], 'g-', label='yy', linewidth=1.0)
    ax1.plot(time, stress[:,2,2], 'b-', label='zz', linewidth=1.0)
    ax1.plot((time[0],time[-1]), (mean_stress[0,0], mean_stress[0,0]), 'r-',
            label=r'$\langle S_{{xx}} \rangle$ = {0:<10.4f}'.format(mean_stress[0,0]))
    ax1.plot((time[0],time[-1]), (mean_stress[1,1], mean_stress[1,1]), 'g-',
            label=r'$\langle S_{{yy}} \rangle$ = {0:<10.4f}'.format(mean_stress[1,1]))
    ax1.plot((time[0],time[-1]), (mean_stress[2,2], mean_stress[2,2]), 'b-',
            label=r'$\langle S_{{zz}} \rangle$ = {0:<10.4f}'.format(mean_stress[2,2]))

    ax1.legend(bbox_to_anchor=(1.05,1), loc=2, borderaxespad=0.)
    if variable_cell:
      ax2.set_ylabel("Cell dimension ($a_0$)")
      ax2.set_xlabel("t (fs)")
      ax2.plot(time, lat[:,0,0], 'r-', label='a', linewidth=1.0)
      ax2.plot(time, lat[:,1,1], 'g-', label='b', linewidth=1.0)
      ax2.plot(time, lat[:,2,2], 'b-', label='c', linewidth=1.0)
      ax2.plot((time[0],time[-1]), (mean_lat[0,0], mean_lat[0,0]), 'r-',
              label=r'$\langle a \rangle$ = {0:<10.4f}'.format(mean_lat[0,0]))
      ax2.plot((time[0],time[-1]), (mean_lat[1,1], mean_lat[1,1]), 'g-',
              label=r'$\langle b \rangle$ = {0:<10.4f}'.format(mean_lat[1,1]))
      ax2.plot((time[0],time[-1]), (mean_lat[2,2], mean_lat[2,2]), 'b-',
              label=r'$\langle c \rangle$ = {0:<10.4f}'.format(mean_lat[2,2]))
      ax2.legend(bbox_to_anchor=(1.05,1), loc=2, borderaxespad=0.)
      fig2.subplots_adjust(hspace=0)
    else:
      ax1.set_xlabel("t (fs)")
    fig2.savefig("stress.pdf", bbox_inches='tight')

  # Plot the rdf
  if opts.rdf:
    pairdist.norm_rdf()
    pairdist.get_coordination()
    pairdist.plot_gr()
    if opts.dump:
      pairdist.dump_gr()

  # Plot the VACF
  if opts.vacf:
    c.norm_vacf()
    c.plot_vacf()
    if opts.dump:
      c.dump_vacf()

  # Plot the MSD
  if opts.msd:
    m.norm_msd()
    m.plot_msd()
    if opts.dump:
      m.dump_msd()
