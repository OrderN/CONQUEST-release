#!/usr/local/bin/python3

import argparse
import re
import scipy as sp
import matplotlib.pyplot as plt
from frame import Frame
from pdb import set_trace

cq_input_file = 'Conquest_input'

# Regular expressions
frame_re = re.compile('frame')
endframe_re = re.compile('end frame')
cell_re = re.compile('cell_vectors(.*?)end cell_vectors', re.M | re.S)
stress_re = re.compile('stress_tensor(.*?)end stress_tensor', re.M | re.S)
position_re = re.compile('positions(.*?)end positions', re.M | re.S)
position_re = re.compile('positions(.*?)end positions', re.M | re.S)
velocity_re = re.compile('velocities(.*?)end velocities', re.M | re.S)
force_re = re.compile('forces(.*?)end forces', re.M | re.S)

# Parsing functions
def parse_cq_input(cq_input_file):
  cq_params = {}
  with open(cq_input_file, 'r') as cqip:
    for line in cqip:
      if (line.strip() and (line.strip()[0] != '#') and (line.strip()[0] != '%')):
        bits = line.strip().split()
        if len(bits[1:]) == 1:
          cq_params[bits[0]] = bits[1:][0]
        else:
          cq_params[bits[0]] = bits[1:]
  return cq_params

def parse_init_config(conf_filename):
  data = {}
  with open(conf_filename, 'r') as infile:
    a = [float(bit) for bit in infile.readline().strip().split()]
    b = [float(bit) for bit in infile.readline().strip().split()]
    c = [float(bit) for bit in infile.readline().strip().split()]
    data['latvec'] = sp.array([a,b,c])
    natoms = int(infile.readline().strip())
    data['natoms'] = natoms
    coords = []
    species = []
    for i in range(natoms):
      x, y, z, spec, cx, cy, cz = infile.readline().strip().split()
      coords.append([float(x), float(y), float(z)])
      species.append(int(spec))
    data['coords'] = sp.array(coords)
    data['species'] = sp.array(species)
  return data

def read_stats(stats_file, nstop):
  nstep = 0
  data = {}
  header = True
  with open(stats_file, 'r') as statfile:
    for line in statfile:
      if nstop != -1:
        if nstep >= nstop:
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
  return data


def parse_frame(buf, f):
  m = re.search(cell_re, buf)
  lines = m.group(1).strip().splitlines()
  for i in range(3):
    bits = lines[i].strip().split()
    for j in range(3):
      f.lat[i,j] = float(bits[j])

  m = re.search(stress_re, buf)
  if m:
    lines = m.group(1).strip().splitlines()
    for i in range(3):
      bits = lines[i].strip().split()
      for j in range(3):
        f.stress[i,j] = float(bits[j])

  m = re.search(position_re, buf)
  lines = m.group(1).strip().splitlines()
  nat = len(lines)
  for i in range(nat):
    bits = lines[i].strip().split()
    bits.pop(0)
    f.species[i] = int(bits.pop(0))
    for j in range(3):
      f.r[i,j] = float(bits[j])

  m = re.search(velocity_re, buf)
  lines = m.group(1).strip().splitlines()
  nat = len(lines)
  for i in range(nat):
    bits = lines[i].strip().split()
    bits.pop(0)
    bits.pop(0)
    for j in range(3):
      f.v[i,j] = float(bits[j])

  m = re.search(force_re, buf)
  lines = m.group(1).strip().splitlines()
  nat = len(lines)
  for i in range(nat):
    bits = lines[i].strip().split()
    bits.pop(0)
    bits.pop(0)
    for j in range(3):
      f.f[i,j] = float(bits[j])

# Command line arguments
parser = argparse.ArgumentParser(description='Analyse a Conquest MD \
        trajectory', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-d', '--dirs', nargs='+', default='.', dest='dirs',
                    action='store', help='Directories to analyse')
parser.add_argument('-f', '--frames', action='store', dest='framesfile',
                    default='Frames', help='MD frames file')
parser.add_argument('-s', '--stats', action='store', dest='statfile',
                    default='Stats', help='MD statistics file')
parser.add_argument('--skip', action='store', dest='nskip', default=0,
                    type=int, help='Number of equilibration steps to skip')
parser.add_argument('--stop', action='store', dest='nstop', default=-1, 
                    type=int, help='Number of last frame in analysis')
parser.add_argument('--equil', action='store', dest='nequil', default=0, 
                    type=int, help='Number of equilibration steps')
parser.add_argument('--vacf', action='store_true', dest='vacf', 
                    help='Compute the velocity autocorrelation function')
parser.add_argument('--msd', action='store_true', dest='msd', 
                    help='Compute the mean squared deviation')
parser.add_argument('--vdist', action='store_true', dest='vdist', 
                    help='Compute the velocity distribution')
parser.add_argument('--landscape', action='store_true', dest='landscape', 
                    help='Generate plot with landscape orientation')
parser.add_argument('--stress', action='store_true', dest='stress', 
                    help='Plot the stress')
parser.add_argument('--nbins', action='store', dest='nbins', default=100,
                    help='Number of histogram bins')

opts = parser.parse_args()
if (opts.vacf or opts.msd or opts.vdist or opts.stress):
  read_frames = True
else:
  read_frames = False

if opts.nskip > 0:
  opts.nequil = opts.nskip

# Parse the md.in parameters file
cq_params = parse_cq_input(cq_input_file)
init_config = parse_init_config(cq_params['IO.Coordinates'])
natoms = init_config['natoms']
dt = float(cq_params['AtomMove.Timestep'])

# Parse the statistics file
data = read_stats(opts.statfile,opts.nstop)
avg = {}
std = {}
for key in data:
  data[key] = sp.array(data[key])
  avg[key] = sp.mean(data[key][opts.nequil:-1])
  std[key] = sp.std(data[key][opts.nequil:-1])
time = [float(s)*dt for s in data['step']]
data['time'] = sp.array(time)

# Plot the statistics
if opts.landscape:
  fig1, ((ax1, ax2), (ax3, ax4)) = plt.subplots(nrows=2, ncols=2, sharex=True, figsize=(11,7))
  plt.tight_layout(pad=6.5)
else:
  fig1, (ax1, ax2, ax3, ax4) = plt.subplots(nrows=4, ncols=1, sharex=True, figsize=(7,10))

ax1.plot(data['time'], data['pe'], 'r-', label='Potential energy')
ax1a = ax1.twinx()
ax1a.plot(data['time'], data['ke'], 'b-', label='Kinetic energy')
if cq_params['MD.Ensemble'][2] == 't':
  if cq_params['MD.Thermostat'] == 'nhc':
    ax1a.plot(data['time'], data['nhc'], 'g-', label='NHC energy')
if cq_params['MD.Ensemble'][1] == 'p':
  if 'mttk' in cq_params['MD.Barostat']:
    ax1a.plot(data['time'], data['box'], 'c-', label='Box energy')
  ax1a.plot(data['time'], data['pV'], 'm-', label='pV')
ax2.plot(data['time'], data['H\''])
ax2.plot((opts.nskip,data['time'][-1]), (avg['H\''],avg['H\'']), '-',
      label=r'$\langle H\' \rangle$ = {0:>12.4f} $\pm$ {1:<12.4f}'.format(avg['H\''], std['H\'']))
ax3.plot(data['time'], data['T'])
ax3.plot((opts.nskip,data['time'][-1]), (avg['T'],avg['T']), '-',
      label=r'$\langle T \rangle$ = {0:>12.4f} $\pm$ {1:<12.4f}'.format(avg['T'], std['T']))
ax4.plot(data['time'], data['P'], 'b-')
ax4.plot((opts.nskip,data['time'][-1]), (avg['P'],avg['P']), 'b--',
      label=r'$\langle P \rangle$ = {0:>12.4f} $\pm$ {1:<12.4f}'.format(avg['P'], std['P']))
if cq_params['MD.Ensemble'][1] == 'p':
  ax4a = ax4.twinx()
  ax4a.plot(data['time'], data['V'], 'r-')
  ax4a.plot((opts.nskip,data['time'][-1]), (avg['V'],avg['V']), 'r--',
        label=r'$\langle V \rangle$ = {0:>12.4f} $\pm$ {1:<12.4f}'.format(avg['V'], std['V']))
ax1.set_ylabel("E")
ax2.set_ylabel("E")
ax3.set_ylabel("T")
ax4.set_ylabel("P", color='b')
if cq_params['MD.Ensemble'][1] == 'p':
  ax4a.set_ylabel("V", color='r')
ax4.set_xlabel("time (fs)")
ax1.legend(loc="upper left")
ax1a.legend(loc="lower right")
ax2.legend()
ax3.legend()
ax4.legend(loc="upper left")
if cq_params['MD.Ensemble'][1] == 'p':
  ax4a.legend(loc="lower right")
plt.xlim((opts.nskip,data['time'][-1]))
fig1.subplots_adjust(hspace=0)
fig1.savefig("stats.pdf", bbox_inches='tight')

if read_frames:
  nframes = 0
  newframe = True
  buf = ""
  time = []
  stress = []
  vacf = []
  msd = []
  sdistr = sp.zeros(opts.nbins, dtype='float')
  vdistr = sp.zeros((opts.nbins,3), dtype='float')
  with open(opts.framesfile, 'r') as framesfile:
    while True:
      line = framesfile.readline()
      if not line:
        break
      if re.match(frame_re, line):
        n = int(line.split()[1])

      if re.match(endframe_re, line):
        newframe = True
        if n < opts.nskip:
          continue
        if n == opts.nstop:
          break
        else:
          nframes += 1
        if nframes == 1:
          f1 = Frame(natoms,n)
          parse_frame(buf,f1)
        f = Frame(natoms, n)
        parse_frame(buf, f)

        time.append(n*dt)
        if opts.stress:
          stress.append(f.stress)
        if opts.vacf:
          vacf.append(f.update_vacf(f1))
        if opts.msd:
          msd.append(f.update_msd(f1))
        if opts.vdist:
          sdistr_tmp, bin_edges = f.update_sdistr(opts.nbins)
          vdistr_tmp, bin_edges = f.update_vdistr(opts.nbins)
          sdistr += sdistr_tmp
          vdistr += vdistr_tmp
        continue
      if newframe:
        buf = ""
        newframe = False
      else:
        buf += line

  time = sp.array(time)
  time = time - time[0]

# Plot the stress
  if opts.stress:
    stress = sp.array(stress)
    mean_stress = sp.zeros((3,3))
    for i in range(3):
      for j in range(3):
        mean_stress[i,j] = sp.mean(stress[:,i,j])
    plt.figure("Stress")
    fig2, (ax1, ax2) = plt.subplots(nrows=2, ncols=1, sharex=True)
    plt.xlabel("t")
    ax1.set_ylabel("Stress")
    ax2.set_ylabel("Pressure")
    plt.xlim((time[0], time[-1]))
    ax1.plot(time, stress[:,0,0], 'r-', label='xx', linewidth=1.0)
    ax1.plot(time, stress[:,1,1], 'g-', label='yy', linewidth=1.0)
    ax1.plot(time, stress[:,2,2], 'b-', label='zz', linewidth=1.0)
    ax1.plot((time[0],time[-1]), (mean_stress[0,0], mean_stress[0,0]), 'r-',
            label=r'$\langle S_{{xx}} \rangle$ = {0:<10.4f}'.format(mean_stress[0,0]))
    ax1.plot((time[0],time[-1]), (mean_stress[1,1], mean_stress[1,1]), 'g-',
            label=r'$\langle S_{{yy}} \rangle$ = {0:<10.4f}'.format(mean_stress[1,1]))
    ax1.plot((time[0],time[-1]), (mean_stress[2,2], mean_stress[2,2]), 'b-',
            label=r'$\langle S_{{zz}} \rangle$ = {0:<10.4f}'.format(mean_stress[2,2]))

    ax2.plot(data['time'], data['P'])

#    ax2.plot(time, stress[:,0,1], 'r-', label='xy')
#    ax2.plot(time, stress[:,1,0], 'r--', label='yx')
#    ax2.plot(time, stress[:,1,2], 'g-', label='yz')
#    ax2.plot((time[0],time[-1]), (mean_stress[0,1], mean_stress[0,1]), 'r-',
#            label=r'$\langle S_{{xy}} \rangle$ = {0:<10.4f}'.format(mean_stress[0,1]))
#    ax2.plot((time[0],time[-1]), (mean_stress[0,2], mean_stress[0,2]), 'g-',
#            label=r'$\langle S_{{xz}} \rangle$ = {0:<10.4f}'.format(mean_stress[0,2]))
#    ax2.plot((time[0],time[-1]), (mean_stress[1,2], mean_stress[1,2]), 'b-',
#            label=r'$\langle S_{{yz}} \rangle$ = {0:<10.4f}'.format(mean_stress[1,2]))
    ax1.legend(bbox_to_anchor=(1.05,1), loc=2, borderaxespad=0.)
#    ax2.legend(bbox_to_anchor=(1.05,1), loc=2, borderaxespad=0.)
    fig2.subplots_adjust(hspace=0)
    plt.setp([a.get_xticklabels() for a in fig1.axes[:-1]], visible=False)
    fig2.savefig("stress.pdf", bbox_inches='tight')


# Plot the VACF
  if opts.vacf:
    vacf = sp.array(vacf)
    plt.figure("VACF")
    plt.xlabel("t")
    plt.ylabel("C(t)")
    plt.xlim((time[0],time[-1]))
    plt.plot(time, vacf)
    plt.plot((0,time[-1]), (0, 0), 'k-')
    plt.savefig("vacf.pdf", bbox_inches='tight')

# Plot the MSD
  if opts.msd:
    msd = sp.array(msd)
    plt.figure("MSD")
    plt.xlabel("t")
    plt.ylabel("MSD(t)")
    plt.xlim((time[0],time[-1]))
    plt.plot(time, msd)
    plt.plot((0,time[-1]), (0, 0), 'k-')
    plt.ylim(ymin=0)
    plt.savefig("msd.pdf", bbox_inches='tight')

# Plot the velocity distribution
  if opts.vdist:
    vdistr = vdistr/float(nframes)
    sdistr = sdistr/float(nframes)
    plt.figure("v_distribution")
    fig3, ((ax1, ax2), (ax3, ax4)) = plt.subplots(nrows=2, ncols=2)
    ax1.hist(sdistr, bins=bin_edges, label='speed')
    ax2.hist(vdistr[:,0], bins=bin_edges, label=r'$v_x$')
    ax3.hist(vdistr[:,1], bins=bin_edges, label=r'$v_y$')
    ax4.hist(vdistr[:,2], bins=bin_edges, label=r'$v_z$')
    plt.savefig("vdistr.pdf", bbox_inches='tight')
