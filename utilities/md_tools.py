#!/usr/local/bin/python3

import scipy as sp
import matplotlib.pyplot as plt
from scipy.linalg import norm
from scipy.integrate import cumtrapz
from scipy.signal import correlate
from scipy import histogram
from math import ceil, pi
from frame import Frame
from pdb import set_trace

bohr2ang = 0.529177249
small = 1.0e-3

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

def disp_mic_npt(pos1, pos2, cell1, cell2):
  """MIC displacement when cell dimensions change"""
  disp = pos2/cell2 - pos1/cell1
  for i in range(3):
    disp[i] -= round(disp[i]/cell2[i])
  return disp

class Pairdist:
  """Object for computing pair distribution functions"""

  def __init__(self, nat, nspec, rcut, binwidth, species, species_count):
    self.nframes = 0
    self.nat = nat
    self.nspec = nspec
    self.rcut = rcut
    self.binwidth = binwidth
    self.nbins = ceil(rcut/binwidth)+1
    self.spec_count = species_count
    self.species = species
    self.bins = []
    for i in range(self.nbins):
      self.bins.append((float(i)*binwidth + binwidth/2.))
    self.bins = sp.array(self.bins)
    self.dt = sp.zeros((self.nat,self.nat), dtype='float')
    self.freq_total = sp.zeros(self.nbins, dtype='int')
    self.freq = sp.zeros((self.nbins,self.nspec,self.nspec), dtype='int')
    self.nfac_total = sp.zeros(self.nbins, dtype='float')
    self.nfac = sp.zeros((self.nbins,self.nspec,self.nspec), dtype='float')
    self.gr_total = sp.zeros(self.nbins, dtype='float')
    self.gr = sp.zeros((self.nbins,self.nspec,self.nspec), dtype='float')

  def update_rdf(self, frame):
    self.nframes += 1
    cell = sp.zeros(3)
    for i in range(3):
      cell[i] = frame.lat[i,i]

    self.volume = cell[0]*cell[1]*cell[2]*bohr2ang**3
    self.rho = float(self.nat)/self.volume

    for i in range(self.nat):
      for j in range(i+1, self.nat):
        diff = diff_mic(frame.r[i,:], frame.r[j,:], cell)*bohr2ang
        self.dt[i,j] = norm(diff)
        self.dt[j,i] = norm(diff)
        if self.dt[i,j] < self.rcut:
          ind = int(round((self.dt[i,j]+self.binwidth)/self.binwidth))-1
          self.freq_total[ind] += 2
          if self.nspec > 1:
            for ispec in range(self.nspec):
              for jspec in range(ispec, self.nspec):
                if (ispec == frame.species[i]-1 and jspec == frame.species[j]-1):
                  self.freq[ind, ispec, jspec] += 2

  def norm_rdf(self):
    """Normalise the RDF"""
    const1 = 4.0*pi*(self.binwidth**3)/3.0
    const2 = self.rho*self.nat*self.nframes
    for i in range(self.nbins):
      vshell = (float(i+1)**3 - float(i)**3)*const1
      self.nfac_total[i] = vshell*const2
      if self.nspec > 1:
        for ispec in range(self.nspec):
          for jspec in range(self.nspec):
            const3 = self.rho*self.spec_count[ispec+1]*self.spec_count[jspec+1]/self.nat
            self.nfac[i,ispec,jspec] = vshell*const3*self.nframes
    self.gr_total = self.freq_total.astype(float)/self.nfac_total
    if self.nspec > 1:
      self.gr = self.freq.astype(float)/self.nfac

  def get_coordination(self):
    """Compute coordination"""
    gxrsq = self.gr_total*self.bins**2
    self.coord_total = sp.zeros(self.nbins, dtype='float')
    self.coord_total[1:] = cumtrapz(gxrsq,self.bins)
    self.coord_total *= 4.*pi*self.rho
    if self.nspec > 1:
      self.coord = sp.zeros((self.nbins,self.nspec,self.nspec), dtype='float')
      for ispec in range(self.nspec):
        for jspec in range(ispec,self.nspec):
          gxrsq = self.gr[:,ispec,jspec]*self.bins**2
          self.coord[1:,ispec,jspec] = cumtrapz(gxrsq[:], self.bins)
          self.coord *= 4.*pi*self.rho # check this

  def plot_gr(self):
    plt.figure("RDF")
    filename = "rdf.pdf"
    if (self.nspec > 1):
      fig3, (axl, ax2) = plt.subplots(nrows=2, ncols=1, sharex=True)
    else:
      fig3, axl = plt.subplots()
    axl.minorticks_on()
    axl.grid(b=True, which='major', axis='x', color='gray', linestyle='-')
    axl.grid(b=True, which='minor', axis='x', color='gray', linestyle='--')
    axl.grid(b=True, which='major', axis='y', color='gray', linestyle='-')
    # axl.grid(b=True, which='minor', axis='y', color='gray', linestyle='--')
    axr = axl.twinx()
    axl.set_ylabel("g(r)", color='b')
    axr.set_ylabel("Coordination", color='r')
    axl.plot(self.bins, self.gr_total,'b-', label="total", linewidth=1.0)
    axr.plot(self.bins, self.coord_total, 'r-', label="total", linewidth=1.0)
    plt.xlim((0,self.rcut))
    axl.set_ylim(bottom=0)
    axr.set_ylim(bottom=axl.get_ylim()[0], top=axl.get_ylim()[1]*10.0)
    if self.nspec > 1:
      ax2.minorticks_on()
      ax2.grid(b=True, which='major', axis='x', color='gray', linestyle='-')
      ax2.grid(b=True, which='minor', axis='x', color='gray', linestyle='--')
      ax2.grid(b=True, which='major', axis='y', color='gray', linestyle='-')
      # ax2.grid(b=True, which='minor', axis='y', color='gray', linestyle='--')
      for ispec in range(self.nspec):
        for jspec in range(ispec,self.nspec):
          pair = "{}-{}".format(self.species[ispec+1], self.species[jspec+1])
          ax2.plot(self.bins, self.gr[:,ispec,jspec], label=pair, linewidth=1.0)
      ax2.set_ylim(bottom=0)
      ax2.set_xlabel("r (A)")
      ax2.set_ylabel("partial g(r)")
      ax2.legend(loc="upper right")
    else:
      axl.set_xlabel("r (A)")
    fig3.savefig(filename, bbox_inches='tight')

  def dump_gr(self):
    header_bit = "{0:>16s}"
    rdf_bit = "{0:>16.6f}"
    header_fmt = "{0:>16s}{1:>16s}{2:>16s}"
    rdf_fmt = "{0:>16.6f}{1:>16.6f}{2:>16.6f}"
    filename = "rdf.dat"

    header = header_fmt.format("r (A)", "total", "coordination")
    if self.nspec > 1:
      for ispec in range(self.nspec):
        for jspec in range(ispec,self.nspec):
          pair = "{}-{}".format(self.species[ispec+1], self.species[jspec+1])
          header += header_bit.format(pair)
    header += "\n"

    with open(filename, 'w') as outfile:
      outfile.write(header)
      for i in range(self.nbins):
        rdf_line = rdf_fmt.format(self.bins[i], self.gr_total[i],
                                  self.coord_total[i])
        if self.nspec > 1:
          for ispec in range(self.nspec):
            for jspec in range(ispec,self.nspec):
              rdf_line += rdf_bit.format(self.gr[i,ispec,jspec])
        rdf_line += "\n"
        outfile.write(rdf_line)


  def get_bondlength(self, bondcut, frame, printall):

    bond_tot = sp.zeros((self.nspec, self.nspec), dtype=float)
    bondsq_tot = sp.zeros((self.nspec, self.nspec), dtype=float)
    bond_avg = sp.zeros((self.nspec, self.nspec), dtype=float)
    bond_sd = sp.zeros((self.nspec, self.nspec), dtype=float)
    bond_min = sp.zeros((self.nspec, self.nspec), dtype=float)
    nbonds = sp.zeros((self.nspec, self.nspec), dtype=int)

    bond_min = bondcut
    for i in range(self.nat):
      for j in range(i+1, self.nat):
        s1 = frame.species[i]-1
        s2 = frame.species[j]-1
        if self.dt[i,j] < bondcut[s1,s2]:
          if self.dt[i,j] > small:
            bond_tot[s1,s2] += self.dt[i,j]
            bondsq_tot[s1,s2] += self.dt[i,j]**2
            nbonds[s1,s2] += 1
            if self.dt[i,j] < bond_min[s1,s2]:
              bond_min[s1,s2] = self.dt[i,j]

            if printall:
              pair = "{}--{}".format(self.species[s1+1], self.species[s2+1])
              print(f'{pair}: {i:>4d}--{j:<4d} {self.dt[i,j]:>8.4f}')

    print("Mean bond lengths:")
    for i in range(self.nspec):
      for j in range(i,self.nspec):
        if nbonds[i,j] > 0:
          bond_avg[i,j] = bond_tot[i,j]/float(nbonds[i,j])
          bond_sd[i,j] = sp.sqrt(bondsq_tot[i,j]/nbonds[i,j] - bond_avg[i,j]**2)
          pair = "{}-{}".format(self.species[i+1], self.species[j+1])
          print(f'{pair}: {bond_avg[i,j]:>8.4f} +/- {bond_sd[i,j]:>8.4f}')

    print("Minimum bond lengths:")
    for i in range(self.nspec):
      for j in range(i,self.nspec):
          pair = "{}-{}".format(self.species[i+1], self.species[j+1])
          print(f'{pair}: {bond_min[i,j]:>8.4f}')

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

class VACF:
  """Velocity autocorrelation function"""

  def __init__(self, nat, dt, init_frame):
    self.nframes = 0
    self.nat = nat
    self.dt = dt
    self.init_v = init_frame.v
    self.vacf = []
    self.steps = []

  def update_vacf(self, step, frame):
    self.nframes += 1
    self.steps.append(step)
    self.vacf.append(0.0)
    for i in range(self.nat):
      self.vacf[-1] += sp.dot(self.init_v[i,:], frame.v[i,:])

  def norm_vacf(self):
    self.vacf = sp.array(self.vacf)/self.nat
    self.time = sp.array(self.steps, dtype='float')*self.dt

  def plot_vacf(self):
    filename = "vacf.pdf"
    plt.figure("VACF")
    plt.xlabel("t (fs)")
    plt.ylabel("C(t)")
    plt.xlim((self.time[0],self.time[-1]))
    plt.plot(self.time, self.vacf)
    plt.plot((0,self.time[-1]), (0, 0), 'k-')
    plt.savefig(filename, bbox_inches='tight')

  def dump_vacf(self):
    vacf_fmt = "{0:>12.4f}{1:>16.6f}\n"
    filename = "vacf.dat"
    with open(filename, 'w') as outfile:
      for i in range(self.nframes):
        outfile.write(vacf_fmt.format(self.time[i], self.vacf[i]))

class MSD:
  """Mean Squared Deviation"""

  def __init__(self, nat, dt, init_frame):
    self.nframes = 0
    self.nat = nat
    self.dt = dt
    self.init_r = init_frame.r
    self.r_prev = sp.copy(self.init_r)
    self.msd = []
    self.steps = []
    self.init_cell = sp.zeros(3, dtype='float')
    self.r_diff = sp.zeros((self.nat,3), dtype='float')
    for i in range(3):
      self.init_cell[i] = init_frame.lat[i,i]

  def update_msd(self, step, frame):
    self.steps.append(step)
    self.nframes += 1
    cell = sp.zeros(3, dtype='float')
    for i in range(3):
      cell[i] = frame.lat[i,i]
    self.msd.append(0.0)
    for i in range(self.nat):
      diff = diff_mic(frame.r[i,:], self.r_prev[i,:], cell) # is this right?
      self.r_diff[i,:] += diff
      self.msd[-1] += sp.sum(self.r_diff[i,:]**2)
    self.r_prev = sp.copy(frame.r)

  def norm_msd(self):
    self.msd = sp.array(self.msd)/self.nat
    self.time = sp.array(self.steps, dtype='float')*self.dt

  def plot_msd(self):
    filename = "msd.pdf"
    plt.figure("MSD")
    plt.xlabel("t (fs)")
    plt.ylabel("MSD")
    plt.xlim((self.time[0],self.time[-1]))
    plt.plot(self.time, self.msd)
    plt.plot((0,self.time[-1]), (0, 0), 'k-')
    plt.ylim(ymin=0)
    plt.savefig(filename, bbox_inches='tight')

  def dump_msd(self):
    msd_fmt = "{0:>12.4f}{1:>16.6f}\n"
    filename = "msd.dat"
    with open(filename, 'w') as outfile:
      for i in range(self.nframes):
        outfile.write(msd_fmt.format(self.time[i], self.msd[i]))
