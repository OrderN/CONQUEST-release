#!/usr/local/bin/python3

import scipy as sp
import matplotlib.pyplot as plt
from scipy.linalg import norm
from scipy.integrate import cumtrapz
from scipy import histogram
from math import ceil, pi
from frame import Frame
from pdb import set_trace

bohr2ang = 0.529177249

def diff_mic(pos1, pos2, cell):
  """Minimum image convenciton relative vector (orthorhombic ell only)"""
  diff = pos2 - pos1
  for i in range(3):
    diff[i] -= round(diff[i]/cell[i])
  return diff

def disp_mic_npt(pos1, pos2, cell1, cell2):
  """MIC displacement when cell dimensions change"""
  disp = pos2/cell2 - pos1/cell1
  for i in range(3):
    disp[i] -= round(disp[i]/cell2[i])
  return disp

class Pairdist:
  """Object for computing pair distribution functions"""

  def __init__(self, nat, nspec, rcut, binwidth, frame):
    self.nframes = 0
    self.nat = nat
    self.nspec = nspec
    self.rcut = rcut
    self.binwidth = binwidth
    self.nbins = ceil(rcut/binwidth)+1
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
    self.spec_count = sp.zeros(self.nspec, dtype='int')
    for i in range(self.nspec):
      for j in range(self.nat):
        if (frame.species[j] == i+1):
          self.spec_count[i] += 1

  def get_dt(self, frame):
    """Generate the distance table"""
    cell = sp.zeros(3)
    for i in range(3):
      cell[i] = frame.lat[i,i]

    self.volume = cell[0]*cell[1]*cell[2]*bohr2ang**3
    self.rho = float(self.nat)/self.volume

    for i in range(self.nat):
      for j in range(i+1, self.nat):
        diff = diff_mic(frame.r[i,:], frame.r[j,:], cell)*bohr2ang
        self.dt[i,j] = norm(diff)
        self.dt[j,i] = self.dt[i,j]

  def update_rdf_dt(self, frame):
    """update the radial distribution function"""
    self.get_dt(frame)
    for i in range(self.nat):
      for j in range(i+1,self.nat):
        if (self.dt[i,j] < self.rcut):
          ind = int(round((self.dt[i,j]+self.binwidth)/self.binwidth))-1
          self.freq_total[ind] += 2
          if self.nspec > 1:
            for ispec in range(self.nspec):
              for jspec in range(self.nspec):
                if (ispec == frame.species[i] and jspec == frame.species[j]):
                  self.freq[ind, ispec, jspec] += 1

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
        d = norm(diff)
        if d < self.rcut:
          ind = int(round((d+self.binwidth)/self.binwidth))-1
          self.freq_total[ind] += 2
          if self.nspec > 1:
            for ispec in range(self.nspec):
              for jspec in range(ispec, self.nspec):
                if (ispec == frame.species[i] and jspec == frame.species[j]):
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
            const3 = self.rho*self.spec_count[ispec]*self.spec_count[jspec]/self.nat
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

  def plot_gr(self):
    plt.figure("RDF")
    filename = "rdf.pdf"
    fig3, axl = plt.subplots()
    plt.grid(b=True, which='major', axis='both', color='gray', linestyle='-')
    plt.grid(b=True, which='minor', axis='both', color='gray', linestyle='--')
    axr = axl.twinx()
    axl.set_xlabel("r (A)")
    axl.set_ylabel("g(r)")
    axr.set_ylabel("Coordination")
    axl.plot(self.bins, self.gr_total,'b-', label="total", linewidth=1.0)
    axr.plot(self.bins, self.coord_total, 'r-', label="total", linewidth=1.0)
    plt.xlim((0,self.rcut))
    plt.minorticks_on()
    axl.set_ylim(bottom=0)
    axr.set_ylim(bottom=0)
    fig3.savefig(filename, bbox_inches='tight')

  def dump_gr(self):
    rdf_fmt = "{0:>16.6f}{1:>16.6f}{2:>16.6f}\n"
    filename = "rdf.dat"
    with open(filename, 'w') as outfile:
      for i in range(self.nbins):
        outfile.write(rdf_fmt.format(self.bins[i], self.gr_total[i],
                      self.coord_total[i]))

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
    mser_min = traj.mser_min()
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