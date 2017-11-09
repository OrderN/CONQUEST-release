import scipy as sp
from pdb import set_trace


class Frame:
  """Stores a frame from a MD trajectory"""

  def __init__(self, nat, step):
    self.step = step
    self.nat = nat
    self.species = sp.zeros(nat)
    self.r = sp.zeros((nat, 3))
    self.v = sp.zeros((nat, 3))
    self.f = sp.zeros((nat, 3))
    self.lat = sp.zeros((3, 3))
    self.stress = sp.zeros((3, 3))
    self.ke = 0.
    self.pe = 0.
    self.E = 0.
    self.T = 0.
    self.P = 0.
    self.vmax = 1.0

  def update_vacf(self, init):
    """Compute the velocity autocorrelation function given the initial frame"""
    return sp.sum(init.v * self.v)/self.nat

  def diff_mic(self, v1, v2):
    diff = v2 - v1
    for i in range(3):
      diff[i] = diff[i] - round(diff[i]/self.lat[i,i])
    return diff/self.nat

  def update_msd(self, init):
    """Compute the mean squared displacement given the initial frame"""
    diff = sp.zeros((self.nat,3))
    for i in range(self.nat):
      diff[i,:] = self.diff_mic(self.r[i,:], init.r[i,:])
    return sp.sum(diff**2)/self.nat

  def update_vdistr(self, nbins):
    """Generate the velocity distribution histogram"""
    vdistr = sp.zeros((nbins,3), dtype='float')
    for i in range(3):
      hist, bin_edges = sp.histogram(self.v[:,i], nbins, range=(0,self.vmax))
      vdistr[:,i] = hist
    return vdistr, bin_edges

  def update_sdistr(self, nbins):
    """Generate the speed distribution histogram"""
    speed = sp.zeros(self.nat, dtype='float')
    for i in range(self.nat):
      speed[i] = sp.sum(self.v[i,:]**2)
    return sp.histogram(speed, nbins, range=(0,self.vmax))
