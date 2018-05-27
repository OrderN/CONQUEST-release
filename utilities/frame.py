import re
import scipy as sp
from pdb import set_trace

# Regular expressions
cell_re = re.compile('cell_vectors(.*?)end cell_vectors', re.M | re.S)
stress_re = re.compile('stress_tensor(.*?)end stress_tensor', re.M | re.S)
position_re = re.compile('positions(.*?)end positions', re.M | re.S)
position_re = re.compile('positions(.*?)end positions', re.M | re.S)
velocity_re = re.compile('velocities(.*?)end velocities', re.M | re.S)
force_re = re.compile('forces(.*?)end forces', re.M | re.S)

class Frame:
  """Stores a frame from a MD trajectory"""

  def __init__(self, nat, step):
    self.step = step
    self.nat = nat
    self.species = sp.zeros(nat)
    self.r = sp.zeros((nat, 3), dtype='float')
    self.v = sp.zeros((nat, 3), dtype='float')
    self.f = sp.zeros((nat, 3), dtype='float')
    self.lat = sp.zeros((3, 3), dtype='float')
    self.stress = sp.zeros((3, 3), dtype='float')
    self.ke = 0.
    self.pe = 0.
    self.E = 0.
    self.T = 0.
    self.P = 0.
    self.vmax = 1.0

  def parse_frame(self, buf):
    """Read frame data from a string buffer"""
    m = re.search(cell_re, buf)
    lines = m.group(1).strip().splitlines()
    for i in range(3):
      bits = lines[i].strip().split()
      for j in range(3):
        self.lat[i,j] = float(bits[j])

    m = re.search(stress_re, buf)
    if m:
      lines = m.group(1).strip().splitlines()
      for i in range(3):
        bits = lines[i].strip().split()
        for j in range(3):
          self.stress[i,j] = float(bits[j])

    m = re.search(position_re, buf)
    lines = m.group(1).strip().splitlines()
    nat = len(lines)
    for i in range(nat):
      bits = lines[i].strip().split()
      bits.pop(0)
      self.species[i] = int(bits.pop(0))
      for j in range(3):
        self.r[i,j] = float(bits[j])

    m = re.search(velocity_re, buf)
    lines = m.group(1).strip().splitlines()
    nat = len(lines)
    for i in range(nat):
      bits = lines[i].strip().split()
      bits.pop(0)
      bits.pop(0)
      for j in range(3):
        self.v[i,j] = float(bits[j])

    m = re.search(force_re, buf)
    lines = m.group(1).strip().splitlines()
    nat = len(lines)
    for i in range(nat):
      bits = lines[i].strip().split()
      bits.pop(0)
      bits.pop(0)
      for j in range(3):
        self.f[i,j] = float(bits[j])