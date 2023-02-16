import re
import scipy as sp
from pdb import set_trace

# Regular expressions
cell_re = re.compile('cell_vectors(.*?)end cell_vectors', re.M | re.S)
stress_re = re.compile('stress_tensor(.*?)end stress_tensor', re.M | re.S)
position_re = re.compile('positions(.*?)end positions', re.M | re.S)
position_re = re.compile('positions(.*?)end positions', re.M | re.S)
# position_re is declared twice.
velocity_re = re.compile('velocities(.*?)end velocities', re.M | re.S)
force_re = re.compile('forces(.*?)end forces', re.M | re.S)

class Frame:
  """Stores a frame from a MD trajectory"""

  def __init__(self, nat, step): 
    #######################################################
    # Initialisation of the class.
    # In the class Frame(nat,step). nat and step are taken.
    # step is the number of the step.
    # nat will define the number of species involved.
    #######################################################
    self.step = step
    self.nat = nat
    self.species = sp.zeros(nat)
    self.r = sp.zeros((nat, 3), dtype='float')    # For the position
    self.v = sp.zeros((nat, 3), dtype='float')    # For the velocity
    self.f = sp.zeros((nat, 3), dtype='float')    # For the force
    self.lat = sp.zeros((3, 3), dtype='float')    # For the cell tensor
    self.stress = sp.zeros((3, 3), dtype='float') # For the stress tensor
    self.ke = 0.
    self.pe = 0.
    self.E = 0.
    self.T = 0.
    self.P = 0.
    self.vmax = 1.0

  def parse_frame(self, buf):
    """Read frame data from a string buffer"""
    # We read only the data about the cell vectors.
    m = re.search(cell_re, buf)
    lines = m.group(1).strip().splitlines()
    for i in range(3):
      bits = lines[i].strip().split()
      for j in range(3):
        self.lat[i,j] = float(bits[j])

    # We read only the data about the stress tensor.
    m = re.search(stress_re, buf)
    if m:
      lines = m.group(1).strip().splitlines()
      for i in range(3):
        bits = lines[i].strip().split()
        for j in range(3):
          self.stress[i,j] = float(bits[j])

    # We read only the data about the positions.
    m = re.search(position_re, buf)
    lines = m.group(1).strip().splitlines()
    nat = len(lines)
    for i in range(nat):
      bits = lines[i].strip().split()
      bits.pop(0)
      self.species[i] = int(bits.pop(0))
      for j in range(3):
        self.r[i,j] = float(bits[j])

    # We read only the data about the velocity.
    m = re.search(velocity_re, buf)
    lines = m.group(1).strip().splitlines()
    nat = len(lines)
    for i in range(nat):
      bits = lines[i].strip().split()
      bits.pop(0)
      bits.pop(0)
      for j in range(3):
        self.v[i,j] = float(bits[j])

    # We read only the data about the forces.
    m = re.search(force_re, buf)
    lines = m.group(1).strip().splitlines()
    nat = len(lines)
    for i in range(nat):
      bits = lines[i].strip().split()
      bits.pop(0)
      bits.pop(0)
      for j in range(3):
        self.f[i,j] = float(bits[j])