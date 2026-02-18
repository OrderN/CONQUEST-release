import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
from matplotlib.ticker import LogFormatter

# Conversion Bohr -> Å
bohr2ang = 0.529177

# Lists to store results
force_tol_list = []
energy_list = []
rCH_list = []

# Loop over folders starting with Force_
for folder in sorted(os.listdir('.')):
    if not folder.startswith('f_') or not os.path.isdir(folder):
        continue

    out_file = os.path.join(folder, "Conquest_out")
    coord_file = os.path.join(folder, "coord_next.dat")

    if not os.path.exists(out_file) or not os.path.exists(coord_file):
        print(f"Skipping {folder}: missing files")
        continue

    # Extract force tolerance from folder name
    try:
        force_tol = float(folder.split("_")[1])
    except:
        print(f"Cannot parse force tolerance from {folder}")
        continue

    # Extract last GeomOpt energy
    with open(out_file) as f:
        geom_lines = [line for line in f if "GeomOpt - Iter:" in line]

    if not geom_lines:
        print(f"No GeomOpt in {folder}")
        continue

    last_geom = geom_lines[-1]
    # Find position of E: and take next number
    pos = last_geom.find("E:")
    if pos == -1:
        print(f"Cannot find E: in {folder}")
        continue
    rest = last_geom[pos+2:].strip()
    energy = float(rest.split()[0])
    energy_list.append(energy)

    # Read coord_next.dat and compute C-H distance 
    with open(coord_file) as f:
        lines = f.readlines()

    # Lattice vectors (Bohr)
    lattice = np.array([[float(x) for x in lines[i].split()] for i in range(3)])

    # Number of atoms
    natoms = int(lines[3])

    # Coordinates and labels
    coords_frac = []
    labels = []
    for line in lines[4:4+natoms]:
        parts = line.split()
        coords_frac.append([float(parts[0]), float(parts[1]), float(parts[2])])
        labels.append(int(parts[3]))
    coords_frac = np.array(coords_frac)
    labels = np.array(labels)

    # Convert fractional to Cartesian
    coords_cart = coords_frac @ lattice

    # Compute average C-H distance
    C_pos = coords_cart[labels==2][0]  # Carbon
    H_pos = coords_cart[labels==1]     # Hydrogens
    distances = np.linalg.norm(H_pos - C_pos, axis=1)
    rCH = np.mean(distances) * bohr2ang
    rCH_list.append(rCH)

    force_tol_list.append(force_tol)

# Sort results by force tolerance
force_tol_list = np.array(force_tol_list)
energy_list = np.array(energy_list)
rCH_list = np.array(rCH_list)

sorted_idx = np.argsort(force_tol_list)
force_tol_list = force_tol_list[sorted_idx]
energy_list = energy_list[sorted_idx]
rCH_list = rCH_list[sorted_idx]

# Plot energy and C-H distance
fig, ax1 = plt.subplots(figsize=(8,6))

color1 = 'tab:blue'
ax1.set_xlabel("Force Tolerance [Ha/Bohr]",size=12)
ax1.set_ylabel("Total Energy [Ha]", size=12, color=color1)
ax1.plot(force_tol_list, energy_list, "o-", color=color1, label="Energy")
ax1.tick_params(axis='y', labelcolor=color1)
ax1.set_xscale('log')
ax1.xaxis.set_major_formatter(LogFormatter(base=10, labelOnlyBase=True))
ax = plt.gca()
ax.yaxis.set_major_formatter(ScalarFormatter(useMathText=False))
ax.ticklabel_format(style="plain", axis="y", useOffset=False)
ax1.grid(True)

ax2 = ax1.twinx()
color2 = 'tab:red'
ax2.set_ylabel("Average C-H Distance [Å]", color=color2,size=12)
ax2.plot(force_tol_list, rCH_list, "s--", color=color2, label="C-H distance")
ax2.tick_params(axis='y', labelcolor=color2)

lines_1, labels_1 = ax1.get_legend_handles_labels()
lines_2, labels_2 = ax2.get_legend_handles_labels()
ax1.legend(lines_1 + lines_2, labels_1 + labels_2, loc="best")

plt.title("Energy and C-H Distance vs Force Tolerance")
ax.yaxis.set_major_formatter(ScalarFormatter(useMathText=False))
ax.ticklabel_format(style="plain", axis="y", useOffset=False)
plt.tight_layout()
plt.savefig("energy_rCH_vs_force_tol.png", dpi=300)
plt.show()

