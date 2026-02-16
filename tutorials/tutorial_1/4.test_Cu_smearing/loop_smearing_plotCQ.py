import os
import re
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter

# Lists for FD smearing
sigma_FD = []
energies_FD = []
cpu_FD = []
cycle_FD = [] 

# Lists for MP smearing
sigma_MP = []
energies_MP = []
cpu_MP = []
cycle_MP = [] 

# Loop over all folders
for folder in sorted(os.listdir('.')):
    if os.path.isdir(folder):
        filepath = os.path.join(folder, "Conquest_out")
        if os.path.isfile(filepath):
            with open(filepath, 'r') as f:
                lines = f.readlines()

            # Find all lines containing "DFT Total Energy"
            dft_lines = [line for line in lines if "DFT Total Energy" in line]
            if not dft_lines:
                continue
            last_line = dft_lines[-1]

            cycle = len(dft_lines)
            
            # Extract numeric value after colon
            match = re.search(r":\s*([-+]?\d*\.\d+|\d+)", last_line)
            if not match:
                continue
            energy = float(match.group(1))

            # Extract CPU time
            time_lines = [line for line in lines if "Total run time was:" in line]
            if time_lines:
                last_time_line = time_lines[-1]
                match_time = re.search(r"Total run time was:\s*([\d.]+)", last_time_line)
                cpu_time = float(match_time.group(1)) if match_time else None
            else:
                cpu_time = None

            # Extract sigma value from folder name (FD_0.1 or MP_0.1)
            sigma_match = re.search(r"_(\d+\.?\d*)$", folder)
            if not sigma_match:
                continue
            sigma_val = float(sigma_match.group(1))

            # Sort into the correct smearing type
            if folder.startswith("fd"):
                sigma_FD.append(sigma_val)
                energies_FD.append(energy)
                cpu_FD.append(cpu_time)
                cycle_FD.append(cycle)
                
            elif folder.startswith("mp"):
                sigma_MP.append(sigma_val)
                energies_MP.append(energy)
                cpu_MP.append(cpu_time)
                cycle_MP.append(cycle)

# Sort datasets by sigma
if sigma_FD:
    sigma_FD, energies_FD, cpu_FD = zip(*sorted(zip(sigma_FD, energies_FD, cpu_FD)))
if sigma_MP:
    sigma_MP, energies_MP, cpu_MP = zip(*sorted(zip(sigma_MP, energies_MP, cpu_MP)))

# Plot Total Energy vs Smearing Width
plt.figure(figsize=(7,5))
if sigma_FD:
    plt.plot(sigma_FD, energies_FD, marker='o', color='blue', label='FD')
if sigma_MP:
    plt.plot(sigma_MP, energies_MP, marker='s', color='red', label='MP')

plt.xlabel("Smearing width (Ha)")
plt.ylabel("DFT Total Energy (Ha)")
plt.title("Total Energy vs Smearing Width")

ax = plt.gca()
ax.yaxis.set_major_formatter(ScalarFormatter(useMathText=False))
ax.ticklabel_format(style="plain", axis="y", useOffset=False)
ax.set_xscale('log')

plt.grid(True)
plt.legend()
plt.tight_layout()
plt.savefig("energy_vs_sigma_FD_MP.png", dpi=300)
plt.show()

# Plot CPU Time vs Smearing Width
plt.figure(figsize=(7,5))
if sigma_FD:
    plt.plot(sigma_FD, cpu_FD, marker='o', color='blue', label='FD')
    for i, txt in enumerate(cycle_FD) : 
        plt.annotate('  '+str(txt)+' SCF cycles', (sigma_FD[i], cpu_FD[i]), color='blue')
        
if sigma_MP:
    plt.plot(sigma_MP, cpu_MP, marker='s', color='red', label='MP')
    for i, txt in enumerate(cycle_MP) : 
        plt.annotate('  '+str(txt)+' SCF cycles', (sigma_MP[i], cpu_MP[i]), color='red')

plt.xlabel("Smearing width (Ha)")
plt.ylabel("CPU Time (s)")
plt.title("CPU Time vs Smearing Width")

ax = plt.gca()
ax.yaxis.set_major_formatter(ScalarFormatter(useMathText=False))
ax.ticklabel_format(style="plain", axis="y", useOffset=False)

plt.grid(True)
plt.legend()
plt.tight_layout()
plt.savefig("cpu_vs_sigma_FD_MP.png", dpi=300)
plt.show()

