import os
import re
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter

# Lists to store cutoff energies, total energies, and CPU times
cutoffs = []
energies = []
cpu_times = []

# basename for directory
base = 'kp'

testfolder = False
# Loop over all folders starting with "E_"
for folder in sorted(os.listdir('.')):
    if folder.startswith(base) and os.path.isdir(folder):
        testfolder = True
        filepath = os.path.join(folder, "Conquest_out")
        if os.path.isfile(filepath):
            with open(filepath, 'r') as f:
                lines = f.readlines()
            
            # Find all lines containing "DFT Total Energy"
            dft_lines = [line for line in lines if "DFT Total Energy" in line]
            if dft_lines:
                last_line = dft_lines[-1]
                match = re.search(r":\s*([-+]?\d*\.\d+|\d+)", last_line)
                if match:
                    energy = float(match.group(1))
                else:
                    energy = None
            else:
                energy = None

            # Find CPU time
            time_lines = [line for line in lines if "Total run time was:" in line]
            if time_lines:
                last_time_line = time_lines[-1]
                match_time = re.search(r"Total run time was:\s*([\d.]+)", last_time_line)
                cpu_time = float(match_time.group(1)) if match_time else None
            else:
                cpu_time = None

            # Extract cutoff value from folder name (E_50 -> 50)
            cutoff = float(folder.split("_")[1])

            cutoffs.append(cutoff)
            energies.append(energy)
            cpu_times.append(cpu_time)

            print(f"{folder} : Energy = {energy} Ha, CPU Time = {cpu_time} s")

if (testfolder):
    
    # Sort by cutoff
    cutoffs, energies, cpu_times = zip(*sorted(zip(cutoffs, energies, cpu_times)))

    # Plot Total Energy vs Cutoff
    plt.figure(figsize=(6,4))
    plt.plot(cutoffs, energies, marker='o', linestyle='-')
    plt.xlabel("Grid Cutoff Energy (Ha)")
    plt.ylabel("DFT Total Energy (Ha)")
    plt.title("Convergence of Total Energy vs Grid Cutoff")
    ax = plt.gca()
    ax.yaxis.set_major_formatter(ScalarFormatter(useMathText=False))
    ax.ticklabel_format(style="plain", axis="y", useOffset=False)

    plt.grid(True)
    plt.tight_layout()
    plt.savefig("energy_vs_cutoff.png", dpi=300,  transparent=True)
    plt.show()

    # Plot CPU Time vs Cutoff
    plt.figure(figsize=(6,4))
    plt.plot(cutoffs, cpu_times, marker='s', color='orange', linestyle='-')
    plt.xlabel("Grid Cutoff Energy (Ha)")
    plt.ylabel("CPU Time (s)")
    plt.title("CPU Time vs Grid Cutoff")
    plt.grid(True)
    plt.tight_layout()
    plt.savefig("cpu_time_vs_cutoff.png", dpi=300,  transparent=True)
    plt.show()
