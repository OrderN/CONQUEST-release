import os
import re
import matplotlib.pyplot as plt

basis_order = ["bs_SZ", "bs_SZP", "bs_DZP", "bs_TZTP"]

basis_labels = []
energies = []
cpu_times = []

testfolder = False
# Loop over the folders in the specified order
for folder in basis_order:
    if os.path.isdir(folder):
        testfolder = True
        filepath = os.path.join(folder, "Conquest_out")
        if os.path.isfile(filepath):
            with open(filepath, 'r') as f:
                lines = f.readlines()

            # Extract DFT total energy
            dft_lines = [line for line in lines if "DFT Total Energy" in line]
            if dft_lines:
                last_line = dft_lines[-1]
                match = re.search(r":\s*([-+]?\d*\.\d+|\d+)", last_line)
                energy = float(match.group(1)) if match else None
            else:
                energy = None

            # Extract CPU time
            time_lines = [line for line in lines if "Total run time was:" in line]
            if time_lines:
                last_time_line = time_lines[-1]
                match_time = re.search(r"Total run time was:\s*([\d.]+)", last_time_line)
                cpu_time = float(match_time.group(1)) if match_time else None
            else:
                cpu_time = None

            basis_labels.append(folder)
            energies.append(energy)
            cpu_times.append(cpu_time)
            print(f"{folder} : Energy = {energy} Ha, CPU Time = {cpu_time} s")
        else:
            print(f"{folder} : Conquest_out not found")
    else:
        print(f"{folder} : Folder not found")

if (testfolder):
    # Plot Total Energy vs Basis Set
    x_indices = range(len(energies))
    plt.figure(figsize=(6,4))
    plt.plot(x_indices, energies, marker='o', linestyle='-', label="Energy")
    plt.xticks(x_indices, basis_labels)
    plt.xlabel("Basis Set")
    plt.ylabel("DFT Total Energy (Ha)")
    plt.title("Total Energy vs Basis Set")
    plt.grid(True)
    plt.tight_layout()
    plt.savefig("energy_vs_basis.png", dpi=300)
    plt.show()
    
    # Plot CPU time vs Basis Set
    plt.figure(figsize=(6,4))
    plt.plot(x_indices, cpu_times, marker='s', color='orange', linestyle='-', label="CPU Time")
    plt.xticks(x_indices, basis_labels)
    plt.xlabel("Basis Set")
    plt.ylabel("CPU Time (s)")
    plt.title("CPU Time vs Basis Set")
    plt.grid(True)
    plt.tight_layout()
    plt.savefig("cpu_time_vs_basis.png", dpi=300)
    plt.show()
