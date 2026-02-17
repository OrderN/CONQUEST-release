import os
import re
import matplotlib.pyplot as plt

# Lists to store energies, total energies, and CPU times
kpoints_labels = []
energies = []
cpu_times = []

# basename for directory
base = 'kp'

testfolder = False
# Loop over all folders
for folder in sorted(os.listdir('.')):
    if folder.startswith(base) and os.path.isdir(folder):
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

            # Extract k-points label
            kpts_label = folder.split("_")[1]

            kpoints_labels.append(kpts_label)
            energies.append(energy)
            cpu_times.append(cpu_time)
            print(f"{folder} : Energy = {energy} Ha, CPU Time = {cpu_time} s")

if (testfolder):
    # Convert labels to numeric value for sorting (first number)
    k_numeric = [int(label.split("x")[0]) for label in kpoints_labels]

    # Sort by numeric value
    sorted_data = sorted(zip(k_numeric, kpoints_labels, energies, cpu_times), key=lambda x: x[0])
    kpoints_labels_sorted = [x[1] for x in sorted_data]
    energies_sorted = [x[2] for x in sorted_data]
    cpu_times_sorted = [x[3] for x in sorted_data]

    # Plot Total Energy vs K-points
    x_indices = range(len(energies_sorted))
    plt.figure(figsize=(6,4))
    plt.plot(x_indices, energies_sorted, marker='o', linestyle='-')
    plt.xticks(x_indices, kpoints_labels_sorted)
    plt.xlabel("K-points grid")
    plt.ylabel("DFT Total Energy (Ha)")
    plt.title("Total Energy vs K-points")
    plt.grid(True)
    plt.tight_layout()
    plt.savefig("energy_vs_kpoints.png", dpi=300)
    plt.show()

# Plot CPU Time vs K-points
    plt.figure(figsize=(6,4))
    plt.plot(x_indices, cpu_times_sorted, marker='s', color='orange', linestyle='-')
    plt.xticks(x_indices, kpoints_labels_sorted)
    plt.xlabel("K-points grid")
    plt.ylabel("CPU Time (s)")
    plt.title("CPU Time vs K-points")
    plt.grid(True)
    plt.tight_layout()
    plt.savefig("cpu_time_vs_kpoints.png", dpi=300)
    plt.show()
