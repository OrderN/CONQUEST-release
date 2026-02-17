import os
import re
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# Birch-Murnaghan Equation of State (3rd order)
def birch_murnaghan(V, E0, V0, B0, B0p):
    """
    V   : volume (Å³)
    E0  : minimum energy (Ha)
    V0  : equilibrium volume (Å³)
    B0  : bulk modulus (Ha/Å³)
    B0p : derivative of bulk modulus
    """
    eta = (V0 / V)**(2/3)
    return E0 + (9*B0*V0/16.0)*(
        (eta - 1)**3 * B0p +
        (eta - 1)**2 * (6 - 4*eta)
    )

# --- Data extraction ---
a_values = []
energies = []

for folder in sorted(os.listdir('.')):
    if folder.startswith("a_") and os.path.isdir(folder):
        filepath = os.path.join(folder, "Conquest_out")
        if os.path.isfile(filepath):
            with open(filepath, 'r') as f:
                lines = f.readlines()

            # Extract total energy (last occurrence)
            dft_lines = [line for line in lines if "DFT Total Energy" in line]
            if dft_lines:
                last_line = dft_lines[-1]
                match = re.search(r":\s*([-+]?\d*\.\d+|\d+)", last_line)
                energy = float(match.group(1)) if match else None
            else:
                energy = None

            # Extract lattice constant from folder name (e.g. a_5.44 -> 5.44 Å)
            try:
                a_val = float(folder.split("_")[1])
            except ValueError:
                continue

            a_values.append(a_val)
            energies.append(energy)
            print(f"{folder} : a = {a_val:.3f} Å, E = {energy:.6f} Ha")

a_values = np.array(a_values)
energies = np.array(energies)

# Sort by lattice constant
sort_idx = np.argsort(a_values)
a_values = a_values[sort_idx]
energies = energies[sort_idx]

# Compute volume for cubic cell (V = a³)
V_values = a_values**3

# Fit Birch-Murnaghan EOS 
E0_guess = min(energies)
V0_guess = V_values[np.argmin(energies)]
B0_guess = 0.5   # Ha/Å³ (initial guess)
B0p_guess = 4.0

popt, pcov = curve_fit(birch_murnaghan, V_values, energies,
                       p0=[E0_guess, V0_guess, B0_guess, B0p_guess])

E0, V0, B0, B0p = popt

# Convert B0 to GPa
# 1 Ha = 27.2114 eV ; 1 Å³ = 1e-30 m³ ; 1 eV/Å³ = 160.21766208 GPa
B0_GPa = B0 * 27.2114 * 160.21766208

print("\n--- Birch-Murnaghan Fit Results ---")
print(f"E0   = {E0:.6f} Ha")
print(f"V0   = {V0:.3f} Å³")
print(f"a0   = {V0**(1/3):.3f} Å")
print(f"B0   = {B0:.6f} Ha/Å³  = {B0_GPa:.2f} GPa")
print(f"B0'  = {B0p:.3f}")

# Generate smooth curve for the fit
V_fit = np.linspace(min(V_values), max(V_values), 200)
E_fit = birch_murnaghan(V_fit, *popt)

# Plot Energy vs Volume
plt.figure(figsize=(6,4))
plt.plot(V_values, energies, 'o', label="DFT data")
plt.plot(V_fit, E_fit, '-', label="Birch-Murnaghan fit")
plt.xlabel(r"Volume ($\AA^3$)")
plt.ylabel("Total Energy (Ha)")
plt.title("Equation of State: Energy vs Volume")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig("EOS_energy_vs_volume.png", dpi=300)
plt.show()
