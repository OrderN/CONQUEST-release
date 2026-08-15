#!/usr/bin/env python3
import sys

# 8-atom diamond cubic unit cell in fractional coordinates
a0 = 10.2625335 # Bohr (5.431 Angstrom)
basis = [
    [0.00, 0.00, 0.00],
    [0.00, 0.50, 0.50],
    [0.50, 0.50, 0.00],
    [0.50, 0.00, 0.50],
    [0.25, 0.25, 0.25],
    [0.25, 0.75, 0.75],
    [0.75, 0.75, 0.25],
    [0.75, 0.25, 0.75],
]

def make_supercell(nx, ny, nz, filename="coords.dat"):
    coords = []
    for ix in range(nx):
        for iy in range(ny):
            for iz in range(nz):
                for b in basis:
                    pos = [
                        (b[0] + ix) / nx,
                        (b[1] + iy) / ny,
                        (b[2] + iz) / nz
                    ]
                    coords.append(pos)
    natoms = len(coords)
    Lx = nx * a0
    Ly = ny * a0
    Lz = nz * a0
    with open(filename, 'w') as f:
        f.write(f'{Lx:14.6f} {0.0:14.6f} {0.0:14.6f}\n')
        f.write(f'{0.0:14.6f} {Ly:14.6f} {0.0:14.6f}\n')
        f.write(f'{0.0:14.6f} {0.0:14.6f} {Lz:14.6f}\n')
        f.write(f'{natoms}\n')
        for c in coords:
            f.write(f'{c[0]:14.6f} {c[1]:14.6f} {c[2]:14.6f}  1 T T T\n')
    print(f'Generated {filename}: {nx}x{ny}x{nz} supercell, {natoms} atoms, N_basis = {natoms*13}')

if __name__ == "__main__":
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 2
    out = sys.argv[2] if len(sys.argv) > 2 else "coords.dat"
    make_supercell(n, n, n, out)
