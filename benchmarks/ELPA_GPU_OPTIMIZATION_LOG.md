# CONQUEST ELPA GPU Optimization & Debugging Log

**Branch:** `feature/elpa-gpu-diag-opt`  
**Target Hardware:** Machine `waltz` (2 × NVIDIA GeForce RTX 2080 Ti, Intel Xeon W-2225 @ 4.10GHz, CUDA 12.9, OpenMPI 4.1.6, GCC 13.3.0)

---

## [2026-08-15 01:20:00 UTC] Point 1: ELPA GPU Compilation, Interface Modernization & Bug Fixes

### 1. Problem Description & Root Cause
1. **Symbol Conflict**: In `src/ELPAModule.f90`, the module variable `integer :: elpa_API` collided with the ELPA library module `elpa_api` imported via `use elpa`, leading to compilation errors (`Error: Symbol 'elpa_module' cannot have a type`).
2. **Typo in Parameter Name**: Parameter setting `blacs_cotext` instead of `blacs_context` caused runtime failure or ignored BLACS context.
3. **Broken / Duplicate GPU Setup**: `ELPAModule.f90` called `elp%setup()` before solver selection, and then called deprecated `elp%setup_gpu()`.
4. **Build Environment on Waltz**: ELPA 2024.05.001 required compilation with CUDA support (`--enable-nvidia-gpu-kernels`, `--with-NVIDIA-GPU-compute-capability=sm_75`, `-lcublas`, `-lcusolver`, `-lcudart`, `-lstdc++`).

### 2. Changes Made
* **`src/ELPAModule.f90` & `src/ELPAModuleDUMMY.f90`**:
  - Renamed `elpa_API` $\to$ `elpa_api_version`.
  - Fixed `blacs_context`.
  - Replaced duplicate setup with unified GPU property configuration (`gpu` / `nvidia-gpu`) before `elp%setup()`.
  - Added support for `GPU` and `NVIDIA_GPU` in `Diag.ELPA2Kernel`.
* **`src/initial_read_module.f90`**:
  - Updated input parser for `elpa_api_version` and valid `Diag.ELPA2Kernel` options.
* **`src/system/system.waltz.make`**:
  - Added build configuration linking against `/home/augustin/local/elpa-gpu/lib/libelpa_openmp.so` and CUDA 12.9 libraries.

### 3. Benchmark Verification
* **Test Case**: `benchmarks/water_64mols` (64 $H_2O$ molecules, $N = 1,088$ basis functions, SZP basis, 2 MPI ranks).
* **Timings per SCF iteration**:

| Method | Pass 1: Eval Diag (`mode='N'`) | Pass 2: Evec Diag (`mode='V'`) | Total Diag Time / SCF Step | Total Run Time (7 SCF steps) |
| :--- | :--- | :--- | :--- | :--- |
| **ScaLAPACK CPU** | $1,805\text{ ms}$ | $3,721\text{ ms}$ | $5,526\text{ ms}$ | $136.5\text{ s}$ |
| **ELPA CPU (1-stage)** | $1,780\text{ ms}$ | $3,650\text{ ms}$ | $5,430\text{ ms}$ | $138.6\text{ s}$ |
| **ELPA GPU (1-stage)** | $1,699\text{ ms}$ | $1,888\text{ ms}$ | $3,587\text{ ms}$ | $113.6\text{ s}$ |

* **Verification Status**: PASSED. Eigenvector solve is **$2.0\times$ faster on GPU**, and total diagonalization time per SCF iteration reduced from **$5.53\text{ s}$ to $3.59\text{ s}$**.

---

## [2026-08-15 10:35:00 UTC] Point 2: Single-Pass `FindEvals` (Eliminating Redundant Pass 1)

### 1. Problem Description & Root Cause
* In original CONQUEST, `FindEvals` executed a two-pass algorithm on every SCF iteration:
  1. **Pass 1 (`mode = 'N'`)**: Performed full generalized diagonalization ($H C = \epsilon S C$) across all k-points just to find the eigenvalue spectrum and compute $E_F$ / occupancies (`findFermi`).
  2. **Pass 2 (`mode = 'V'`)**: Re-executed full generalized diagonalization from scratch across all k-points to compute both eigenvalues AND eigenvectors $C$, discarding the Pass 1 eigenvalues.
* Generalized solvers (`pzhegvx` and `ELPA_zhegv`) with `mode = 'V'` already compute both eigenvalues and eigenvectors simultaneously in a single call.
* As a consequence, $\approx 50\%$ of total diagonalization time and data distribution was completely redundant.

### 2. Changes Made
* **`src/DiagModule.f90`**:
  - Expanded module variable `z` from 3D `(row_size, col_size, nspin)` to 4D `(row_size, col_size, nkpoints_max, nspin)` so eigenvectors for all k-points within a process group are preserved.
  - Updated `initDiag` and `endDiag` to allocate and track `z(row_size, col_size, nkpoints_max, nspin)`.
  - In `FindEvals`: changed the first loop to `call distrib_and_diag(spin, i, 'V', .true.)`, computing both eigenvalues and eigenvectors in a single pass.
  - Replaced the redundant second diagonalization loop with a direct assembly pass: reading precomputed `z(:,:,i,spin)` into `DistributeSC_to_ref` and `buildK`.
  - Updated DeltaSCF local excitation routines to reuse `z(:,:,1,spin)`.
* **`src/Makefile`**:
  - Fixed `sed` delimiter in dependency rule from `/` to `|` (`s|BBB|...|`) to support git branch names containing slashes (e.g. `feature/elpa-gpu-diag-opt`).

### 3. Benchmark Verification & Scaling Results

#### 3.1 Water 64 Molecules ($N_{\text{basis}} = 1,088$, SZP Basis, $\Gamma$-point, 2 MPI ranks)

| Configuration | Solver Time / SCF Step | Assembly (`buildK`) / Step | Total Diag Time / Step | Speedup vs Baseline |
| :--- | :--- | :--- | :--- | :--- |
| **Baseline ScaLAPACK (2-Pass)** | $5,526\text{ ms}$ (Pass 1: 1805ms, Pass 2: 3721ms) | N/A (included in Pass 2) | $5,526\text{ ms}$ | $1.00\times$ (ref) |
| **ScaLAPACK (Single-Pass)** | $3,390\text{ ms}$ | $377\text{ ms}$ | $3,767\text{ ms}$ | **$1.47\times$** |
| **ELPA GPU (Single-Pass)** | **$1,570\text{ ms}$** | **$376\text{ ms}$** | **$1,946\text{ ms}$** | **$2.84\times$** |

#### 3.2 Bulk Silicon Supercells ($\Gamma$-point, DZP Basis, 13 basis fns/atom, 2 MPI ranks)

| System | $N_{\text{basis}}$ | ScaLAPACK CPU Diag Time | ELPA GPU Diag Time | GPU Solver Speedup |
| :--- | :--- | :--- | :--- | :--- |
| **Bulk Si 64 atoms** ($2\times 2\times 2$) | $832$ | $1,547\text{ ms}$ | **$1,129\text{ ms}$** | **$1.37\times$** |
| **Bulk Si 216 atoms** ($3\times 3\times 3$) | $2,808$ | $51,401\text{ ms}$ ($51.4\text{ s}$) | **$7,320\text{ ms}$ ($7.3\text{ s}$)** | **$7.02\times$** |

* **Key Takeaway**: As matrix size grows from $N = 832 \to 1,088 \to 2,808$, the GPU speedup escalates dramatically from **$1.37\times \to 2.15\times \to 7.02\times$** due to the $O(N^3)$ computational density on GPU tensor cores.
* **Verification Status**: PASSED. Single-pass output on standard testsuite is bit-for-bit identical to reference.

---

