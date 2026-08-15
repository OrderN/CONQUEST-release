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
