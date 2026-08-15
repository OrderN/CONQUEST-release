#!/bin/bash
set -e

DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$DIR"

CONQUEST_BIN="/home/augustin/Projects/CONQUEST-develop/bin/Conquest"

run_case() {
    local size_label=$1
    local nx=$2
    local np=$3
    local solver_mode=$4 # "scalapack" or "elpa_gpu"

    python3 gen_supercell.py $nx coords.dat

    local use_elpa="F"
    local use_gpu="F"
    if [ "$solver_mode" == "elpa_gpu" ]; then
        use_elpa="T"
        use_gpu="T"
    fi

    cat << EOF > Conquest_input
AtomMove.TypeOfRun static
IO.Coordinates coords.dat
IO.Iprint 4
IO.WriteOutToFile F
Grid.GridCutoff 80

DM.SolutionMethod diagon
Diag.MPMesh F
Diag.UseELPA $use_elpa
Diag.ELPA_GPU $use_gpu
Diag.ELPASolver ELPA1

minE.SCTolerance 1.0e-5
minE.MaxSCSCycles 2

General.NumberOfSpecies 1
%block ChemicalSpeciesLabel
 1 28.086 Si
%endblock
EOF

    local out_file="${size_label}_${solver_mode}_np${np}.out"
    echo "=========================================================="
    echo "Running $size_label ($solver_mode, np=$np)..."
    echo "=========================================================="
    mpirun -np $np $CONQUEST_BIN > $out_file 2>&1
    echo "--- Summary for $out_file ---"
    grep "Time taken for" $out_file || true
    grep "Total run time" $out_file || true
    echo ""
}

echo "Starting Bulk Silicon Benchmark Suite..."
# Run 1: 64 atoms (N_basis = 832)
run_case "Si_64" 2 2 "scalapack"
run_case "Si_64" 2 2 "elpa_gpu"

# Run 2: 216 atoms (N_basis = 2808)
run_case "Si_216" 3 2 "scalapack"
run_case "Si_216" 3 2 "elpa_gpu"
