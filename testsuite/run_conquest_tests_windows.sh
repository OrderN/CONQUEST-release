#!/bin/bash

# Windows (MSYS2 UCRT64) version of run_conquest_tests.sh.
#
# Runs every test_[0-9][0-9][0-9]_* directory with MS-MPI's mpiexec for
# each "NPxNT" mode given on the command line (NP = MPI processes,
# NT = OpenMP threads), then checks the results with pytest.  Outputs
# are kept in results_windows/<mode>/ so the modes can be compared.
#
# Usage (from an MSYS2 UCRT64 shell, in this directory):
#   ./run_conquest_tests_windows.sh                  # 1x1 4x1 1x4 2x2
#   ./run_conquest_tests_windows.sh 2x1 1x2          # custom modes

MODES=${*:-"1x1 4x1 1x4 2x2"}

# MS-MPI runtime (installed by msmpisetup.exe / winget Microsoft.msmpi)
MSMPI_BIN=${MSMPI_BIN:-"/c/Program Files/Microsoft MPI/Bin"}
export PATH="$(cygpath -u "$MSMPI_BIN"):$PATH"
export OMP_STACKSIZE=100M
export OPENBLAS_NUM_THREADS=1

CONQUEST=$(cd ../bin && pwd)/Conquest.exe
if [ ! -x "$CONQUEST" ]; then
    echo "Cannot find $CONQUEST - build it first with: (cd ../src; make SYSTEM=msys2)"
    exit 1
fi

summary=""
for mode in $MODES
do
    NP=${mode%x*}
    NT=${mode#*x}
    export OMP_NUM_THREADS=$NT
    out=results_windows/$mode
    mkdir -p "$out"
    echo "=== Mode $mode: $NP MPI processes x $NT OpenMP threads ==="
    for dn in $(ls -d test_[0-9][0-9][0-9]_*)
    do
        # Conquest needs at least one atom per MPI process, so cap the
        # number of processes at the atom count (line 4 of the coordinates)
        coords=$(awk 'tolower($1)=="io.coordinates" {print $2}' "$dn/Conquest_input")
        natoms=$(sed -n 4p "$dn/$coords" | awk '{print $1}')
        np=$NP
        note=""
        if [ -n "$natoms" ] && [ "$natoms" -lt "$NP" ]; then
            np=$natoms
            note=" (capped at $np processes: only $natoms atoms)"
        fi
        start=$(date +%s)
        (cd "$dn"; rm -f Conquest_out; mpiexec -n "$np" "$CONQUEST" > stdout.log 2>&1)
        rc=$?
        echo "  $dn: exit $rc, $(( $(date +%s) - start )) s$note"
        cp "$dn/Conquest_out" "$out/$dn.Conquest_out" 2>/dev/null
        mv "$dn/stdout.log" "$out/$dn.stdout.log"
    done
    python -m pytest -q test_check_output.py > "$out/pytest.log" 2>&1
    result=$(tail -1 "$out/pytest.log")
    echo "  pytest: $result"
    summary="$summary\n  $mode: $result"
done

echo -e "\n=== Summary ===$summary"
