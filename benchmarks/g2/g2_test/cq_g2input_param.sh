#!/bin/bash
#
# **< lat >** 24Apr2015 
# 
boxsize=15.0
gridcutoff=60.0
functional=101

function cq_input() {

celldm=$boxsize

cat > $cqinput << EOF
IO.Title       $header 
IO.Coordinates $cqcoord  
IO.FractionalAtomicCoords F
IO.WriteOutToFile         F
IO.Iprint                 1
IO.BackTraceOn            F

## General Parameters
General.DifferentFunctional  T
General.FunctionalType       $functional 
General.NumberOfSpecies      $nspe 
General.PseudopotentialType  haman 
General.PartitionMethod      Hilbert
General.PAOFromFiles         T
General.NPartitionsX  1
General.NPartitionsY  1
General.NPartitionsZ  1

## Density matrix
DM.SolutionMethod     diagon 

##Integration Grid
Grid.GridCutoff       $gridcutoff 

## Moving Atoms
AtomMove.TypeOfRun static
AtomMove.TestAllForces F

EXX.GridSpacing  0.5
EXX.IntegRadius 15.0 

## Basis Sets
Basis.BasisSet         PAOs
Basis.SymmetryBreaking T

## Energy Minimisation
minE.SelfConsistent   T
minE.VaryBasis        F
minE.EnergyTolerance  1.0e-7
minE.LTolerance       1.0e-7
minE.SCTolerance      1.0e-7

SC.KerkerPreCondition T
SC.LinearMixingSC     T
SC.LinearMixingFactor 0.5
SC.LinearMixingEnd    1.0e-8
SC.KerkerFactor       0.01
SC.MaxIters           40
SC.MaxPulay           5

EOF

eval "$1=$celldm"

}


function cq_input_spin() {

celldm=$boxsize

cat > $cqinput << EOF
IO.Title       $header 
IO.Coordinates $cqcoord  
IO.FractionalAtomicCoords F
IO.WriteOutToFile         F
IO.Iprint                 1 
IO.BackTraceOn            F

## General Parameters
General.DifferentFunctional  T
General.FunctionalType       $functional 
General.NumberOfSpecies      $nspe 
General.PseudopotentialType  haman 
General.PartitionMethod      Hilbert
General.PAOFromFiles         T
General.NPartitionsX  1
General.NPartitionsY  1
General.NPartitionsZ  1

## Density matrix
DM.SolutionMethod    diagon 

##Integration Grid
Grid.GridCutoff      $gridcutoff 

## Moving Atoms
AtomMove.TypeOfRun       static
AtomMove.TestAllForces   F

## Spin Polarisation
Spin.SpinPolarised T    
Spin.FixSpin       F  
#Spin.NeUP          $spinup 
#Spin.NeDN          $spindn
Spin.Magn           $magn

EXX.GridSpacing   0.5
EXX.IntegRadius  15.0

## Basis Sets
Basis.BasisSet         PAOs
Basis.SymmetryBreaking T

## Energy Minimisation
minE.SelfConsistent   T
minE.VaryBasis        F
minE.EnergyTolerance  1.0e-7
minE.LTolerance       1.0e-7
minE.SCTolerance      1.0e-7

SC.KerkerPreCondition T
SC.LinearMixingSC     T
SC.LinearMixingFactor 0.5
SC.LinearMixingEnd    1.0e-8
SC.KerkerFactor       0.01
SC.MaxIters           40
SC.MaxPulay           5

EOF

eval "$1=$celldm"

}

