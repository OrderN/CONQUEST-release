#!/bin/sh
#QSUB2 queue qD
#QSUB2 core 192
#QSUB2 mpi 192
#QSUB2 smp 1
#QSUB2 wtime 1:00:00
#PBS -N cSiSP-9k
cd $PBS_O_WORKDIR

nmpi=192
nsmp=1

source /etc/profile.d/modules.sh

#module load comp1 fftw
source ~/.module_load
cd $PBS_O_WORKDIR

# Set CQ
CQdir1=/home/jianbo/apps/CQ/develop/branch_f-mlff/CONQUEST-release_vec_r/bin/Conquest
CQdir=${CQdir1}

# Initial coordinates
nf_coord_in=coord.in_Si1000_perfect
cp ${nf_coord_in} coord.in

i=1
nf_temp=Conquest_input.Temp
comment=.Initial
type_pot=Si9000K

cp ${nf_temp}${comment} ${nf_temp}

#for nmpi in 2 4 8 16 32 48
#for nmpi in 48 32 24 16 12 8 4
for nmpi in ${nmpi}  #96 #48 #24 #4 8 12 16 24 48
do
  for temp in 9000 #9000 #9000 #300 900 1500 4000
  do
    #rm md.*
    #rm trajectory*
    #rm fort.*
    sed -e "s/@@@Temp/${temp}/g" Conquest_input.Temp > Conquest_input
    mpijob -mpi ${nmpi} -smp ${nsmp} ${CQdir} >& output.log.mpi${nmpi}.${temp}K
    mv Conquest_out   Conquest_out.md${i}${comment}.mpi${nmpi}.${temp}K.pot${type_pot}

    cp coord_next.dat coord.in.md${i}${comment}.mpi${nmpi}.${temp}K.pot${type_pot}
    cp coord_next.dat coord.in
    cp md.stats       md.stats.md${i}${comment}.mpi${nmpi}.${temp}K.pot${type_pot}
    mv trajectory.xyz trajectory.md${i}${comment}.mpi${nmpi}.${temp}K.xyz.pot${type_pot}
  done
done
