#!/bin/sh
#QSUB2 queue qA
#QSUB2 core 48
#QSUB2 mpi 48
#QSUB2 smp 1
#QSUB2 wtime 1:00:00
#PBS -N cSi5-test
cd $PBS_O_WORKDIR

nmpi=24
nsmp=1

source /etc/profile.d/modules.sh

#module load comp1 fftw
source ~/.module_load

cd $PBS_O_WORKDIR

# Set CQ
# With Energy and Stress
CQdir1=/home/jianbo/apps/CQ/develop/branch_f-mlff/CONQUEST-release_vec_r/bin/Conquest
CQdir=${CQdir1}
# Stable without Stress
#CQdir=/home/jianbo/apps/CQ/develop/CONQUEST-release_mlff_v3.0/bin/Conquest

i=1
save_dir=./diff_mpi
save_dir=./

mkdir ${save_dir}

numsteps=100
cq_input=Conquest_input.Temp.numsteps

#for nmpi in 1 2 4 8 16 32 48
for nmpi in 48
do
  for temp in 4000 #300 900 1500 4000
  do
    #rm md.*
    #rm trajectory.*
    #rm fort.*
    sed -e "s/@@@Temp/${temp}/g" -e "s/@@@numsteps/${numsteps}/g" ${cq_input} > Conquest_input
    mpijob -mpi ${nmpi} -smp ${nsmp} ${CQdir} >& output.log.mpi${nmpi}.${temp}K.numsteps${numsteps}
    mv Conquest_out ${save_dir}/Conquest_out.mpi${nmpi}.${temp}K.numsteps${numsteps}
    cp coord_next.dat coord.in.mpi${nmpi}.${temp}K.numsteps${numsteps}
    cp md.stats ${save_dir}/md.stats.mpi${nmpi}.${temp}K.numsteps${numsteps}
    mv trajectory.xyz ${save_dir}/trajectory.mpi${nmpi}.${temp}K.xyz
    #rm out*log
    #rm *matrix*
  done
done
