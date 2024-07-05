#!/bin/sh
#QSUB2 queue qD
#QSUB2 core 48 
#QSUB2 mpi 48 
#QSUB2 smp 1 
#QSUB2 wtime 1:00:00
#PBS -N cSi8
cd $PBS_O_WORKDIR

npcores=48
nnode=1
nsmp=1
let nmpi=npcores*nnode/nsmp

nmpi=1

source /etc/profile.d/modules.sh

#module load comp1 fftw
source ~/.module_load
cd $PBS_O_WORKDIR

# Set CQ
CQdir1=/home/jianbo/apps/CQ/develop/branch_f-mlff/CONQUEST-release_vec_r/bin/Conquest
CQdir=${CQdir1}

# Set coordinate file
coord_dir='./coords'
#coord_type='shift1Atom'
coord_type='perfect'
sys_name="Si8_cubic_${coord_type}"

# Set potential
flag_MLFF=T #F
pp_dir='./pp'
pot_Split="Si_all_Split2b3b" 
pot_BP2b="Si_all_BP2b"
pot_type=$pot_Split

# bohr unit of la=5.4309 Angstrom
la=10.262912871373274

#for la_ratio in $(seq -50 5 50)
for la_ratio in $(seq -5 1 5)
do
  loc_ratio=$(echo "scale=2; 1.00 + 0.01 * $la_ratio" | bc | awk '{printf("%.2f\n",$1)}')
  loc_la=$(echo "scale=12; $la * $loc_ratio" | bc)
  echo $loc_ratio $loc_la
  for pot_type in $pot_type
  do
    cp ${pp_dir}/ml_pot_${pot_type}.txt  ml_pot.txt

    for runtype in static # static sqnm cg lbfgs
    do
      for optcell in F # T # F
      do
        if [ "${flag_MLFF}" == "T" ]; then
          loc_cmt="${runtype}_${coord_type}_Optcell${optcell}_ML${flag_MLFF}_pot${pot_type}"
        else
          loc_cmt="${runtype}_${coord_type}_Optcell${optcell}"
        fi
#       loc_log_cmt=${loc_cmt}_la${loc_ratio}
        loc_log_cmt=la${loc_ratio}
        nf_log=Conquest_out.${loc_log_cmt} 
  
        ## update title IO.Title
        loc_title="Si8_cubic_${runtype}_Optcell${optcell}"
        sed -i "/IO.Title/ s/.*/IO.Title ${loc_title}/" Conquest_input
        sed -i "/AtomMove.TypeOfRun/ s/.*/AtomMove.TypeOfRun ${runtype}/" Conquest_input
        sed -i "/AtomMove.OptCell/ s/.*/AtomMove.OptCell  ${optcell}/" Conquest_input
        sed -i "/General.MLFF/ s/.*/General.MLFF ${flag_MLFF}/" Conquest_input

        echo ${nf}
        # cp coord.in.${sys_name} coord.in
        sed -e "s/${la}/${loc_la}/g" ${coord_dir}/coord.in.${sys_name}_tmp > coord.in
        
        # run cq
        mpijob -mpi ${nmpi} -smp ${nsmp} ${CQdir} >& output_${loc_log_cmt}.log

        # save results
        cp Conquest_out   Conquest_out.${loc_log_cmt}
        cp coord_next.dat coord.in.${loc_log_cmt}
        mv md.stats       md.stats.${loc_log_cmt}
        mv trajectory.xyz trajectory.${loc_log_cmt}.xyz
  
      done
    done
  done
done
