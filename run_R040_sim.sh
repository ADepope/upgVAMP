#!/bin/bash

# SBATCH --reservation=hpcgrp_45
#SBATCH --job-name=run_VAMP_sim_June2023_R040
#SBATCH --mem-per-cpu=4gb
#SBATCH --time=1-00:00:00
#SBATCH --ntasks 12
#SBATCH --cpus-per-task 2
#SBATCH --output=VAMP_June2023_R040_sim_rho09_19.log

module purge
ml gcc openmpi boost
module list 

# printing out the details of the job
bed_name=ukb22828_UKB_EST_v3_ldp08_fd_R040_pheno_200k_train
phen_name=ukb22828_UKB_EST_v3_ldp08_fd_R040_pheno_200k_train
Mt=521208
Mt=326760 

echo "bed_name = " ${bed_name} 
echo "phen_name = " ${phen_name} 

export OMP_NUM_THREADS=2

vloc=/nfs/scistore13/robingrp/human_data/adepope_preprocessing/VAMPJune2023/gVAMP
bed_file_loc=/nfs/scistore13/robingrp/human_data/adepope_preprocessing/geno_testVAMP

mpic++ ${vloc}/sim.cpp $vloc/vamp.cpp $vloc/utilities.cpp $vloc/data.cpp $vloc/options.cpp -march=native -Ofast -g -fopenmp -lstdc++fs -D_GLIBCXX_DEBUG -o  ${vloc}/sim.exe

time mpirun -np 12 ${vloc}/sim.exe --bed-file  ${bed_file_loc}/${bed_name}.bed \
                                --N 200000 \
                                --Mt ${Mt} \
                                --out-dir /nfs/scistore13/robingrp/human_data/adepope_preprocessing/VAMPJune2023/gVAMP/estimates/ \
                                --out-name x1_hat_HT_${bed_name}_sim \
                                --iterations 25 \
                                --h2 0.5 \
                                --CV 50000 \
                                --num-mix-comp 6 \
                                --CG-max-iter 10 \
                                --probs 0.8,0.1,0.03,0.03,0.02,0.02 \
                                --vars 0,0.000001,0.00001,0.0001,0.001,0.01 \
                                --model linear \
                                --run-mode infere \
                                --rho 0.9 \
                                --store-pvals 1
