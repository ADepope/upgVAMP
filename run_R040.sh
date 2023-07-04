#!/bin/bash

# SBATCH --reservation=hpcgrp_45
#SBATCH --job-name=run_VAMP_June2023_R040
#SBATCH --mem-per-cpu=2gb
#SBATCH --time=0-06:00:00
#SBATCH --ntasks 10
#SBATCH --cpus-per-task 2
#SBATCH --output=VAMP_June2023_R040_9.log

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

mpic++ $vloc/main_real.cpp $vloc/vamp.cpp $vloc/utilities.cpp $vloc/data.cpp $vloc/options.cpp -march=native -Ofast -g -fopenmp -lstdc++fs -D_GLIBCXX_DEBUG -o  $vloc/main_real.exe

time mpirun -np 10 $vloc/main_real.exe --bed-file  ${bed_file_loc}/${bed_name}.bed \
                                --phen-files ${bed_file_loc}/${phen_name}.phen \
                                --N 200000 \
                                --Mt ${Mt} \
                                --out-dir /nfs/scistore13/robingrp/human_data/adepope_preprocessing/VAMPJune2023/gVAMP/estimates/ \
                                --out-name x1_hat_HT_${bed_name}_8 \
                                --iterations 25 \
                                --num-mix-comp 6 \
                                --CG-max-iter 10 \
                                --probs 0.8,0.1,0.03,0.03,0.02,0.02 \
                                --vars 0,0.000001,0.00001,0.0001,0.001,0.01 \
                                --model linear \
                                --run-mode infere \
                                --rho 0.6 \
                                --store-pvals 1
