#!/bin/bash

# SBATCH --exclusive
# SBATCH --reservation=mondegrp_35
#SBATCH --job-name=run_VAMP_ukb22828_UKB_EST_v3_ldp08_fd_217_pheno_train
#SBATCH --time=2-00:00:00
# SBATCH --mem 0
#SBATCH --mem-per-cpu=8gb
#SBATCH --ntasks 30
#SBATCH --cpus-per-task 2
# SBATCH --constraint=delta
# SBATCH --exclude=bjoern[48-56]
# SBATCH --array=1
#SBATCH --output=VAMP_ukb22828_UKB_EST_v3_ldp08_fd_217_pheno_train_gam2_from_Gibbs.log

module purge
ml gcc openmpi boost
module list 

# printing out the details of the job
bed_name=ukb22828_UKB_EST_v3_ldp08_fd_train
phen_name=ukb_ldp08_fd_R050_HT_train
Mt=2174071

# previously it was this:
# --probs 0.8,0.10,0.05,0.03,0.02 \
# --vars 0,0.00001,0.001,0.01,0.1 \

echo "bed_name = " ${bed_name} 
echo "phen_name = " ${phen_name} 

export OMP_NUM_THREADS=2

vloc=/nfs/scistore13/robingrp/human_data/adepope_preprocessing/VAMPJune2023/gVAMP
bed_file_loc=/nfs/scistore13/robingrp/human_data/adepope_preprocessing/geno_testVAMP

mpic++ $vloc/main_real.cpp $vloc/vamp.cpp $vloc/utilities.cpp $vloc/data.cpp $vloc/options.cpp -march=native -Ofast -g -fopenmp -lstdc++fs -D_GLIBCXX_DEBUG -o  $vloc/main_real.exe

time mpirun -np 30 $vloc/main_real.exe --bed-file  ${bed_file_loc}/${bed_name}.bed \
                                --phen-files ${bed_file_loc}/${phen_name}.phen \
                                --N 399055 \
                                --Mt ${Mt} \
                                --out-dir /nfs/scistore13/robingrp/human_data/adepope_preprocessing/output_testVAMP/estimates/ \
                                --out-name x1_hat_HT_${bed_name}_after8novars_gam2_from_Gibbs \
                                --iterations 35 \
                                --num-mix-comp 5 \
                                --CG-max-iter 15 \
                                --probs 0.97,0.01,0.01,0.005,0.005 \
                                --vars 0,0.0000001,0.000001,0.00001,0.0001 \
                                --model linear \
                                --run-mode infere \
                                --use-lmmse-damp 0 \
                                --store-pvals 1 \
                                --rho 0.5
