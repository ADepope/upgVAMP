#!/bin/bash

# SBATCH --reservation=hpcgrp_45
#SBATCH --job-name=run_VAMP_ukb22828_UKB_EST_v3_ldp08_fd_217_pheno_test
#SBATCH --time=0-00:40:00
#SBATCH --mem 20gb
#SBATCH --ntasks 2
#SBATCH --cpus-per-task 2
#SBATCH --output=VAMP_217_test.log
# SBATCH --exclude=leonid[64],bjoern[42]
# SBATCH --partition=gpu

module purge
ml gcc openmpi boost
module list 

# printing out the details of the job
bed_name=ukb22828_UKB_EST_v3_ldp08_fd_test
phen_name=ukb_ldp08_fd_R050_HT_test

echo "bed_name = " ${bed_name} 
echo "phen_name = " ${phen_name} 

export OMP_NUM_THREADS=2

vloc=/nfs/scistore13/robingrp/human_data/adepope_preprocessing/VAMPBirtyhday/gVAMP/cpp_vamp
bed_file_loc=/nfs/scistore13/robingrp/human_data/adepope_preprocessing/geno_testVAMP

mpic++ $vloc/main_real.cpp $vloc/vamp.cpp $vloc/utilities.cpp $vloc/data.cpp $vloc/options.cpp -march=native -Ofast -g -fopenmp -lstdc++fs -D_GLIBCXX_DEBUG -o  $vloc/main_real.exe

time mpirun -np 2 $vloc/main_real.exe --bed-file-test  ${bed_file_loc}/${bed_name}.bed \
                                --phen-files-test ${bed_file_loc}/${phen_name}.phen \
                                --N-test 15000 \
                                --N 399055 \
                                --Mt-test 2174071 \
                                --test-iter-range 1,35 \
                                --estimate-file /nfs/scistore13/robingrp/human_data/adepope_preprocessing/VAMPJune2023/gVAMP/estimates/x1_hat_HT_ukb22828_UKB_EST_v3_ldp08_fd_R040_pheno_200k_train_sim_it_6.bin \
                                --run-mode test