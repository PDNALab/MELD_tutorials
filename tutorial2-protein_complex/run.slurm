#!/bin/bash
#SBATCH --job-name=meld_complex
#SBATCH --output=out.out
#SBATCH --error=err.err
##SBATCH --mail-type=ALL
##SBATCH --mail-user=g.jokent@gmail.com
#SBATCH --time=100:00:00
#SBATCH --ntasks=30
#SBATCH --gpus-per-task=1
#SBATCH --cpus-per-gpu=1
#SBATCH --distribution=cyclic:cyclic
#SBATCH --mem-per-cpu=2000mb
#SBATCH --qos=alberto.perezant
#SBATCH --partition=gpu
#SBATCH --constraint=2080ti

##############################################

module purge

source ~/.load_OpenMMv8_cuda12-4
srun --mpi=pmix_v3 launch_remd --debug
