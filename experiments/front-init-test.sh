#!/bin/bash
#SBATCH --nodes=1
#SBATCH --gpus-per-node=1
#SBATCH --time=3:00:00
#SBATCH --job-name=front-initialisation-test
#SBATCH --output=../scratch/logs/front-growth/front.txt
module load cuda/11.0.3

cd ~/Project

# Ri=0.3 Ro=0.6 1000x64x100 took 15 minutes to do 10/f
./../julia-1.8.5/bin/julia -t 8 --project="env" -- src-fronts/simulation.jl ../scratch/Project/front-init-test 0 "1e-4" 100 1000 64 100 "0.6" "0" "0" "10" "2" "0.5" 100000
