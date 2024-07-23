#!/bin/bash
#SBATCH --nodes=1
#SBATCH --gpus-per-node=1
#SBATCH --time=17:00:00
#SBATCH --job-name=NIechoes-Ek
#SBATCH --output=../scratch/logs/NIechoes-Ek.txt
#module load cuda/11.7
module load julia/1.10.4
export JULIA_DEPOT_PATH=$SCRATCH/.julia-mist
export JULIA_SCRATCH_TRACK_ACCESS=0
cd ~/Project

# Ri=0.3 Ro=0.6 1000x64x100 took 15 minutes to do 10/f
#./../julia-1.8.5/bin/julia -t 8 --project="env" -- src-fronts/simulation.jl ../scratch/Project/front-init-test 0 "1e-4" 100 1024 64 128 "0.6" "0" "0" "10" "2" "1" 1000000
julia -t 8 -- src-fronts/simulation.jl ../scratch/Project/NIechoes-Ek 0 "1e-4" 100 1024 512 128 "0.8" "0.0" "0" "0.001" "0" "2" "1" 500000