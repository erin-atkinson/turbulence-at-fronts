#!/bin/bash
#SBATCH --nodes=1
#SBATCH --gpus-per-node=1
#SBATCH --time=2:30:00
#SBATCH --job-name=vortex
#SBATCH --output=../scratch/logs/front-growth/vortex.txt
module load cuda/11.0.3

cd ~/Project
# output_folder run_time a f H Nx Ny Nz Ro Ri α damping_rate damping_frac

./../julia-1.8.5/bin/julia -t 8 --project="env" -- src-vortices/simulation.jl ../scratch/Project/vortex-strain 200000 6000 "1e-4" 100 1000 1000 100 "0.6" "2" "1e-5" "3e-3" "1"
