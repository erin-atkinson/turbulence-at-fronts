#!/bin/bash
#SBATCH --nodes=1
#SBATCH --gpus-per-node=1
#SBATCH --time=24:00:00
#SBATCH --job-name=Init
#SBATCH --output=../scratch/logs/Init.txt

export JULIA_DEPOT_PATH=$SCRATCH/.julia-trig
export JULIA_SCRATCH_TRACK_ACCESS=0
cd ~/Project

output_path=../scratch/turbulence-at-fronts/Strain
run_time="1e5"
f="1e-4"
H="100"
Nx=2048
Ny=128
Nz=256
Ro="0.1"
Ri="2.0"
alpha="1e-5"
Q="100"
c="0.5"
start_time="-4e5"
save_time="1e3"
s="1.0173"

julia-1.10.10/bin/julia -t 24 -- src-fronts/simulation.jl $output_path $run_time $f $H $Nx $Nh $Ny $Nz $Ro $Ri $alpha $Q $c $damping_width $start_time $save_time
