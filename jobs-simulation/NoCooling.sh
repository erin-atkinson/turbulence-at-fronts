#!/bin/bash
#SBATCH --nodes=1
#SBATCH --gpus-per-node=1
#SBATCH --time=10:30:00
#SBATCH --job-name=NoCooling
#SBATCH --output=../scratch/logs/NoCooling.txt
#module load cuda/11.7
module load julia/1.10.4
export JULIA_DEPOT_PATH=$SCRATCH/.julia-mist
export JULIA_SCRATCH_TRACK_ACCESS=0
cd ~/turbulence-at-fronts

output_path=../scratch/turbulence-at-fronts/NoCooling
run_time="8e5"
f="1e-4"
H="100"
Nx=1024
Ny=32
Nz=128
Ro="0.4"
Ri="2.0"
alpha="1e-5"
Q="0"
c="0.5"
start_time="0"
save_time="1e3"
s="1.05"

julia -t 8 -- src-simulation/simulation.jl $output_path $run_time $f $H $Nx $Ny $Nz $Ro $Ri $alpha $Q $c $start_time $save_time $s "no cooling"
