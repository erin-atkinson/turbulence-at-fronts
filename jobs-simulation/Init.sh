#!/bin/bash
#SBATCH --nodes=1
#SBATCH --gpus-per-node=1
#SBATCH --time=24:00:00
#SBATCH --job-name=Init
#SBATCH --output=../scratch/logs/Init.txt
#module load cuda/11.7
module load julia/1.10.4
export JULIA_DEPOT_PATH=$SCRATCH/.julia-mist
export JULIA_SCRATCH_TRACK_ACCESS=0
cd ~/Project

output_path=../scratch/Project/StrainDamping2
run_time="4e5"
f="1e-4"
H="100"
Nx=1350
Nh=1000
Ny=128
Nz=128
Ro="0.4"
Ri="2.0"
alpha="1e-5"
Q="100"
c="0.5"
damping_width="3"
start_time="-4e5"
save_time="1e3"

julia -t 8 -- src-fronts/simulation.jl $output_path $run_time $f $H $Nx $Nh $Ny $Nz $Ro $Ri $alpha $Q $c $damping_width $start_time $save_time
