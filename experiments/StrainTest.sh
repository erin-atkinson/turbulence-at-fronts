#!/bin/bash
#SBATCH --nodes=1
#SBATCH --gpus-per-node=1
#SBATCH --time=24:00:00
#SBATCH --job-name=Strain
#SBATCH --output=../scratch/logs/StrainTest.txt
#module load cuda/11.7
module load julia/1.10.4
export JULIA_DEPOT_PATH=$SCRATCH/.julia-mist
export JULIA_SCRATCH_TRACK_ACCESS=0
cd ~/Project

output_path=../scratch/Project/StrainTest2
run_time="3e2"
f="1e-4"
H="0.1"
Nx=384
Nh=256
Ny=256
Nz=8
Ro="100"
Ri="0.0"
alpha="1e-2"
Q="0"
c="0.5"
damping_width="3"
start_time="0"
save_time="1e0"

julia -t 8 -- src-fronts/simulation.jl $output_path $run_time $f $H $Nx $Nh $Ny $Nz $Ro $Ri $alpha $Q $c $damping_width $start_time $save_time
