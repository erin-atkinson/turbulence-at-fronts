#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=40
#SBATCH --time=3:30:00
#SBATCH --job-name=ppNoStrain
#SBATCH --output=../scratch/logs/ppNoStrain.txt

module load julia/1.10.4
export JULIA_DEPOT_PATH=~/.julia-niagara
export JULIA_SCRATCH_TRACK_ACCESS=0
# Path to ramdisk to avoid too many writes on parallel filesystem
export RAM=/dev/shm/Project
mkdir $RAM

cd ~/Project

# Location of output.jld2
export SIM_OUTPUT_FOLDER=../scratch/Project/NoCooling

#julia -t 40 -- src-analysis/post-process.jl $SIM_OUTPUT_FOLDER DFM $RAM o
julia -t 40 -- src-analysis/post-process.jl $SIM_OUTPUT_FOLDER PV $RAM o
#julia -t 40 -- src-analysis/post-process.jl $SIM_OUTPUT_FOLDER ENERGY $RAM o
#julia -t 40 -- src-analysis/post-process.jl $SIM_OUTPUT_FOLDER PVSILLY $RAM o

rm $RAM -rf
