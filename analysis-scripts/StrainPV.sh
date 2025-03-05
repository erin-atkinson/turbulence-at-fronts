#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=40
#SBATCH --time=2:00:00
#SBATCH --job-name=ppStrainPV
#SBATCH --output=../scratch/logs/ppStrainPV.txt

module load julia/1.10.4
export JULIA_DEPOT_PATH=~/.julia-niagara
export JULIA_SCRATCH_TRACK_ACCESS=0
# Path to ramdisk to avoid too many writes on parallel filesystem
mkdir /dev/shm/Project
export RAM=/dev/shm/Project

cd ~/Project

# Location of output.jld2
export SIM_OUTPUT_FOLDER=../scratch/Project/Strain

julia -t 40 -- src-analysis/post-process.jl $SIM_OUTPUT_FOLDER PV $RAM o
rm $RAM -rf
