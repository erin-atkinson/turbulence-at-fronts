#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=40
#SBATCH --time=0:45:00
#SBATCH --job-name=slices
#SBATCH --output=../scratch/logs/slices.txt

module load julia/1.10.4
export JULIA_DEPOT_PATH=~/.julia-niagara
export JULIA_SCRATCH_TRACK_ACCESS=0
# Path to ramdisk to avoid too many writes on parallel filesystem
export RAM=/dev/shm/turbulence-at-fronts
rm $RAM -rf
mkdir $RAM

cd ~/turbulence-at-fronts

julia -t 40 -- src-fig/save.jl

rm $RAM -rf
