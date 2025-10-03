#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=192
#SBATCH --time=13:00:00
#SBATCH --job-name=ppTKE
#SBATCH --output=../scratch/logs/ppTKE.txt

module load julia/1.10.10
export JULIA_DEPOT_PATH=$SCRATCH/.julia-tri
export JULIA_SCRATCH_TRACK_ACCESS=0
# Path to ramdisk to avoid too many writes on parallel filesystem
export RAM=/dev/shm/turbulence-at-fronts
rm $RAM -rf
mkdir $RAM

cd ~/turbulence-at-fronts

# Location of output.jld2
export SIM_OUTPUT_FOLDER=../scratch/turbulence-at-fronts/Strain
julia -t 192 -- src-analysis/postprocess.jl $SIM_OUTPUT_FOLDER TKE $RAM

rm $RAM -rf
