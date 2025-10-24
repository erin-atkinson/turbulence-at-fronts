#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=192
#SBATCH --time=0:45:00
#SBATCH --job-name=ppDFM
#SBATCH --output=../scratch/logs/ppDFM.txt

# Copy installation to RAM disk
export RAM=/dev/shm/turbulence-at-fronts
cp -r $HOME/.julia-tri $RAM

# Launch from RAM disk
export JULIA_DEPOT_PATH=$RAM/.julia-tri
export JULIA_SCRATCH_TRACK_ACCESS=0
cd ~/turbulence-at-fronts

# Location of output.jld2
export SIM_OUTPUT_FOLDER=../scratch/turbulence-at-fronts/Strain
julia -t 192 -- src-analysis/postprocess.jl $SIM_OUTPUT_FOLDER DFM $RAM

rm $RAM -rf
