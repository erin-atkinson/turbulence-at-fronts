#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=192
#SBATCH --time=6:00:00
#SBATCH --job-name=ppSLICES-Ly
#SBATCH --output=../scratch/logs/ppSLICES-Ly.txt

module load julia/1.10.10

# Copy installation to RAM disk
echo "Copying installation to RAM disk"
export RAM=/dev/shm/turbulence-at-fronts
mkdir $RAM

cp -r $HOME/turbulence-at-fronts/.julia-tri $RAM

# Launch from RAM disk
echo "Running..."
export JULIA_DEPOT_PATH=$RAM/.julia-tri
export JULIA_SCRATCH_TRACK_ACCESS=0
cd ~/turbulence-at-fronts

# Location of output.jld2
export SIM_OUTPUT_FOLDER=../scratch/turbulence-at-fronts/Strain-512
julia -t 192 -- src-analysis/postprocess/postprocess.jl $SIM_OUTPUT_FOLDER SLICES $RAM

#export SIM_OUTPUT_FOLDER=../scratch/turbulence-at-fronts/Strain-256
#julia -t 192 -- src-analysis/postprocess/postprocess.jl $SIM_OUTPUT_FOLDER SLICES $RAM

#export SIM_OUTPUT_FOLDER=../scratch/turbulence-at-fronts/Strain-128
#julia -t 192 -- src-analysis/postprocess/postprocess.jl $SIM_OUTPUT_FOLDER SLICES $RAM

#export SIM_OUTPUT_FOLDER=../scratch/turbulence-at-fronts/Strain-64
#julia -t 192 -- src-analysis/postprocess/postprocess.jl $SIM_OUTPUT_FOLDER SLICES $RAM


rm $RAM -rf
