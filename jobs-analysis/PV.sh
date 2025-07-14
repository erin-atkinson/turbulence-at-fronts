#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=40
#SBATCH --time=2:00:00
#SBATCH --job-name=ppPV
#SBATCH --output=../scratch/logs/ppPV.txt

module load julia/1.10.4
export JULIA_DEPOT_PATH=~/.julia-niagara
export JULIA_SCRATCH_TRACK_ACCESS=0
# Path to ramdisk to avoid too many writes on parallel filesystem
export RAM=/dev/shm/turbulence-at-fronts
rm $RAM -rf
mkdir $RAM

cd ~/turbulence-at-fronts

# Location of output.jld2
export SIM_OUTPUT_FOLDER=../scratch/turbulence-at-fronts/Strain
julia -t 40 -- src-analysis/post-process.jl $SIM_OUTPUT_FOLDER PV $RAM

#export SIM_OUTPUT_FOLDER=../scratch/Project/StrainQ1a2
#julia -t 40 -- src-analysis/post-process.jl $SIM_OUTPUT_FOLDER TKE $RAM o
#export SIM_OUTPUT_FOLDER=../scratch/Project/StrainQ2

#julia -t 40 -- src-analysis/post-process.jl $SIM_OUTPUT_FOLDER TKE $RAM o
#export SIM_OUTPUT_FOLDER=../scratch/Project/StrainQ3

#julia -t 40 -- src-analysis/post-process.jl $SIM_OUTPUT_FOLDER TKE $RAM o

rm $RAM -rf
