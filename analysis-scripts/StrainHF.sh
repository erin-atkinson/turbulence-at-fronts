#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=40
#SBATCH --time=8:30:00
#SBATCH --job-name=ppStrainHF
#SBATCH --output=../scratch/logs/ppStrainHF.txt

module load julia/1.10.4
export JULIA_DEPOT_PATH=~/.julia-niagara
export JULIA_SCRATCH_TRACK_ACCESS=0
# Path to ramdisk to avoid too many writes on parallel filesystem
mkdir /dev/shm/Project
export RAM=/dev/shm/Project

cd ~/Project

# Location of output.jld2
export SIM_OUTPUT_FOLDER=../scratch/Project/StrainHF

julia -t 40 -- src-analysis/post-process.jl $SIM_OUTPUT_FOLDER DIVSTEADY $RAM o
julia -t 40 -- src-analysis/post-process.jl $SIM_OUTPUT_FOLDER KSTEADY $RAM o
julia -t 40 -- src-analysis/post-process.jl $SIM_OUTPUT_FOLDER PVSTEADY $RAM o
julia -t 40 -- src-analysis/post-process.jl $SIM_OUTPUT_FOLDER ENERGY $RAM o
julia -t 40 -- src-analysis/post-process.jl $SIM_OUTPUT_FOLDER DFM $RAM o
julia -t 40 -- src-analysis/post-process.jl $SIM_OUTPUT_FOLDER DIV $RAM o
