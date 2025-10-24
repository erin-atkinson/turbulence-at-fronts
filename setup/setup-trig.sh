#!/bin/bash
# Run this in a login node to install Julia for Trillium (GPU)

module load julia/1.10.10
export JULIA_DEPOT_PATH=$HOME/.julia-trig

julia -t 24 -- setup/setup-trig.jl
