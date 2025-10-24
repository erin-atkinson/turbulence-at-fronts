#!/bin/bash
# Run this in a login node to install Julia for Trillium

module load julia/1.10.10
export JULIA_DEPOT_PATH=$HOME/.julia-tri

julia -t 24 -- setup/setup-tri.jl
