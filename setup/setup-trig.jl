using Pkg
@info "Setting up Julia depot for Trillium (GPU) at $DEPOT_PATH"
Pkg.add("CUDA")

Pkg.add("Oceananigans"; version="0.101.0")
Pkg.pin("Oceananigans")

Pkg.add("JLD2"; version="0.6.2")
Pkg.pin("JLD2")

Pkg.add("OffsetArrays")
Pkg.add("SpecialFunctions")

Pkg.precompile()

Pkg.status()
