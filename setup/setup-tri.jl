using Pkg
@info "Setting up Julia depot for Trillium at $DEPOT_PATH"

Pkg.add("Oceananigans"; version="0.99.1")
Pkg.pin("Oceananigans")

Pkg.add("JLD2"; version="0.6.2")
Pkg.pin("JLD2")

Pkg.add("OffsetArrays")
Pkg.add("SpecialFunctions")

Pkg.precompile()

Pkg.status()
