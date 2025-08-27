# Turbulence at fronts
A component of my PhD thesis at The University of Toronto: The pathway to frontal arrest in the submesoscale mixed layer

## Simulation
LES simulations of fronts in the ocean using `Oceananigans.jl`. The aim is to be able to represent a submesoscale mixed layer with forcing via cooling and a background deformation flow.

## Science
`Write a little about the why, after we have figured out a story to tell`

## File setup
### src-simulation
Scripts to create a simulation in Oceananigans. Split into components to make life easier.
### src-analysis
Post-processing scripts and helper functions to generate data to be plotted once simulations are run. Post-processing scripts generate all the data to be plotted for ease of recreation later.
### src-fig
### experiments
Shell scripts that run each case.
### analysis
Notebooks used for visualisation, testing and general messing around. May also have some `.jl` files just for figure setup
