#!/usr/bin/julia

using Logging
using LinearAlgebra
include("./functions.jl")

# Setup

if length(ARGS) != 6
    println("How to use this program?
Type the following: \$ julia ./ising2D_metro.jl arg1 arg2 arg3 arg4 arg5 arg6
Where:
· arg1 = topology (string, choose: \"square\" / \"triangular\" / \"hexagonal\")
· arg2 = beta (float)
· arg3 = number of spins per side in your square lattice (int)
· arg4 = number of samples to be stashed due to thermalization (int)
· arg5 = number of Monte Carlo samples (int)
· arg6 = path to file where to save data (string)")
    exit()
else
    UserInput = ARGS
    Topology = UserInput[1]						  # Lattice topology
    Beta = parse(Float64,UserInput[2])            # Temperature        
    L = parse(Int64,UserInput[3])                 # Size
    ThermN = parse(Int64,UserInput[4])            # Thermalization
    N = parse(Int64,UserInput[5])                 # Monte Carlo
    FilePath = UserInput[6]                       # Data file
end

# Simulation

function main()
    @info "Working on a $Topology lattice with L=$L, β=$Beta, N=$N"
    @time begin
        # Set up lattice
        LatticeConfiguration = SetLattice(L)
        ExpansionProbability = 1 - exp(-2*Beta)					# Cluster expansion probability
        AllNeighbours = PreCalculateNeighbours(Topology, L)		# Speedup
        
        Energy = GetEnergy(Topology, L, LatticeConfiguration)
        @info "Initial energy = $Energy"
        Magnetization = GetMagnetization(L, LatticeConfiguration)
        @info "Initial magnetization = $Magnetization"
        
        # Naked run: stash thermalization without any measurement, just update,
        # don't measure thermal acceptance
        for _ in 1:ThermN
            NextClusterStep!(L, LatticeConfiguration, ExpansionProbability, AllNeighbours)
        end
        
        DataFile = open(FilePath,"a")
        
        # Actual run: simulate 2D-Ising
        for i in 1:N

            # Calculate energy after N^2 updates not to have steps too short
            # in energy space
            
            for j in 1:1 # (L^2) # TODO CHANGE 
                NextClusterStep!(L, LatticeConfiguration, ExpansionProbability, AllNeighbours)
            end
            
            # Calculate and store:
            # Energy and magnetization
            Energy = GetEnergy(Topology, L, LatticeConfiguration)
            Magnetization = GetMagnetization(L, LatticeConfiguration)
            
            # |Magnetization|, Magnetization^2, Magnetization^4
            AbsMag = abs(Magnetization)
            Mag2 = Magnetization^2
            Mag4 = Magnetization^4
            
            # Write on file
            write(DataFile,"$Energy, $AbsMag, $Mag2, $Mag4\n")
        end
        
        close(DataFile)  
    end             # End for @time, printing the timestamp  
end                 # End main

main()
