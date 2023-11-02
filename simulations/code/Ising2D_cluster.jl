using Logging
using LinearAlgebra
include("./functions.jl")

# Setup

if length(ARGS) != 5
    println("How to use this program?
Type the following: \$ julia ./ising2D_metro.jl arg1 arg2 arg3 arg4 arg5
Where:
· arg1 = beta (float)
· arg2 = number of spins per side in your square lattice (int)
· arg3 = number of samples to be stashed due to thermalization (int)
· arg4 = number of Monte Carlo samples (int)
· arg5 = path to file where to save data (string)")
    exit()
else
    UserInput = ARGS
    Beta = parse(Float64,UserInput[1])                      # Temperature        
    LatticeSize = parse(Int64,UserInput[2])                 # Size
    ThermSamples = parse(Int64,UserInput[3])                # Thermalization
    MonteCarloSamples = parse(Int64,UserInput[4])           # Monte Carlo
    FilePath = UserInput[5]                                 # Data file
end

# Simulation

LatticeConfiguration = SetLattice(LatticeSize)
InitialSite = GetRandomSite(LatticeSize)
# Make sure the site index is a Vector
Cluster = CreateCluster(LatticeConfiguration, Beta, InitialSite, [], [])
FlipCluster(LatticeConfiguration, Cluster)
display(LatticeConfiguration)


function main()
    DataFile = open(FilePath,"w")
    @time begin
        # Set up lattice
        LatticeConfiguration = SetLattice(LatticeSize)
        
        Energy = GetEnergy(LatticeConfiguration)
        @info "Initial energy = $Energy"
        Magnetization = GetMagnetization(LatticeConfiguration)
        @info "Initial magnetization = $Magnetization"
        
        # Naked run: stash thermalization without any measurement, just update
        for ThermStep in 1:ThermSamples
            LatticeConfiguration = NextClusterStep(LatticeConfiguration,Beta)
        end
        # Actual run: simulate 2D-Ising
        for Step in 1:MonteCarloSamples

            # Calculate energy after N^2 updates not to have steps too short
            # in energy space
            for _ in 1:(LatticeSize^2)
                LatticeConfiguration = NextClusterStep(LatticeConfiguration,Beta)
            end
            # Calculate and collect energy
            Energy = GetEnergy(LatticeConfiguration)
            Magnetization = GetMagnetization(LatticeConfiguration)
            write(DataFile,"$Energy, $Magnetization\n")
        end
        
        close(DataFile)  
    end             # End for @time, printing the timestamp  
end                 # End main

main()
