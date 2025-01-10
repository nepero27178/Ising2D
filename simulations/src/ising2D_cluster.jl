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
    Beta = parse(Float64,UserInput[1])            # Temperature        
    L = parse(Int64,UserInput[2])                 # Size
    ThermN = parse(Int64,UserInput[3])            # Thermalization
    N = parse(Int64,UserInput[4])                 # Monte Carlo
    FilePath = UserInput[5]                       # Data file
end

# Simulation

function main()
    @info "Working on L=$L, β=$Beta, N=$N"
    DataFile = open(FilePath,"w")
    write(DataFile,"# Energy, Magnetization, Magnetization2, Magnetization4\n")
    @time begin
        # Set up lattice
        LatticeConfiguration = SetLattice(L)
        ExpansionProbability = 1 - exp(-2*Beta)
        
        Energy = GetEnergy(L, LatticeConfiguration)
        @info "Initial energy = $Energy"
        Magnetization = GetMagnetization(L, LatticeConfiguration)
        @info "Initial magnetization = $Magnetization"
        
        # Naked run: stash thermalization without any measurement, just update,
        # don't measure thermal acceptance
        for _ in 1:ThermN
            NextClusterStep(L, LatticeConfiguration, ExpansionProbability)
        end
        
        # Actual run: simulate 2D-Ising
        for i in 1:N

            # Calculate energy after N^2 updates not to have steps too short
            # in energy space
            
            for j in 1:(L^2) 
                NextClusterStep(L, LatticeConfiguration, ExpansionProbability)
            end
            
            # Calculate and store:
            # Energy and magnetization
            Energy = GetEnergy(L, LatticeConfiguration)
            Magnetization = GetMagnetization(L, LatticeConfiguration)
            
            # Magnetization^2, Magnetization^4
            Mag2 = Magnetization^2
            Mag4 = Magnetization^4
            write(DataFile,"$Energy, $Magnetization, $Mag2, $Mag4\n")
        end
        
        close(DataFile)  
    end             # End for @time, printing the timestamp  
end                 # End main

main()
