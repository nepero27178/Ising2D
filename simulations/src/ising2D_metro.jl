#!/usr/bin/julia

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
    write(DataFile,"# Energy, Magnetization, Magnetization2, Magnetization4 (per site)\n")
    @time begin
        # Set up lattice, with all spins up (+1)
        LatticeConfiguration = SetLattice(L)
        
        # Initialize accepted steps counters (to be normalized at the end).
        # Note: it's an array to make it modifiable by the function NextMetropolisStep
        # (Integers are immutable)
        AcceptedSteps = [0]
        AcceptedStepsThermalization = [0]
        
        Energy = GetEnergy(L, LatticeConfiguration)
        @info "Initial energy per site = $Energy"
        Magnetization = GetMagnetization(L, LatticeConfiguration)
        @info "Initial magnetization per site = $Magnetization"

        # Pre-calculate acceptation probabilities for every possibile state
        StepAcceptabilities = GetStepAcceptabilities(Beta)
        
        # Naked run: stash thermalization 
        for _ in 1:ThermN
            NextMetropolisStep!(L, LatticeConfiguration,
                                StepAcceptabilities, AcceptedStepsThermalization)
        end

        ThermalizationAcceptance = AcceptedStepsThermalization[1] /  (ThermN * L^2)
        @info "Thermalization acceptance = $ThermalizationAcceptance"
        
        # Real run: simulate 2D-Ising
        for i in 1:N
            # Calculate energy after N^2 updates not to have steps too short
            # in energy space
            for j in 1:(L^2)
                NextMetropolisStep!(L, LatticeConfiguration,StepAcceptabilities,AcceptedSteps)
            end
            # Calculate and store:
            # Energy and magnetization
            Energy = GetEnergy(L, LatticeConfiguration)
            Magnetization = GetMagnetization(L, LatticeConfiguration)
            
            # Energy^2, Energy^4, Magnetization^2, Magnetization^4
            # Magnetization^2, Magnetization^4
            Mag2 = Magnetization^2
            Mag4 = Magnetization^4
            write(DataFile,"$Energy, $Magnetization, $Mag2, $Mag4\n")
        end
        
        close(DataFile)

        Acceptance = AcceptedSteps[1] /  (N * L^2)

        @info "Acceptance = $Acceptance"        
    end 			# End for @time, printing the timestamp  
end 				# End main

main()
