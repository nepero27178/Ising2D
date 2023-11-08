#!/usr/bin/julia

using Logging
using LinearAlgebra
include("./include-jl/functions.jl")

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
    FilePath = UserInput[5]		    		    # Data file
end

# Simulation

function main()
    DataFile = open(FilePath,"w")
    @time begin
        # Set up lattice
        LatticeConfiguration = SetLattice(LatticeSize)
        
        # Initialize acceptance counters (to be normalized at the end)
        Acceptance = 0
        ThermAcceptance = 0
        
        Energy = GetEnergy(LatticeConfiguration)
        @info "Initial energy = $Energy"
        Magnetization = GetMagnetization(LatticeConfiguration)
        @info "Initial magnetization = $Magnetization"

        # Pre-calculate acceptation probabilities for every possibile state
        StepAcceptability = GetStepAcceptability(Beta)
        
        # Naked run: stash thermalization
        for ThermStep in 1:ThermSamples
            Output = NextMetropolisStep(LatticeConfiguration,Beta,
                                        StepAcceptability)
            
            # The first output is boolean: was the step accepted?             
            if Output[1]
                LatticeConfiguration = Output[2]
                ThermAcceptance+=1
            end
        end
        ThermAcceptance /= (ThermSamples * LatticeSize^2)
        @info "Thermal acceptance = $ThermAcceptance"
        
        # Real run: simulate 2D-Ising
        for Step in 1:MonteCarloSamples

            # Calculate energy after N^2 updates not to have steps too short
            # in energy space
            for _ in 1:(LatticeSize^2)
                Output = NextMetropolisStep(LatticeConfiguration,Beta,
                                            StepAcceptability)
		
                if Output[1]
                    LatticeConfiguration = Output[2]
                    Acceptance+=1
                end
            end
            # Calculate and store energy
            Energy = GetEnergy(LatticeConfiguration)
            Magnetization = GetMagnetization(LatticeConfiguration)
            write(DataFile,"$Energy, $Magnetization\n")
        end
        
        close(DataFile)
        Acceptance /= (MonteCarloSamples * LatticeSize^2)
        @info "Acceptance = $Acceptance"        
    end 			# End for @time, printing the timestamp  
end 				# End main

main()
