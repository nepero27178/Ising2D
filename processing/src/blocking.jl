#!/usr/bin/julia

using DelimitedFiles

include("./functions.jl")

# Setup

if length(ARGS) != 4
    println("How to use this program?
Type the following: \$ julia ./ising2D_metro.jl arg1 arg2 arg3 arg4
Where:
· arg1 = beta (float)
· arg2 = number of spins per side in your square lattice (int)
· arg3 = path to file with the original data (string)
· arg4 = path to file where to save the blocked data (string)")
    exit()
else
    UserInput = ARGS
    Beta = parse(Float64,UserInput[1])                      # Temperature        
    LatticeSize = parse(Int64,UserInput[2])                 # Size
    FilePath = UserInput[3]	                                # Data file
    FilePathOut = UserInput[4]	                            # Output file
end

function main()

    Data = readdlm(FilePath, ',', Float64, comments=true)
    BlockedData = BlockData(Data, 200)

    open(FilePathOut, "a") do io
        write(io, "# Energy, Magnetization, Magnetization2, Magnetization4\n")
        writedlm(io, BlockedData, ',')
    end
    
end

main()
