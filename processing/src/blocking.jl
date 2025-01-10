#!/usr/bin/julia

using DelimitedFiles
using Statistics

# Setup

if length(ARGS) != 4
    println("How to use this program?
Type the following: \$ julia ./ising2D_metro.jl arg1 arg2 arg3 arg4
Where:
· arg1 = program setting 0/1 (see below)
· arg2 = lattice size
· arg3 = beta
· arg4 = block lengths (setting 0) / length (setting 1)
Note: input and output paths are generated automatically based on L
and beta.
---------------------------------------------------------------------
Program setting 0
An iterative search for different block lengths is performed, and for
each standard deviation of the magnetization is computed. The results
for given L, beta, and k (block length) are stored in the file

	../data/blocking_std_dev.txt
	
To set the different block lengths to check, use arg4 as a list of
strings:

e.g. arg4 = \"10 20 30 40\"

(do not forget the quotation marks!)
---------------------------------------------------------------------
Program setting 1
Data coming from the simulation file are blocked in blocks of length
equal to arg4. Quotation marks are not needed in this settings. Data
are saved in the corresponding folder inside processing/data/. 
")
    exit()
else
    UserInput = ARGS
    UserSelection = parse(Int64, UserInput[1])	# User selection
    L = UserInput[2]							# Lattice Size
    Beta = UserInput[3]							# Temperature
    UserLength = UserInput[4]					# Block length
end

# Blocking function

function BlockData(Data::Matrix{Float64}, BlockLength::Int64)
    "
    Reshapes the data by dividing them in blocks, then computes the mean value
    of each block and creates a new Matrix to store these mean values in.
    Input
        - Data: Matrix
        - BlockLength: Length of the block
    Output
        - Matrix
    "
    BlockNumber = Int(floor(size(Data, 1) / BlockLength)) # Number of blocks
    BlockedData = zeros(BlockNumber, size(Data, 2))
    for i in 1:BlockNumber
        iStart = (i - 1) * BlockLength + 1
        iEnd = i * BlockLength
        BlockedData[i, :] = mean(Data[iStart:iEnd,:], dims = 1)
    end
    return BlockedData
end

function main()
	
	FilePathIn = "../simulations/data/L=" * L * "/beta=" * Beta * ".txt"
	
    if Bool(UserSelection)
    	Data = readdlm(FilePathIn, ',', Float64, comments=true)
    	k = parse(Int64, UserLength)
    	BlockedData = BlockData(Data, k)
    	
    	FilePathOut = "./data/L=" * L * "/beta=" * Beta * ".txt"
		open(FilePathOut, "w") do io
            write(io, "# Energy, Magnetization, Magnetization2, Magnetization4\n")
            writedlm(io, BlockedData, ',')
    	end
    else
    	Data = readdlm(FilePathIn, ',', Float64, comments=true)
    	BlockLengths = [parse(Int64, i) for i in split(UserLength)]
    	
    	FilePathOut = "./data/blocking_std_dev.txt"
    	for k in BlockLengths
    		BlockedData = BlockData(Data, k)
    		MagStdDev = std(BlockedData[:, 2], corrected=true)
    		open(FilePathOut, "a") do io
    			writedlm(io, [L Beta k MagStdDev], ',')
    		end
    	end
    end
end

main()
