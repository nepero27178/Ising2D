#!/usr/bin/julia

using DelimitedFiles
using Statistics
using Dates

ErrorMessage = "How to use this program?
Type the following: \$ julia ./ising2D_metro.jl arg1 arg2 arg3 arg4
Where:
· arg1 = lattice size
· arg2 = beta
· arg3 = block length (program setting --use-optimal) / lengths (program setting --try) (see below)
· arg4 = input data filepath (raw)
· arg5 = output data filepath
· arg6 = program setting: --use-optimal / --try (see below)
Note: input and output paths are generated automatically based on L
and beta.
---------------------------------------------------------------------
Program setting 0
An iterative search for different block lengths is performed, and for
each standard deviation of the magnetization is computed. The results
for given L, beta, and k (block length) are stored in the file

	../data/std-dev-analysis/blocking_std_dev.txt
	
To set the different block lengths to check, use arg3 as a list of
strings:

e.g. arg3 = \"10 20 30 40\"

(do not forget the quotation marks!)
---------------------------------------------------------------------
Program setting 1
Data coming from the simulation file are blocked in blocks of length
equal to arg4. Quotation marks are not needed in this settings. Data
are saved in the corresponding folder inside processing/data/. 
"

# Setup

if length(ARGS) != 6
    println(ErrorMessage)
    exit()
else
    UserInput = ARGS
    L = UserInput[1]							# Lattice Size
    Beta = UserInput[2]							# Temperature
    UserLengthString = UserInput[3]				# Block length (use optimal) or trial block lengths (use trial)
    InputDataFilepath = UserInput[4]  			# Input data filepath
    OutputDataFilepath = UserInput[5]  			# Output data filepath
   	UserSelection = UserInput[6]				# User selection
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
		
    if UserSelection == "--use-optimal"
    
    	"""
    	Use optimal mode: here we assume already to have computed the optimal
    	block size, thus UserLengthString is a simple Int64 input.
    	"""
    	
    	Data = readdlm(InputDataFilepath, ',', Float64, comments=true)
    	OptimalBlockLength = parse(Int64, UserLengthString)
    	BlockedData = BlockData(Data, OptimalBlockLength)
    	
		open(OutputDataFilepath, "w") do io
            write(io, "# Energy, Magnetization, Magnetization2, Magnetization4\n")
            writedlm(io, BlockedData, ',')
    	end
    
    elseif UserSelection == "--try"
    	
    	"""
    	Use trial mode: here we need to compute many different trial lengths.
    	User input is a string of different trial lengths, delimited in prompt
    	by quotation marks \" length1 length2 ... \". Thus BlockLength input
    	is a string type object, 
    	"""
    	
    	Data = readdlm(InputDataFilepath, ',', Float64, comments=true)
    	BlockLengths = [parse(Int64, i) for i in split(UserLengthString)]
    	
    	for TrialBlockLength in BlockLengths
    
    		BlockedData = BlockData(Data, TrialBlockLength)
    		EnergyStdDev = std(BlockedData[:, 1], corrected=true) / sqrt(size(BlockedData, 1))
    		MagStdDev = std(BlockedData[:, 2], corrected=true) / sqrt(size(BlockedData, 1))
    		Mag2StdDev = std(BlockedData[:, 3], corrected=true) / sqrt(size(BlockedData, 1))
    		Mag4StdDev = std(BlockedData[:, 4], corrected=true) / sqrt(size(BlockedData, 1))
    
    		open(OutputDataFilepath, "a") do io
    			writedlm(io, [L Beta TrialBlockLength EnergyStdDev MagStdDev Mag2StdDev Mag4StdDev], ',')
    		end
    
    	end
    else
    	print(ErrorMessage)
    	exit()
    end
end

main()
