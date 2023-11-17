using DelimitedFiles
using Statistics


function ReadData(Filename)
    "
    Reads the data and returns a Matrix
    "
    Result = readdlm(Filename)
    return Result 
end

# ---------------
# --- Physics ---
# ---------------

# Binder's cumulant

function GetBinderCumulant(Mag2::Array{Float64}, Mag4::Array{Float64})
	# Computes Binder's cumulant <m^4>/<m^2>^2 starting from data
	MeanMag2 = mean(Mag2)
	MeanMag4 = mean(Mag4)
	return MeanMag4 / ( (MeanMag2) * (MeanMag2) )
end

# (Artificial) magnetization variance

function GetMagnetizationVariance(Mag::Array{Float64}, Mag2::Array{Float64})
	# Computes magnetization "variance" (<m^2>-<|m|>^2) starting from data
	MeanAbsMag = mean(abs.(Mag))
	MeanMag2 = mean(Mag2)
	return MeanMag2 - ( MeanAbsMag * MeanAbsMag )

# ----------------
# --- Blocking ---
# ----------------

function BlockData(Data::Matrix{Float64}, BlockLenght::Int64)
    "
    Reshapes the data by dividing them in blocks, then computes the mean value
    of each block and creates a new Matrix to store these mean values in.
    Input
        - Data: Matrix
        - BlockLenght: Length of the block
    Output
        - Matrix
    "
    BlockNumber = Int(floor(size(Data, 1)/BlockLenght))
    BlockedData = zeros(BlockNumber, size(Data, 2))
    for BlockIndex in 1:BlockNumber
        BlockStart = (BlockIndex-1)*BlockLenght +1
        BlockEnd = BlockIndex*BlockLenght
        BlockMean = mean(Data[BlockStart:BlockEnd,:],dims = 1)
        BlockedData[BlockIndex, :] = BlockMean
    end
    return BlockedData
end

# ------------------
# --- Resampling ---
# ------------------
  
function ResampleData(Data::Matrix{Float64})
		      
    """
    Resampling function:
    Input
    - Data matrix
    	- Rows number (size 1) = number of data point per observable
    	- Columns number (size 2) = number of observables
    Output
    - Array of resampled data, with repetition
    """
    ResampledData = zeros(size(Data)) # New matrix, same size as Data
    RowsNumber = size(Data,1)    
    
    for ResampledDataIndex in 1:RowsNumber
        RandomIndex = floor(Int,RowsNumber*rand()+1)
        ResampledData[ResampledDataIndex,:] = Data[RandomIndex,:]
    end
    
    return ResampledData
end
