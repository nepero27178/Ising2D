using DelimitedFiles

# -----------------
# --- Bootstrap ---
# -----------------
  
function ResampleData(Data::Matrix{Float64})
		      
    """
    Resampling function:
    Input
    - Data matrix
    	- Rows number (size 1) = number of data point per observable
    	- Columns number (size 2) = number of observables
    Output
    - Array of resampled blocks, with repetition
    """
    ResampledData = zeros(size(Data)) # New matrix, same size as Data
    RowsNumber = size(Data,1)    
    
    for ResampledDataIndex in 1:RowsNumber
        RandomIndex = floor(Int,RowsNumber*rand()+1)
        ResampledData[ResampledDataIndex,:] = Data[RandomIndex,:]
    end
    
    return ResampledData
end
