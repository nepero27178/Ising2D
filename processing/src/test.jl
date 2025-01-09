#!/usr/bin/julia

using DelimitedFiles
using Statistics

FilePathIn = "/home/nepero27178/Documents/University/Metodi-Numerici/Mod1/Ising2D/simulations/data/L=10/beta=0.25.txt"
FilePathOut = "/home/nepero27178/Documents/University/Metodi-Numerici/Mod1/Ising2D/processing/data/L=10/beta=0.25.txt"
BlockLength = 100

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


Data = readdlm(FilePathIn, ',', Float64, comments=true)
BlockedData = BlockData(Data, BlockLength)

open(FilePathOut, "w") do io
    write(io, "# Energy, Magnetization, Magnetization2, Magnetization4\n")
    writedlm(io, BlockedData, ',')
end
