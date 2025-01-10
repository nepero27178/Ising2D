#!/usr/bin/julia

using DelimitedFiles
using Statistics

if length(ARGS) != 3
    println("How to use this program?
Type the following: \$ julia ./code/secondary-observables.jl arg1 arg2 arg3
Where:
· arg1 = lattice size
· arg2 = beta
· arg3 = number of resamples
Note that L is the linear size of lattice. 
This code's output is a single file, named L=<your-L>, and saved in 

		\"Ising2D/analysis/data/L=<your-L>.txt\"

The file is made of five columns: (1) Beta, (2) Binder's cumulant, (3) error
on Binder's cumulant, (4) magnetic susceptibility, (5) error on magnetic 
susceptibility.
")
    exit()
#if CurrentDirectory[end-9:end] != "processing"
#	println("Wrong directory! Please run this program inside 
#		\"Ising2D/processing\"
#")
#	exit()
else
    UserInput = ARGS
    L = parse(Int64,UserInput[1])	# Lattice Size
    Beta = UserInput[2]				# Temperature
    R = parse(Int64,UserInput[3])   # Resampling steps
end

# Functions to extract secondary observables and their errors

function GetSecondaryObservables(Data::Matrix{Float64})
	AvgAbsMag = mean(broadcast(abs,Data[:,2]))
	AvgMag2 = mean(Data[:,3])
	AvgMag4 = mean(Data[:,4])
	
	# Compute Binder's cumulant and magnetic susceptibility
	BinderCumulant = AvgMag4 / (AvgMag2 * AvgMag2)
	MagnetizationVariance = AvgMag2 - (AvgAbsMag * AvgAbsMag)
	return BinderCumulant, MagnetizationVariance
end

function ResampleData(Data::Matrix{Float64})
    FakeData = zeros(size(Data))
    for i in 1:size(Data,1)
        RandomIndex = rand(1:size(Data,1))
        FakeData[i,:] = Data[RandomIndex,:]
    end
    return FakeData
end

function GetBootstrapErrors(Data::Matrix{Float64}, R::Int64)
	# Get errors using bootstrap algorithm
	# R = Resampling steps
	FakeObservables = zeros(R,2)    # R rows, two columns
	
	for i in 1:R
		# Resample data to produce matrix of fake data
		FakeData = ResampleData(Data)
		FakeObservables[i,:] .= GetSecondaryObservables(FakeData)
	end
	
    # Compute std-dev over direction 1 (y direction)
    # First: Binder error. Second: suscpetibility error
	BootstrapErrors = std(FakeObservables,dims=1)
	return BootstrapErrors
end

# Main run

function main()
	""""
	Note: we take the data from the simulation folder, and not from the
	processing folder, because we noticed that there is no need for blocking
	(i.e. k=1 is enough).
	"""
    FilePathIn = "../simulations/data/L=$L/beta=" * Beta * ".txt"
    FilePathOut = "../analysis/data/L=$L.txt"
		
	"""
	U = Binder's Cumulant
	eU = error on Binder's Cumulant
	Chi = Magnetic Susceptibility
	eChi = error on Magnetic Susceptibility
	"""
	
	Data = readdlm(FilePathIn, ',', Float64, comments=true)
	U, MagnetizationVariance = GetSecondaryObservables(Data)
	
    # Susceptibility gets a factor (beta) L^D, ignore (beta)
	Chi = L^2 * MagnetizationVariance
	eU, eChi = GetBootstrapErrors(Data,R)

    open(FilePathOut, "a") do io
		writedlm(io, [Beta U eU Chi eChi], ',')
    end
end

main()
