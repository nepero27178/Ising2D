include("./functions.jl")
CurrentDirectory = pwd()

if length(ARGS) != 1
    println("How to use this program?

Run this program inside:
		\"Ising2D/processing\"
since here will be stored all blocked data.

Type the following: \$ julia ./code/secondary-observables.jl arg
Where:
· arg = path to blocked data files folder
Note that L is the linear size of lattice. This program can take as input the
entire folder <whatever-is-before>/L=<your-L>. For example, to analyze the
entire L=10 folder, just type
\$ julia ./code/secondary-observables.jl ./data/L=10

This code's output is a single file, named L=<your-L>, and saved in 
		\"Ising2D/analysis/data/L=<your-L>.txt\"

The file is made of five columns: (1) Beta, (2) Binder's cumulant, (3) error
on Binder's cumulant, (4) magnetic susceptibility, (5) error on magnetic 
susceptibility.
")
    exit()
if CurrentDirectory[end-9:end] != "processing"
	println("Wrong directory! Please run this program inside 
		\"Ising2D/processing\"
")
	exit()
else
    InputDirectory = ARGS # Folder to data files
end

# Functions to extract secondary observables and their errors

function GetSecondaryObservables(BlockedData::Matrix{Float64})
	Mag = BlockedData[:,2]
	Mag2 = BlockedData[:,3]
	Mag4 = BlockedData[:,4]
	
	# Compute Binder's cumulant and magnetic susceptibility
	BinderCumulant = GetBinderCumulant(Mag2,Mag4)
	MagneticSusceptibility = GetMagnetizationVariance(Mag,Mag2)
	return BinderCumulant, MagnetizationVariance
end

function GetBootstrapErrors(BlockedData::Matrix{Float64},
							ResamplingSteps::Int64)
	# Get errors using bootstrap algorithm
	R = ResamplingSteps
	FakeObservables = zeros(R,2) # Two columns
	BootstrapErrors = zeros(2)
	
	for i in 1:R
		# Resample data to produce matrix of fake data
		FakeData = ResampleData(BlockedData)
		
		# Select the correct columns
		Mag = FakeData[:,2]
		Mag2 = FakeData[:,3]
		Mag4 = FakeData[:,4]
		
		# Calculate the fake observables per each resample
		FakeObservables[i,1] = GetBinderCumulant(Mag2,Mag4)
		FakeObservables[i,2] = L^2 * GetMagnetizationVariance(Mag,Mag2)
	end
	
	BootstrapErrors = var(FakeObservables,dims=1)
	return BootstrapErrors
end

# Main run

function main()

	# Change directory to input, and list files inside it
	L = parse(Int64, InputDirectory[end-2:end])
	cd(InputDirectory)
	BetaFiles = readdir()
	
	# File path to output data: change "processing" to "analysis" and add ".txt"	
	FilePathOut = replace(pwd(), "processing" => "analysis")
	FilePathOut = FilePathOut * ".txt"
	FileOut = open(FilePathOut,"w")
	write(FileOut, "# Beta, Binder, eBinder, Susceptibility, eSusceptibility\n")
	
	for i in 1:lenght(BetaFiles)
	
		# Extract single Beta
		BetaFile = BetaFiles[i]
		Beta = parse(Float64, BetaFile[6:end])
		
		"""
		BC = Binder's Cumulant
		eBC = error on Binder's Cumulant
		MS = Magnetic Susceptibility
		eMS = error on Magnetic Susceptibility
		"""
		
		DataIn = open(BetaFile,"r")
		BC, MagnetizationVariance = GetSecondaryObservables(DataIn)
		MS = L^2 * MagnetizationVariance
		
		R=100 # Resampling steps, change?
		eBC, eMS = GetBootstrapErrors(DataIn,R)
		
		write(FileOut, "# $Beta, $BC, $eBC, $MS, $eMS\n")
	end
	
	close(FileOut)
	cd("../..")
end

main()
