# ------------------------------------------------------------------------------
# PART 1: Setup
# ------------------------------------------------------------------------------

# Set initial lattice

function SetLattice(L::Int64)::Matrix{Int8}
    
    """
    Create a lattice represented as a matrix of size (L,L), initialized to +1.
    """
        
    LatticeConfiguration = ones(Int8, L,L)
    return LatticeConfiguration
end

# Pre-calculate step acceptability in Markov Chain approach

function GetStepAcceptabilities(Topology::String, Beta::Float64)::Array{Float64}
    
    if Topology=="square"
    	StepAcceptabilities = GetStepAcceptabilitiesSquare(Beta)
    	
    elseif Topology=="triangular"
    	StepAcceptabilities = GetStepAcceptabilitiesTriangular(Beta)
    	
    elseif Topology=="hexagonal"
    	StepAcceptabilities = GetStepAcceptabilitiesHexagonal(Beta)
    
    else
    	throw(DomainError(Topology, "\nPlease choose one of the following possible topologies: \"square\", \"triangular\", \"hexagonal\".\n"))
    	
	end
	
	return StepAcceptabilities
end

# Pre-calculate neighbors for each site, based on topology

function PreCalculateNeighbours(Topology::String, L::Int64)::Array{Vector{Tuple{Int64, Int64}}}
    
    if Topology=="square"
    	AllNeighbours = PreCalculateNeighboursSquare(L)
    	
    elseif Topology=="triangular"
    	AllNeighbours = PreCalculateNeighboursTriangular(L)
    	
    elseif Topology=="hexagonal"
    	AllNeighbours = PreCalculateNeighboursHexagonal(L)
    
    else
    	throw(DomainError(Topology, "\nPlease choose one of the following possible topologies: \"square\", \"triangular\", \"hexagonal\".\n"))
    	
	end
	
	return AllNeighbours
end

# ------------------------------------------------------------------------------
# PART 2: Lattices with different topologies, geometrical properties
# ------------------------------------------------------------------------------

# Square lattice ---------------------------------------------------------------

function GetNeighboursSquare(L::Int64, Site::Tuple{Int64, Int64})
    
    """
    Get the indexes of the neighbouring sites to a given site, for the square 
    lattice, in clockwise order, starting from Up.
    """
    
    # Define neighbours. Remember that moving to the right, (x,y) → (x+1, y), means selecting 
    # the next column, (i,j) → (i,j+1). Therefore x is associated to columns.
    # Moving up, (x,y) → (x, y+1), means going to the previous row (i,j) → (i-1,j).

    i, j = Site

    Neighbours = [
        (mod(i-2, L) + 1, j),                    # Up (x, y+1)
        (i, mod(j, L) + 1),                      # Right (x+1, y)
        (mod(i, L) + 1, j),                      # Down (x, y-1)
        (i, mod(j-2, L) + 1)                     # Left (x-1, y)
    ]
    
    # Note that in modular algebra you have e.g. mod(3,3)=0, and 
    # since the indexing of Julia starts by 1 we have to indicate 
    # the position '(x+1)%N,y' as 'mod(x,N)+1,y'. Instead of 'mod(x-2,N)+1,y'
    # we use 'mod(x-2+N,N)+1,y' to indicate the position '(x-1)%N,y'
    # because modular algebra gets in trouble with negative numbers.

    return Neighbours
end

function PreCalculateNeighboursSquare(L::Int64)::Array{Vector{Tuple{Int64, Int64}}}

    """
    Pre-calculate the neighbors for all sites in the square lattice.
    Returns a 2D array where each element contains a list of neighboring sites 
    for that lattice site.
    """
    
    # Create an empty 2D array where each site has an array of neighbors
    AllNeighbours = Array{Vector{Tuple{Int64, Int64}}}(undef, L, L)
    
    # Populate the array 
    for i in 1:L
        for j in 1:L
            AllNeighbours[i, j] = GetNeighboursSquare(L, (i,j))
        end
    end
    
    return AllNeighbours
end

function GetStepAcceptabilitiesSquare(Beta::Float64)::Array{Float64}

	"""
    Since confrontation with an exponential value is necessary for every 
    Metropolis step, being the possible exponents (dE, i.e. delta energy) part 
    of a finite set (-4, -2, 0, 2, 4), we can calculate the exponential once and 
    for all and pick the correct value instead of calculating again. Also, being
    all the values (-4, -2, 0) automatically accepted, we only calculate the
    remaining ones.
    """
	
    PossibleEnergies = [2 4]
    StepAcceptabilities = exp.(-2 * Beta * PossibleEnergies)
    return StepAcceptabilities

end

# Triangular lattice -----------------------------------------------------------

function GetNeighboursTriangular(L::Int64, Site::Tuple{Int64, Int64})

    """
    Get the indexes of the neighbouring sites to a given site,
    for the triangular lattice, in clockwise order, starting from Up.
    """

    i, j = Site

    Neighbours = [
        (mod(i-2, L) + 1, j),                    # Up (x, y+1)
        (mod(i-2, L) + 1, mod(j, L) + 1),        # Upper Right (x+1, y+1)
        (i, mod(j, L) + 1),                      # Right (x+1, y)
        (mod(i, L) + 1, j),                      # Down (x, y-1)
        (mod(i, L)+ 1, mod(j-2, L) + 1),         # Lower Left (x-1, y-1)
        (i, mod(j-2, L) + 1)                     # Left (x-1, y)
    ]

    return Neighbours
end

function PreCalculateNeighboursTriangular(L::Int64)::Array{Vector{Tuple{Int64, Int64}}}

    """
    Pre-calculate the neighbors for all sites in the triangular lattice.
    Returns a 2D array where each element contains a list of neighboring sites
    for that lattice site.
    """
    
    # Create an empty 2D array where each site has an array of neighbors
    AllNeighbours = Array{Vector{Tuple{Int64, Int64}}}(undef, L, L)
    
    # Populate the array 
    for i in 1:L
        for j in 1:L
            AllNeighbours[i, j] = GetNeighboursTriangular(L, (i,j))
        end
    end
    
    return AllNeighbours
end

function GetStepAcceptabilitiesTriangular(Beta::Float64)::Array{Float64}

    PossibleEnergies = [2 4 6] # hex lattice
    StepAcceptabilities = exp.(-2 * Beta * PossibleEnergies)
    return StepAcceptabilities

end

# Hexagonal lattice ------------------------------------------------------------

function GetNeighboursHexagonal(L::Int64, Site::Tuple{Int64, Int64})

    """
    Get the indexes of the neighbouring sites to a given site, for the hexagonal
    lattice. Even rows store the spins of the first sublattice. Odd rows store
    the spins of the second sublattice.
    Note: L must be even.
    """

    i, j = Site

    if mod(i, 2) == 0 # even row
        Neighbours = [
            (mod(i-2, L) + 1, j),                    # Up (x, y+1)
            (mod(i, L) + 1, j),                      # Down (x, y-1)
            (mod(i-2, L) + 1, mod(j-2,L)+1),         # Upper Left (x+1, y+1)
        ]
    else # odd row
        Neighbours = [
            (mod(i-2, L) + 1, j),                    # Up (x, y+1)
            (mod(i, L) + 1, mod(j, L) + 1),          # Lower Right (x+1, y-1)
            (mod(i, L) + 1, j),                      # Down (x, y-1)
        ]
    end

    return Neighbours
end

function PreCalculateNeighboursHexagonal(L::Int64)::Array{Vector{Tuple{Int64, Int64}}}

    """
    Pre-calculate the neighbors for all sites in the hexagonal lattice.
    Returns a 2D array where each element contains a list of neighboring sites
    for that lattice site.
    """
    
    # Create an empty 2D array where each site has an array of neighbors
    AllNeighbours = Array{Vector{Tuple{Int64, Int64}}}(undef, L, L)
    
    # Populate the array 
    for i in 1:L
        for j in 1:L
            AllNeighbours[i, j] = GetNeighboursHexagonal(L, (i,j))
        end
    end
    
    return AllNeighbours
end

function GetStepAcceptabilitiesHexagonal(Beta::Float64)::Array{Float64}

    PossibleEnergies = [1 3] # hex lattice
    StepAcceptabilities = exp.(-2 * Beta * PossibleEnergies)
    return StepAcceptabilities
end

# ------------------------------------------------------------------------------
# PART 3: Metropolis and cluster update functions (all topology-indepenent)
# ------------------------------------------------------------------------------

function NextMetropolisStep!(L::Int64,
                             LatticeConfiguration::Matrix{Int8},
                             StepAcceptabilities::Array{Float64},
                             AcceptedSteps::Array{Int64},
                             AllNeighbours::Array{Vector{Tuple{Int64, Int64}}})
    
    """
    Implementation of the single Metropolis update algorithm.
    """

    # Extraction of a random site of the lattice
    Site = (rand(1:L), rand(1:L))

    # Reading of the current spin state of the site and the neighbours
    SiteSpin = LatticeConfiguration[Site...]
    Neighbours = AllNeighbours[Site...]

    TotalNeighboursSpin = 0
    
    for NeighbourSite in Neighbours
    	# Sum to total spin the i-th entry of the neighbouring spin to site x,y
    	TotalNeighboursSpin += LatticeConfiguration[NeighbourSite...]
    end
    
    DeltaEnergy = 2.0 * SiteSpin * TotalNeighboursSpin

    if DeltaEnergy <= 0 || rand() < StepAcceptabilities[Int(DeltaEnergy/4)] 
        LatticeConfiguration[Site...] = -1 * SiteSpin
        AcceptedSteps[1] += 1 
    end
    
    # Note: no return, since the function updates LatticeConfiguration
    
end

function NextClusterStep!(L::Int64,
						  LatticeConfiguration::Matrix{Int8}, 
						  ExpansionProbability::Float64,
						  AllNeighbours::Array{Vector{Tuple{Int64, Int64}}})
    
    """
    Next cluster step: select a random site and start growing the cluster.
    Note: the variable LatticeConfiguration gets updated when calling
    GrowCluster!.
    """
    
    # Extraction of a random site of the lattice
    StartingSite = (rand(1:L), rand(1:L))
    
    # Recursive cluster expansion from the initial site
    GrowCluster!(L, LatticeConfiguration, ExpansionProbability, StartingSite, 
    			 AllNeighbours)
    			 
    # Note: no return, since the function updates LatticeConfiguration
    			 
end

function GrowCluster!(L::Int64,
                      LatticeConfiguration::Matrix{Int8},
                      ExpansionProbability::Float64,
                      Site::Tuple{Int64, Int64},
                      AllNeighbours::Array{Vector{Tuple{Int64, Int64}}})
    """
    Key idea: we start at some site with spin +1 (for example).
    1) Save site spin as comparing variable for next steps;
    2) Update site spin (minimal cluster): now it is -1;
    3) Read neighbouring spin and compare them with old site spin value (+1).
       If they're equal and the step is accepted, expand the cluster to the new
       site repeating from step 1 with the new position.
    Note that since in the second recursive call LocalState is set to +1, having
    we already updated the first spin to -1, the cluster expansion will not
    touch it again. This mechanism prevents the cluster from going all over it
    already went.
    """
    
    SiteSpin = LatticeConfiguration[Site...]
    Neighbours = AllNeighbours[Site...]
    LatticeConfiguration[Site...] *= -1					# Flip the starting site
    
    for NeighbourSite in Neighbours
        if ( LatticeConfiguration[NeighbourSite...] == SiteSpin 
             && rand() < ExpansionProbability )
             # These conditions already imply NeighbourSite ∉ Cluster
        
             GrowCluster!(L, LatticeConfiguration, 
        		           ExpansionProbability, NeighbourSite, AllNeighbours)
        
        end
    end
    
    # Note: no return, since the function updates LatticeConfiguration

end

# ------------------------------------------------------------------------------
# PART 4: Physics
# ------------------------------------------------------------------------------

# Extract magnetization from lattice configuration (topology-independent)

function GetMagnetization(L::Int64, LatticeConfiguration::Matrix{Int8})::Float64
    
    """
    The magnetization is the sum of the matrix elements per unit volume.
    """
    
    return sum(LatticeConfiguration)/(L^2)
end

# Extract energy from lattice configuration

function GetEnergy(Topology::String, L::Int64, LatticeConfiguration::Matrix{Int8})::Float64

    if Topology=="square"
    	LatticeEnergy = GetEnergySquare(L, LatticeConfiguration)
    	
    elseif Topology=="triangular"
    	LatticeEnergy = GetEnergyTriangular(L, LatticeConfiguration)
    	
    elseif Topology=="hexagonal"
    	LatticeEnergy = GetEnergyHexagonal(L, LatticeConfiguration)
    
    else
    	throw(DomainError(Topology, "\nPlease choose one of the following possible topologies: \"square\", \"triangular\", \"hexagonal\".\n"))
    	
	end
	
	return LatticeEnergy
    
end

# Square lattice ---------------------------------------------------------------

function GetEnergySquare(L::Int64, LatticeConfiguration::Matrix{Int8})::Float64

    """
    In order to calculate the energy of the Ising lattice, we iterate over the
    sites and compute the sum. Note that S is the sum only over forward ones:
    this way we avoid repetition.
    """

    LatticeEnergy = 0.0
    for i in 1:L
        for k in 1:L
            j = L + 1 - k
            S0 = LatticeConfiguration[i,j]
            S = ( LatticeConfiguration[mod(i, L) + 1, j] +                  # Down (x, y-1)
                  LatticeConfiguration[i, mod(j-2, L) + 1] )                # Left (x-1, y)
                  
            NeighboursProduct = -S*S0
            LatticeEnergy += NeighboursProduct
        end
    end
    LatticeEnergy /= (L^2)
    return LatticeEnergy
end

# Triangular lattice -----------------------------------------------------------

function GetEnergyTriangular(L::Int64, LatticeConfiguration::Matrix{Int8})::Float64

    """
    In order to calculate the energy of the Ising lattice, we iterate over the
    sites and compute the sum. Note that S is the sum only over forward ones:
    this way we avoid repetition.
    """

    LatticeEnergy = 0.0
    for i in 1:L
        for k in 1:L
            j = L + 1 - k
            S0 = LatticeConfiguration[i,j]
            S = ( LatticeConfiguration[mod(i, L) + 1, j] +                  # Down (x, y-1)
                  LatticeConfiguration[mod(i, L)+ 1, mod(j-2, L) + 1 ] +    # Lower Left (x-1, y-1)
                  LatticeConfiguration[i, mod(j-2, L) + 1] )	            # Left (x-1, y)
                  
            NeighboursProduct = -S*S0
            LatticeEnergy += NeighboursProduct
        end
    end
    LatticeEnergy /= (L^2)
    return LatticeEnergy
end

# Hexagonal lattice ------------------------------------------------------------

function GetEnergyHexagonal(L::Int64, LatticeConfiguration::Matrix{Int8})::Float64

    """
    To calculate the energy of the hexagonal lattice, we simply sweep across 
    all sites of one sublattice, calculating its interaction with each spin with
    its 3 neighbours. 
    """

    LatticeEnergy = 0.0
    for i in 2:2:L # only even rows (only one sublattice)
        for j in 1:L
            S0 = LatticeConfiguration[i,j]
            S = ( LatticeConfiguration[mod(i-2, L) + 1, j] +                  # Up (x, y+1)
                  LatticeConfiguration[mod(i, L) + 1, j] +                    # Down (x, y-1)
                  LatticeConfiguration[mod(i-2, L) + 1, mod(j-2,L)+1] )       # Upper Left (x+1, y+1)
            
            # TODO Possibly equivalent: GetNeighboursHexagonal
            
            NeighboursProduct = -S*S0
            LatticeEnergy += NeighboursProduct
        end
    end
    LatticeEnergy /= (L^2)
    return LatticeEnergy
end
