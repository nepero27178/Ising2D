# ----------------- #
# --- (1) Setup --- #
# ----------------- #

# Set initial lattice

function SetLattice(L::Int64)::Matrix{Int8}
    """Create a lattice represented as a matrix of size (L,L), initialized to +1."""
    LatticeConfiguration = ones(L,L)
    return LatticeConfiguration
end

# Pre-calculate step acceptability in Markov Chain approach

function GetStepAcceptabilities(Beta::Float64)
    """
    Since confrontation with an exponential value is necessary for every 
    Metropolis step, being the possible exponents (dE, i.e. delta energy) part 
    of a finite set (-4, -2, 0, 2, 4), we can calculate the exponential once and 
    for all and pick the correct value instead of calculating again. Also, being
    all the values (-4, -2, 0) automatically accepted, we only calculate the
    remaining ones.
    """
    # TODO: allow user to choose between square, [2 4], and triangular, [2 4 6]
    PossibleEnergies = [2 4 6] # triangular lattice
    StepAcceptabilities = exp.(-2 * Beta * PossibleEnergies)
    return StepAcceptabilities
end

# --------------------------------- #
# --- (2) Steps in Markov Chain --- #
# --------------------------------- #

# Select neighbours index on lattice

function GetNeighboursSquare(L::Int64, Site::Tuple{Int64, Int64})
    """
    Get the indexes of the neighbouring sites to a given site,
    for the square lattice, in clockwise order, starting from Up.
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

function GetNeighboursTriangle(L::Int64, Site::Tuple{Int64, Int64})
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

function PrecalculateNeighboursTriangle(L::Int64)
    """
    Pre-calculate the neighbors for all sites in the triangular lattice.
    Returns a 2D array where each element contains a list of neighboring sites for that lattice site.
    """
    
    # Create an empty 2D array where each site has an array of neighbors
    AllNeighbours = Array{Vector{Tuple{Int64, Int64}}}(undef, L, L)
    
    # Populate the array 
    for i in 1:L
        for j in 1:L
            AllNeighbours[i, j] = GetNeighboursTriangle(L, (i,j))
        end
    end
    
    return AllNeighbours
end

# Single Metropolis update

function NextMetropolisStep!(L::Int64,
                            LatticeConfiguration::Matrix{Int8},
                            StepAcceptabilities::Array{Float64},
                            AcceptedSteps::Array{Int64})
    """
    Implementation of the Metropolis algorithm.
    """

    # Extraction of a random site of the lattice
    Site = (rand(1:L), rand(1:L))

    # Reading of the current spin state of the site
    SiteSpin = LatticeConfiguration[Site...]

    TotalNeighboursSpin = 0
    NeighboursPositions = GetNeighboursTriangle(L,Site)
    
    for NeighbourSite in NeighboursPositions
    	# Sum to total spin the i-th entry of the neighbouring spin to site x,y
    	TotalNeighboursSpin += LatticeConfiguration[NeighbourSite...]
    end
    
    DeltaEnergy = 2.0 * SiteSpin * TotalNeighboursSpin

    if DeltaEnergy <= 0 || rand() < StepAcceptabilities[Int(DeltaEnergy/4)] 
        LatticeConfiguration[Site...] = -1 * SiteSpin
        AcceptedSteps[1] += 1 
        # # Remember to divide by volume (energy per site)
        # DeltaEnergy /= (L^2)
    end
end

# Cluster Metropolis update (Wolff's algorithm)

function NextClusterStep!(
            L::Int64,
            LatticeConfiguration::Matrix{Int8}, 
			ExpansionProbability::Float64,
            AllNeighbours::Array{Vector{Tuple{Int64, Int64}}})
    """Next cluster step. Note: the variable LatticeConfiguration gets updated when calling GrowCluster"""
    
    StartingSite = (rand(1:L), rand(1:L))
    
    GrowCluster!(L, LatticeConfiguration, ExpansionProbability, 
    			  StartingSite, AllNeighbours)
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
    # Neighbours = GetNeighboursTriangle(L, Site) # old version (before precalculation)
    Neighbours = AllNeighbours[Site...]
    LatticeConfiguration[Site...] *= -1 # flip
    
    for NeighbourSite in Neighbours
        if ( LatticeConfiguration[NeighbourSite...] == SiteSpin 
             && rand() < ExpansionProbability 
             # these conditions imply NeighbourSite ∉ Cluster already
             )
             GrowCluster!(L, LatticeConfiguration, 
        		           ExpansionProbability, NeighbourSite, AllNeighbours)
        end
    end
    # Note: no return, since the function updates LatticeConfiguration
end

# ------------------- #
# --- (3) Physics --- #
# ------------------- #

# Extract energy from lattice configuration

function GetEnergy(L::Int64, LatticeConfiguration::Matrix{Int8})::Float64

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
                  LatticeConfiguration[i, mod(j-2, L) + 1]                  # Left (x-1, y)
                  )
            NeighboursProduct = -S*S0
            LatticeEnergy += NeighboursProduct
        end
    end
    LatticeEnergy /= (L^2)
    return LatticeEnergy
end

# Extract magnetization from lattice configuration

function GetMagnetization(L::Int64, LatticeConfiguration::Matrix{Int8})::Float64
    """
    The magnetization is the sum of the matrix elements per unit volume.
    """
    return sum(LatticeConfiguration)/(L^2)
end
