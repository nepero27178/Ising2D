# ----------------- #
# --- (1) Setup --- #
# ----------------- #

# Set initial lattice

function SetLattice(L::Int64)
    """Create a lattice represented as a matrix of size (L,L), initialized to +1."""
    LatticeConfiguration = ones(L,L)
    return LatticeConfiguration
end

# Pre-calculate step acceptability in Markov Chain approach

function GetStepAcceptability(Beta::Float64)
    """
    Since confrontation with an exponential value is necessary for every 
    Metropolis step, being the possible exponents (dE, i.e. delta energy) part 
    of a finite set (-4, -2, 0, 2, 4), we can calculate the exponential once and 
    for all and pick the correct value instead of calculating again. Also, being
    all the values (-4, -2, 0) automatically accepted, we only calculate the
    remaining ones.
    """
    PossibleEnergies = [2 4]
    StepAcceptability = exp.(-2 * Beta * PossibleEnergies)
    # Remember: index Acceptances as Acceptances[dE/2]!
    return StepAcceptability
end

# --------------------------------- #
# --- (2) Steps in Markov Chain --- #
# --------------------------------- #

# Choose random site on lattice

function GetRandomSite(LatticeSize::Int64) # useless? No
    """
    Creates a Vector that indicates a random site of the lattice.
    Note that without [ ...] it would create a 1 X 2 Matrix.
    """
    RandomSite = [rand((1:LatticeSize),(1,2))...]
    return RandomSite
end

# Select neighbours index on lattice

function GetNeighboursSquare(LatticeConfiguration::Matrix{Float64},
		       Site::Array{Int64})
    """
    Get the indexes of the neighbouring sites to a given site,
    for the square lattice, in clockwise order, starting from Right.
    """
    L = size(LatticeConfiguration,1)
    
    x,y = Site
    Neighbours = [ [x,mod(y,L)+1],	# 1: one column after  (x+1)
    		   [mod(x-2+L,L)+1,y],	# 2: one row above     (y+1)
    		   [x,mod(y-2+L,L)+1],  # 3: one column before (x-1)
    		   [mod(x,L)+1,y] ]	# 4: one row below     (y-1)
    """
    List of all positions of neighbouring sites to given position. Note that in
    modular algebra you have e.g. mod(3,3)=0, and since the indexing of Julia
    starts by 1 we have to indicate the position '(x+1)%N,y' as 'mod(x,N)+1,y'.
    Instead of 'mod(x-2,N)+1,y' we use 'mod(x-2+N,N)+1,y' to indicate the
    position '(x-1)%N,y' because modular algebra gets in trouble with negative
    numbers.
    """    
    return Neighbours
end

function GetNeighboursTriangle(L::Int, site::Tuple{Int, Int})
    """
    Get the indexes of the neighbouring sites to a given site,
    for the triangular lattice, in clockwise order, starting from Up.
    """
    i, j = site

    # Define neighbors. Remember that moving to the right, (x,y) → (x+1, y), means selecting 
    # the next column, (i,j) → (i,j+1). Therefore x is associated to columns.
    # Moving up, (x,y) → (x, y+1), means going to the previous row (i,j) → (i-1,j).
    neighbors = [
        [mod(i-2, L) + 1, j],                    # Up (x, y+1)
        [mod(i-2, L) + 1, mod(j, L) + 1],        # Upper Right
        [i, mod(j, L) + 1],                      # Right (x+1, y)
        [mod(i, L) + 1, j],                      # Down (x, y-1)
        [mod(i, L)+ 1, mod(j-2, L) + 1],         # Lower Left
        [i, mod(j-2, L) + 1]                     # Left (x-1, y)
    ]
    return neighbors
end

# Single Metropolis update

function NextMetropolisStep(LatticeConfiguration::Matrix{Float64}, 
                            Beta::Float64, 
                            Acceptances::Array{Float64})
    """
    Implementation of the Metropolis-Rosenbluth-Teller algorithm.
    """
    L = size(LatticeConfiguration,1)
    Accepted = false

    # Extraction of a random site of the lattice
    Site = GetRandomSite(L)

    # Reading of the current spin state of the site
    SiteSpin = LatticeConfiguration[ Site[1], Site[2] ]

    TotalNeighboursSpin = 0
    NeighboursPositions = GetNeighbours(LatticeConfiguration,Site)
    
    for i in 1:4
    	# Sum to total spin the i-th entry of the neighbouring spin to site x,y
    	LocalX,LocalY = NeighboursPositions[i]
    	TotalNeighboursSpin += LatticeConfiguration[ LocalX, LocalY ]
    end
    
    DeltaEnergy = 2.0 * SiteSpin * TotalNeighboursSpin

    if DeltaEnergy <= 0 || rand() < Acceptances[Int(DeltaEnergy/4)] 
        LatticeConfiguration[ Site[1], Site[2] ] = -1 * SiteSpin
        DeltaMagnetization = -2 * SiteSpin
        Accepted = true

        # Remember to divide by volume in order to have energy and 
        # magnetization per unit volume
        DeltaEnergy /= (L^2)
        DeltaMagnetization /= (L^2)
    else 
        DeltaEnergy = 0
        DeltaMagnetization = 0
    end

    return Accepted, LatticeConfiguration
end

# Cluster Metropolis update (Wolff's algorithm)

function NextClusterStep(LatticeConfiguration::Matrix{Float64}, 
			 ExpansionProbability::Float64)
    
    L = size(LatticeConfiguration, 1)
    StartingSite = GetRandomSite(L)
    
    Cluster = GrowCluster(LatticeConfiguration, ExpansionProbability, 
    			  StartingSite, [])
    return LatticeConfiguration
end

function GrowCluster(LatticeConfiguration::Matrix{Float64},
                     ExpansionProbability::Float64,
                     Site::Vector{Int64},
                     Cluster::Vector{Any})
   
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
    
    L = size(LatticeConfiguration, 1)
    SiteSpin = LatticeConfiguration[Site...]
    Neighbours = GetNeighbours(Site,L)
    
    LatticeConfiguration[Site...] *= -1
    push!(Cluster, Site) 
    
    for NeighSite ∈ Neighbours
        if ( LatticeConfiguration[NeighSite...] == SiteSpin 
             && rand() < AcceptanceProbability )
             Cluster = GrowCluster(LatticeConfiguration, 
        		           ExpansionProbability, NeighSite, Cluster)
        end
    end
    return  Cluster
end

# ------------------- #
# --- (3) Physics --- #
# ------------------- #

# Extract energy from lattice configuration

function GetEnergy(LatticeConfiguration::Matrix{Float64})

    """
    In order to calculate the energy of the Ising lattice, we iterate over the
    sites and compute the sum. Note that S is the sum only over forward ones:
    this way we avoid repetition.
    """

    L = size(LatticeConfiguration,1)
    LatticeEnergy = 0.0
    for i in 1:L
        for j in 1:L
            S0 = LatticeConfiguration[i,j]
            S = ( LatticeConfiguration[mod(i,L)+1,j] +
                  LatticeConfiguration[i,mod(j,L)+1] )
            NeighboursProduct = -S*S0
            LatticeEnergy += NeighboursProduct
        end
    end
    LatticeEnergy /= (L^2)
    return LatticeEnergy
end

# Extract magnetization from lattice configuration

function GetMagnetization(LatticeConfiguration::Matrix{Float64})
    """
    The magnetization is the sum of the matrix elements per unit volume.
    """
    L = size(LatticeConfiguration,1)
    return sum(LatticeConfiguration)/(L^2)
end
