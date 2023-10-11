function SetLattice(LatticeSize::Int64) # useless?
    LatticeConfiguration = ones(LatticeSize,LatticeSize)
    return LatticeConfiguration
end

function GetRandomSite(LatticeSize::Int64) # useless?
    RandomSitePosition = rand((1:LatticeSize),(1,2))
    return RandomSitePosition
end

function GetEnergy(LatticeConfiguration::Matrix{Float64})
    """
    In order to calculate the energy of the Ising lattice, we iterate over the
    sites and compute the sum. Note that S is the sum only over forward ones:
    this way we avoid repetition.
    """
    N = size(LatticeConfiguration,1)
    LatticeEnergy = 0.0
    for i in 1:N
        for j in 1:N
            S0 = LatticeConfiguration[i,j]
            S = ( LatticeConfiguration[mod(i,N)+1,j] +
                  LatticeConfiguration[i,mod(j,N)+1] )
            NeighboursProduct = -S*S0
            LatticeEnergy += NeighboursProduct
        end
    end
    LatticeEnergy /= (N^2)
    return LatticeEnergy
end

function GetMagnetization(LatticeConfiguration::Matrix{Float64})
    """
    The magnetization is the sum of the matrix elements per unit volume.
    """
    N = size(LatticeConfiguration,1)
    return sum(LatticeConfiguration)/(N^2)
end

# Exponentiation (once every beta)

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

# Single metropolis step

function NextMetropolisStep(LatticeConfiguration::Matrix{Float64}, 
                           Beta::Float64, 
                           Acceptances::Array{Float64})
    """
    Implementation of the Metropolis-Rosenbluth-Teller algorithm.
    """
    N = size(LatticeConfiguration,1)
    Accepted = false

    # Extraction of a random site of the lattice
    x,y = GetRandomSite(N)

    # Reading of the current spin state of the site
    SiteSpin = LatticeConfiguration[x,y]

    """
    Calculation of total spins of neighbouring vertexes. Remember that in
    modular algebra you have e.g. mod(3,3)=0, and since the indexing of Julia
    starts by 1 we have to indicate the position '(x+1)%N,y' as 'mod(x,N)+1,y'.
    Instead of 'mod(x-2,N)+1,y' we use 'mod(x-2+N,N)+1,y' to indicate the
    position '(x-1)%N,y' because modular algebra gets in trouble with negative
    numbers.
    """
    NeighboursSpin = ( LatticeConfiguration[mod(x,N)+1, y] 
                     + LatticeConfiguration[mod(x-2+N,N)+1, y]
                     + LatticeConfiguration[x, mod(y,N)+1]
                     + LatticeConfiguration[x, mod(y-2+N,N)+1] )
    
    DeltaEnergy = 2.0*SiteSpin*NeighboursSpin

    if DeltaEnergy <= 0 || rand() < Acceptances[Int(DeltaEnergy/4)] 
        LatticeConfiguration[x,y] = -1*SiteSpin
        DeltaMagnetization = -2*SiteSpin
        Accepted = true

        # Remember to divide by volume in order to have energy and 
        # magnetization per unit volume
        DeltaEnergy /= (N^2)
        DeltaMagnetization /= (N^2)
    else 
        DeltaEnergy = 0
        DeltaMagnetization = 0
    end

    return Accepted, LatticeConfiguration
end
