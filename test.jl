using Logging
using LinearAlgebra

include("./simulations/code/functions.jl")

function sborrami(N)
    """ fà N stepz montecarlo, partendo da +1 """
    L = 20
    β = log(3)/4
    Acc = GetStepAcceptability(β)
    LatticeConfiguration = SetLattice(L)

    sborra = 0

    for i in 1:N
        a, LatticeConfiguration = NextMetropolisStep(L, LatticeConfiguration, β, Acc)
        if a
            sborra +=1
        end
        i += 1
    end

    return sborra, LatticeConfiguration
end

N = 100

n, Lattice = sborrami(N)

# println("dopo $N Metropolis, ho accettato $n volte \n")

# println("Ora mi sento così:")
# for i in 1:20
#     println(Lattice[i,:])
# end





function GrowSborra(L::Int64,
    LatticeConfiguration::Matrix{Float64},
    ExpansionProbability::Float64,
    Site::Tuple{Int64, Int64})
        
    SiteSpin = LatticeConfiguration[Site...]
    Neighbours = GetNeighboursTriangle(L, Site)
    LatticeConfiguration[Site...] *= -1 # flip

    #push!(Cluster, Site)

    for NeighbourSite in Neighbours
    if ( LatticeConfiguration[NeighbourSite...] == SiteSpin 
    && rand() < ExpansionProbability 
    # these conditions imply NeighbourSite ∉ Cluster already
    )
    GrowCluster(L, LatticeConfiguration, 
            ExpansionProbability, NeighbourSite)
    end
    end
end

function NextSborra(
            L::Int64,
            LatticeConfiguration::Matrix{Float64}, 
			ExpansionProbability::Float64)
    
    StartingSite = (3, 3)
    
    GrowSborra(L, LatticeConfiguration, ExpansionProbability, 
    			  StartingSite)
end

function PrintBello(Lattice)
for i in 1:L
    println(join([x == 1 ? "*" : " " for x in Lattice[i,:]], " "))
end
println("-------------------------------------------")
end

L = 5
β = log(3)/4
LatticeConfiguration = rand([-1.0, 1.0], L, L)
ExpansionProbability = 1 - exp(-2*β)

PrintBello(LatticeConfiguration)

NextSborra(L,LatticeConfiguration,ExpansionProbability)

PrintBello(LatticeConfiguration)

# for i in 1:L
#     println(join([x == 1 ? "*" : " " for x in OutputLattice[i,:]], " "))
#  end

# Clus = OutputLattice - LatticeConfiguration

# for i in 1:L
#     println(join([Clus[i,j] != 0 ? "c" : (LatticeConfiguration[i,j] == 1 ? "*" : " ") for j in 1:L], " "))
# end