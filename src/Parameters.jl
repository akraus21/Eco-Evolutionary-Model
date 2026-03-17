using Distributions

# --- Simulation Parameters --- 

"""
    Set the drug schedule, i.e. the time evolution of the drug concentration. 
"""
function logistic_schedule(t)
    return logistic(t, 100.0, 13.502, 0.746)
end

"""
Function that calculates a Mutation Kernel once, returning a nbins x nbins matrix 
where entries contain the probability to jump from bin_i to bin_j. 
"""

function mutation_kernel_discrete(bins; μ=0.0, σ=0.15)
    nb = length(bins)
    log_bins = log10.(bins)
    K = zeros(Float64, nb, nb)
    dist = Normal(μ, σ)

    for i in 1:nb
        Δ = log_bins .- log_bins[i]
        K[i, :] .= pdf.(dist, Δ)
        K[i, :] ./= sum(K[i, :])
    end
    return K
end

mutable struct Params
    """
    Defines the parameters constructor. 
    Contains all relevant parameters used in the simulation. 
    """
    
    # model structure
    no_species::Int 
    no_wells::Int 
    nbins::Int  # number of R50 bins

    # R50 bins
    R50_min::Float64 
    R50_max::Float64 

    # time
    Δt::Float64
    Δt_phys::Float64   # Physical time simulated in minutes 

    # evolution
    mutation_rate::Float64 # probability per timestep
    k::Float64 

    # growth
    µ_dt::Vector{Float64}  # species-specific doubling times in minutes
    α_int::Float64  # Minimum Growth Interaction
    K_int::Float64

    # death 
    d0::Float64  # Global death rate 

    # ecology
    K::Float64 
    dilution_factor::Float64  #Fraction of cells that survive dilution

    # experiment
    sim_days::Int  # Number of days being simulated

    # derived constants (cached)
    bins::Vector{Float64}
    µ::Vector{Float64}
    timesteps::Int
    bottleneck_interval::Int
    K_mut::Matrix{Float64}
end

function const_params(;
    no_species = 1,
    no_wells = 1,
    nbins = 200,
    R50_min = 1.0,
    R50_max = 100.0,
    Δt = 1.0,
    Δt_phys = 1.0,
    mutation_rate = 5e-6,
    k = 5.0,
    µ_dt = [52.0, 40.0],
    α_int = 0.35,
    K_int = 5e-7,
    d0 = 0.002,
    K = 1e8,
    dilution_factor = 1 / 300,
    sim_days = 15
)
    """
    Constructs a Params constructor from the input parameters.
    Additionally calculates the growth rates, timesteps, bottleneck interval and the mutation kernel. 
    Should be used to generate the Params constructor. 
    const_params() returns the default values. 
    """
    bins = 10 .^ range(log10(R50_min), log10(R50_max), length=nbins)
    µ = log(2) ./ µ_dt
    timesteps = round(Int, sim_days * 1440 / Δt_phys - 1)
    bottleneck_interval = round(Int, 1440 / Δt_phys)
    K_mut = mutation_kernel_discrete(bins; μ=0.0, σ=0.15)

    return Params(
        no_species, no_wells, nbins,
        R50_min, R50_max,
        Δt, Δt_phys,
        mutation_rate, k,
        µ_dt, α_int, K_int,
        d0,
        K, dilution_factor, 
        sim_days,
        bins, µ, timesteps, bottleneck_interval, K_mut
    )
end

# --- Default initial counts --- 
"""
The default initial counts for three conditions: 
    - monoculture 
    - co-culture in a single well 
    - co-culture in a metapopulation with 96 wells
"""

focal_init_R50 = 88 # Set initial R_50 of focal strain to ~7.5
co_strain_init_R50 = 200 # Set the initial R_50 of co-strain to max value 

# Initialize population: counts[bin, species]
mono_counts0 = zeros(Int, const_params().nbins, 1, 1)
mono_counts0[focal_init_R50, 1, 1] = 1e5

co_counts0 = zeros(Int, const_params().nbins, 2, 1)
co_counts0[focal_init_R50, 1, 1] = 1e5
co_counts0[200, 2, 1] = 1e5

pool_counts0 = zeros(Int, const_params().nbins, 2, 96)
pool_counts0[focal_init_R50, 1, :] .= 1e5
pool_counts0[200, 2, :] .= 1e5
