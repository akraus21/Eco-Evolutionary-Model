using Random, Distributions, LinearAlgebra

# --- Simulation Functions --- 

function logistic(x, cc_max, x_0, k)
    """
    Logistic function used for the drug schedule. 
    """
    cc_max / (1 + k^(x/24.0 - x_0))
end

function mutate_newborns!(counts, p, d_eff::Vector{Float64}, b_eff::Float64)
    """
    Applies the precomputed mutation kernel to a Vector (!) of Counts
    and mutates and distributes them to another array. 
    """
    
    new_counts = zeros(Int, p.nbins) 

    for i in 1:p.nbins
        n = counts[i]
        if n == 0
            continue
        end

        # how many mutate?
        n_mut = rand(Binomial(n, p.mutation_rate))

        # effective rates for this bin

        p_est = clamp(1 - d_eff[i] / b_eff, 0.0, 1.0)

        n_est = rand(Binomial(n_mut, p_est))

        # survivors
        n_surv = n - n_mut
        new_counts[i] += n_surv

        # distribute mutants
        if n_est > 0
            m = rand(Multinomial(n_est, view(p.K_mut, i, :)))
            new_counts .+= m
        end
    end

    return new_counts
end

function dilute!(::NoPooling, counts, p)
    for i in 1:p.nbins
        for s in 1:p.no_species
            n = counts[i, s, 1]
            # Expected survivors = n * dilution_factor
            # Draw from a Binomial to get stochastic survivors
            counts[i, s, 1] = rand(Binomial(round(Int,n), p.dilution_factor))
        end
    end
end

function dilute!(::Pooling, counts, p)
    # Initialise an array to store the pooled counts
    pooled_counts = zeros(p.nbins, p.no_species) ### Replace with preallocated array!

    # Treat each species seperately 
    for s in 1:p.no_species

        # Pool the counts of the species from all wells into one array
        @views pooled_counts[:,s] .= sum(counts[:,s,:], dims=2)

        # Distribute a random sampling from the pooled counts into each well again
        for i in 1:p.nbins

            n = pooled_counts[i,s]

            # dilution first
            n_surv = rand(Binomial(n, p.dilution_factor))

            # redistribute across wells
            if n_surv > 0
                m = rand(Multinomial(n_surv, fill(1/p.no_wells, p.no_wells)))
                counts[i,s,:] .= m
            else
                counts[i,s,:] .= 0
            end
        end
    end
end

"""
Logistic function takes in an IC50 value and a drug concentration C and returns a death rate. 
    Maximum of the rate is dmax and steepness is h. 
"""
function drug_death(ic50::Float64, C::Float64; dmax::Float64=0.02, h::Float64=2.0)
    return dmax * (C^h / (C^h + ic50^h))
end

"""
Function quantifying the effect of crowding on growth. 
The function is a Hill-type function, steepness is given by h. 
N: the whole population size of a species
K: carrying capacity in absolute numbers. 
Returns a number between 0.0 and 1.0.  
"""
function crowd_growth(N; K = 1e8, h = 2.0)
    return 1/(1 + (N / K)^h)
end

function is_bottleneck(t, p)
    return t % p.bottleneck_interval == 0 
end

growth_interaction(::NoInteraction, N_i, N_other, p) = 1.0 

growth_interaction(::Interaction, N_i, N_other, p) = 
    (N_i + p.α_int * (N_other * p.K_int)) / (N_i + (N_other * p.K_int) + eps())

growth_rate(i, s, w, counts, config, p) =
    p.µ[s] * p.Δt * config.crowding_fn(sum(counts[:, s, w]); K = p.K) 

death_rate(d_drug, d0, Δt) = 
    # Calculate the effective death rate as base death + drug death + interaction death
    clamp((d0 + d_drug) .* Δt, 0.0, 1.0)

other_population(::NoInteraction, state, s, w) = 0

function other_population(::Interaction, state, s, w)
    j = 3 - s
    return sum(state.counts[:, j, w])
end

# --- Main Growth-Death-Mutate-Function 

function step_well!(state, w, t, config, p, newborn_counts, N_other)

    # Array containing the drug death rate for each IC50 value
    d_drug = drug_death.(p.bins, state.drug_conc)

    for s in 1:p.no_species

        N_other = other_population(config.interaction, state, s, w)

        # 1. Death and Birth 

        for i in eachindex(p.bins)

            if state.counts[i, s, w] == 0
                continue 
            end

            # Draw cells from bin that will die 
            deaths = rand(Binomial(state.counts[i, s, w], config.death_fn(d_drug[i], p.d0, p.Δt_phys)))

            # Subtract cells from population 
            state.counts[i, s, w] -= deaths

            # Calculate effective birth count 
            λ = state.counts[i, s, w] * config.growth_fn(i, s, w, state.counts, config, p) * config.interaction_fn(config.interaction, state.counts[i, s, w], N_other, p)

            # Draw number of newborns using a negative binomial distribution 
            newborn_counts[i, s, w] = rand(NegativeBinomial(p.k, p.k / (p.k + λ)))
        end

        # 2. Apply mutation

        # mutate the newborns only
        mutated_newborns = config.mutation_fn(newborn_counts[:, s, w], p, config.death_fn.(d_drug, p.d0, p.Δt), p.µ[s] * p.Δt)

        # Add mutated newborns to counts (survivors already included)
        state.counts[:, s, w] .+= mutated_newborns

        # Re-initialize newborns and deaths 
        fill!(newborn_counts, 0)
    end
end







