using Random, Distributions, LinearAlgebra

# --- Simulation Functions --- 

function logistic(x, cc_max, x_0, k)
    """
    Logistic function used for the drug schedule.

    Parameters:
    - x       : time (hours)
    - cc_max  : maximum concentration
    - x_0     : midpoint (inflection point)
    - k       : steepness/base of exponential

    Returns:
    - Drug concentration at time x
    """
    cc_max / (1 + k^(x/24.0 - x_0))
end

function mutate_newborns!(counts, p, d_eff::Vector{Float64}, b_eff::Float64)
    """
    Applies mutation process to a vector (!) of newborn counts.

    Process:
    1. Determine how many mutate
    2. Apply survival probability to mutants
    3. Redistribute surviving mutants across bins via mutation kernel

    Parameters:
    - counts : Vector of newborn counts per bin (modified conceptually, not in-place)
    - p      : parameter struct (contains mutation_rate, nbins, K_mut, etc.)
    - d_eff  : effective death rates per bin
    - b_eff  : effective birth rate

    Returns:
    - new_counts : updated counts after mutation
    """
    
    new_counts = zeros(Int, p.nbins) 

    for i in 1:p.nbins
        n = counts[i]
        if n == 0
            continue # skip empty bins for efficiency
        end

        # --- Mutation step ---
        # Number of individuals that mutate
        n_mut = rand(Binomial(n, p.mutation_rate))

        # effective rates for this bin

        # --- Survival of mutants ---
        # Estimate probability mutant survives (fitness-based)
        p_est = clamp(1 - d_eff[i] / b_eff, 0.0, 1.0)

        n_est = rand(Binomial(n_mut, p_est))

        # --- Non-mutated survivors remain in same bin ---
        n_surv = n - n_mut
        new_counts[i] += n_surv

        # --- Redistribute surviving mutants ---
        if n_est > 0
            # Multinomial redistribution using mutation kernel
            m = rand(Multinomial(n_est, view(p.K_mut, i, :)))
            new_counts .+= m
        end
    end

    return new_counts
end

function dilute!(::NoPooling, counts, p)
    """
    Dilution without pooling:
    Each well is diluted independently.

    For each bin and species:
    - Apply Binomial thinning with dilution_factor
    """
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
    """
    Dilution with pooling:
    1. Pool all wells for each species
    2. Apply dilution
    3. Redistribute across wells uniformly
    """
    
    # Temporary storage for pooled counts
    pooled_counts = zeros(p.nbins, p.no_species) 

    # Treat each species seperately 
    for s in 1:p.no_species

        # Pool the counts of the species from all wells into one array
        @views pooled_counts[:,s] .= sum(counts[:,s,:], dims=2)

        # Distribute a random sampling from the pooled counts into each well again
        for i in 1:p.nbins

            n = pooled_counts[i,s]

            # --- Apply dilution ---
            n_surv = rand(Binomial(n, p.dilution_factor))

            # --- Redistribute survivors ---
            if n_surv > 0
                # Uniform redistribution across wells
                m = rand(Multinomial(n_surv, fill(1/p.no_wells, p.no_wells)))
                counts[i,s,:] .= m
            else
                counts[i,s,:] .= 0
            end
        end
    end
end

function drug_death(ic50::Float64, C::Float64; dmax::Float64=0.02, h::Float64=2.0)
    """
    Drug-induced death rate using Hill function.

    Parameters:
    - ic50 : concentration where effect is half-maximal
    - C    : current drug concentration
    - dmax : maximum death rate
    - h    : Hill coefficient (steepness)

    Returns:
    - death probability contribution
    """
    return dmax * (C^h / (C^h + ic50^h))
end

function crowd_growth(N; K = 1e8, h = 2.0)
    """
    Density-dependent growth suppression (Hill-type).

    Parameters:
    - N : population size
    - K : carrying capacity
    - h : steepness

    Returns:
    - scaling factor ∈ (0,1]
    """
    return 1/(1 + (N / K)^h)
end

function is_bottleneck(t, p)
    """
    Checks if current timestep is a bottleneck event.
    """
    return t % p.bottleneck_interval == 0 
end

growth_interaction(::NoInteraction, N_i, N_other, p) = 1.0 

growth_interaction(::Interaction, N_i, N_other, p) = 
    # Interaction function returning the interaction factor ∈ [α_int,1]
    # eps() prevents division by zero
    (N_i + p.α_int * (N_other * p.K_int)) / (N_i + (N_other * p.K_int) + eps())

growth_rate(i, s, w, counts, config, p) =
    # Growth scaled by crowding
    p.µ[s] * p.Δt * config.crowding_fn(sum(counts[:, s, w]); K = p.K) 

death_rate(d_drug, d0, Δt) = 
    # Calculate the effective death rate as base death + drug death + interaction death
    clamp((d0 + d_drug) .* Δt, 0.0, 1.0)

other_population(::NoInteraction, state, s, w) = 0

function other_population(::Interaction, state, s, w)
    """
    Returns population of competing species in same well.
    Assumes 2-species system (index trick: 3 - s).
    """
    j = 3 - s
    return sum(state.counts[:, j, w])
end

# --- Main Growth-Death-Mutate-Function 

function step_well!(state, w, t, config, p, newborn_counts, N_other)

    # Array containing the drug death rate for each R_50 value
    d_drug = drug_death.(p.bins, state.drug_conc)

    for s in 1:p.no_species

        # Get competing population
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

        # Re-initialize newborns  
        fill!(newborn_counts, 0)
    end
end







