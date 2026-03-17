# --- Recorder and Metrics ---
"""
Defines the recorder functions and the metric structs that can be recorded. 
Each metric consists of four parts: 
    - The mutable constructor defining the objects where the recorded data is stored
    - A function returning the default state of the constructor 
    - A reset!() function resetting the constructor to its default state 
    - A record!() function recording the metrics of interest and storing them. 
The metrics must be passed as an array, containing each metric one wants to record as a tuple.
The tuple contains the metric constructor and the time interval at which the metric should be recorded. 
"""

function collect_metrics(metrics)
    """
    Creates a tuple out of the metrics given in an array, which makes it easier to reset them. 
    """
    NamedTuple{Tuple(Symbol(nameof(typeof(m.metric))) for m in metrics)}(
        Tuple(m.metric for m in metrics)
    )
end

function init_metrics!(metrics)
    """
    Loop over all metrics and initialise each using their respective reset!() function 
    """
    for m in eachindex(metrics)
        reset!(metrics[m])
    end
end

function record!(m::Every, state, t, p)
    """
    Record each metric using the respective record!() function and defined time interval. 
    """
    if t % m.Δt ≈ 0 
        record!(m.metric, state, p)
    end
end

# --- Population Size Metric ---
"""
A metric that records the population size of the focal strain and the time. 
"""

mutable struct PopSizeMetric <: Metric
    pop1::Vector{Int}
    time::Vector{Float64}
end

function PopSizeMetric()
    PopSizeMetric(
        Vector{Int}(undef, 0),
        Vector{Float64}(undef, 0)
    )
end

reset!(m::PopSizeMetric) = begin
    empty!(m.pop1)
    empty!(m.time)
end

function record!(m::PopSizeMetric, state::SimState, p)
    pop1 = sum(state.counts[:, 1, p.no_wells])

    push!(m.pop1, pop1)
    push!(m.time, state.time)
end

# --- Resistant Metric --- 

mutable struct ResistantMetric <: Metric 
    pop_focal::Vector{Int}
    pop_res::Vector{Int}
    shannon::Vector{Float64}
    ic50::Vector{Float64}
    time::Vector{Float64}
    drug_conc::Vector{Float64}
end

function ResistantMetric()
    ResistantMetric(
        Vector{Int}(undef, 0),
        Vector{Int}(undef, 0),
        Vector{Float64}(undef, 0),
        Vector{Float64}(undef, 0),
        Vector{Float64}(undef, 0),
        Vector{Float64}(undef, 0),
    )
end

reset!(m::ResistantMetric) = begin
    empty!(m.pop_focal)
    empty!(m.pop_res)
    empty!(m.shannon)
    empty!(m.ic50)
    empty!(m.time)
    empty!(m.drug_conc)
end

function record!(m::ResistantMetric, state::SimState, p)
    pop_focal = sum(state.counts[:, 1, p.no_wells])
    pop_res = sum(state.counts[p.bins .>= state.drug_conc, 1, p.no_wells])
    shannon = shannon_diversity(state.counts[:, 1, p.no_wells])
    ic50 = geometric_mean_ic50(p.bins, state.counts[:, 1, p.no_wells])

    push!(m.pop_focal, pop_focal)
    push!(m.pop_res, pop_res)
    push!(m.shannon, shannon)
    push!(m.ic50, ic50)
    push!(m.time, state.time)
    push!(m.drug_conc, state.drug_conc)
end

# --- Linearge Survival Metric ---

mutable struct LineageSurvivalMetric <: Metric 
    birth_time::Matrix{Union{Float64, Nothing}}
    death_time::Matrix{Union{Float64, Nothing}}
    alive::BitMatrix
end

function LineageSurvivalMetric(p::Params)
    LineageSurvivalMetric(
        fill(nothing, p.nbins, p.no_wells),
        fill(nothing, p.nbins, p.no_wells),
        falses(p.nbins, p.no_wells)
    )
end

function reset!(m::LineageSurvivalMetric)
    fill!(m.birth_time, nothing)
    fill!(m.death_time, nothing)
    fill!(m.alive, false)
end

function record!(m::LineageSurvivalMetric, state::SimState, p::Params)
    t = state.time
    for w in 1:p.no_wells
    counts = state.counts[:, 1, w] 
        for i in eachindex(counts)
            n = counts[i]
            if n > 0 && !m.alive[i, w] && isnothing(m.birth_time[i, w])
                m.alive[i, w] = true
                m.birth_time[i, w] = t
            end
            
            if n == 0 && m.alive[i, w] && isnothing(m.death_time[i, w]) 
                m.alive[i, w] = false
                m.death_time[i, w] = t
            end
        end
    end
end

function lin_survival_times(m::LineageSurvivalMetric)
    τ = Float64[]
    for i in eachindex(m.birth_time)
        if m.birth_time[i] !== nothing && m.death_time[i] !== nothing 
            push!(τ, m.death_time[i] - m.birth_time[i])
        end
    end
    return τ
end

p_est(m::LineageSurvivalMetric) =
(length(findall(lin_survival_times(m) .> 24.0)) + sum(m.alive)) / 
length(lin_survival_times(m)) 
