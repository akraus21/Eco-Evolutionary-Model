# --- Types and Structs --- 

abstract type InteractionModel end
struct NoInteraction <: InteractionModel end
struct Interaction <: InteractionModel end

abstract type Dilution end
struct NoPooling <: Dilution end
struct Pooling <: Dilution end 

abstract type Metric end 

mutable struct SimState
    counts::Array{Int, 3}
    drug_conc::Float64
    day::Int 
    time::Float64
    params::Params
    extinction_time::Union{Nothing, Float64}
end

struct Every{M<:Metric} <: Metric
    metric::M
    Δt::Float64
end

Base.@kwdef struct ModelConfig
    growth_fn
    death_fn
    interaction::InteractionModel
    interaction_fn
    dilution::Dilution
    mutation_fn
    crowding_fn
    drug_schedule
    record_fn
    metrics::Vector
    params::Params
end
