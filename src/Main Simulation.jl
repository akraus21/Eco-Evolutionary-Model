include("Parameters.jl")
include("Types and Structs.jl")
include("Simulation Functions.jl")
include("Metric Functions.jl")
include("Metrics and Recorder.jl")

using ProgressMeter


# --- Main Simulation Function --- 

function simulate!(counts, config::ModelConfig, metrics)

    p = config.params

    state = SimState(
            copy(counts), 
            0.0,
            1,
            0.0, 
            p,
            nothing
        )

    N_other = 0

    # Create array to store the newborn cells 
    newborn_counts = zeros(Int, size(state.counts))

    @showprogress for t in 1:p.timesteps

        if isnothing(state.extinction_time) && 
            sum(state.counts[:, 1, :]) == 0 
            
            state.extinction_time = state.time

            for metric in config.metrics
                config.record_fn(metric, state::SimState, t, p)
            end

            break
        end

        for w in 1:p.no_wells
            step_well!(state, w, t, config, p, newborn_counts, N_other)
        end

        if is_bottleneck(t, p)
            dilute!(config.dilution::Dilution, state.counts, p)

            # Set the drug concentration for the next day 
            state.drug_conc = config.drug_schedule(t/60*p.Δt_phys)

            # Increase day count by 1 
            state.day += 1
        end

        for metric in config.metrics
            config.record_fn(metric, state::SimState, t, p)
        end
        state.time += p.Δt_phys/60.0
    end

    if isnothing(state.extinction_time)
        state.extinction_time = 24 * p.sim_days 
    end

    return merge(metrics, (; 
    extinction_time = state.extinction_time))
end

function run_simulation(counts0::Array{Int, 3}, config::ModelConfig)

    run_metrics = collect_metrics(config.metrics)
    init_metrics!(run_metrics)

    results = simulate!(counts0, config, run_metrics)

    return results
end







