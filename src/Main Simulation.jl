include("Parameters.jl")
include("Types and Structs.jl")
include("Simulation Functions.jl")
include("Metric Functions.jl")
include("Metrics and Recorder.jl")

using ProgressMeter

# --- Main Simulation Function --- 

function simulate!(counts, config::ModelConfig, metrics)
    """ 
    This function is the heart of the simulation. 
    It combines and calls the respective functions from the model configuration and records the designated metrics. 

    Parameters: 
        counts (array): The initial R_50 distribution for each species given as a (number of species) x (number of bins) x (number of wells) array
        config (ModelConfig): The model configuration given as a ModelConfig Struct (see "Types and Structs.jl")
        metrics (array): The metris that should be recorded. Can be passed from the config as config.metrics or passed separately. 

    Returns: 
        named tuple: containing the recorded metrics and the extinction time of the focal population in hours. 
    
    """

    p = config.params    # Load the parameters from the config struct 
    
    # Initialise the state of the system (see "Types and Structs.jl")
    state = SimState(
            copy(counts), 
            0.0,
            1,
            0.0, 
            p,
            nothing
        )

    N_other = 0    # Initialise the population size of the other strain 

    newborn_counts = zeros(Int, size(state.counts))    # Create array to store the newborn cells 

    @showprogress for t in 1:p.timesteps    # Loop over all timesteps given in p

        # Update and record the extinction time once the focal population dies out 
        if isnothing(state.extinction_time) && 
            sum(state.counts[:, 1, :]) == 0 
            
            state.extinction_time = state.time

            for metric in config.metrics
                config.record_fn(metric, state::SimState, t, p)
            end

            break
        end

        # For each well of the system apply the bin transformations (death, birth, mutation).
        for w in 1:p.no_wells
            step_well!(state, w, t, config, p, newborn_counts, N_other)
        end

        # Dilute the system if t equals the dilution interval. 
        if is_bottleneck(t, p)
            dilute!(config.dilution::Dilution, state.counts, p)

            # Set the drug concentration for the next day 
            state.drug_conc = config.drug_schedule(t/60*p.Δt_phys)

            # Increase day count by 1 
            state.day += 1
        end

        # Record the metrics 
        for metric in config.metrics
            config.record_fn(metric, state::SimState, t, p)
        end
        state.time += p.Δt_phys/60.0    # Update the state time
    end

    # After looping, if the focal population did not die out, set the extinction time to the duration of the simulation. 
    if isnothing(state.extinction_time)
        state.extinction_time = 24 * p.sim_days 
    end

    # Return the recorded metrics and the extinction time
    return merge(metrics, (; 
    extinction_time = state.extinction_time))
end

function run_simulation(counts0::Array{Int, 3}, config::ModelConfig)
    """
    Initialises the metrics and runs the simulation in one function. 
    This is important for running batch simulations as the metrics have to initialised between runs. 

    Parameters: 
        counts0 (array): The initial R_50 distribution for each species given as a (number of species) x (number of bins) x (number of wells) array
        config (ModelConfig): The model configuration given as a ModelConfig Struct (see "Types and Structs.jl")

    Returns: 
        named tuple: containing the recorded metrics and the extinction time of the focal population in hours. 
    """

    run_metrics = collect_metrics(config.metrics)
    init_metrics!(run_metrics)

    results = simulate!(counts0, config, run_metrics)

    return results
end







