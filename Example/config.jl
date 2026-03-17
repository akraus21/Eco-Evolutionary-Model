include(raw"C:\Users\SurfacePro8\Documents\Studium\Master Thesis\Theoretical Model\src\Model 2.0\Main Simulation.jl")

# --- Mono-Culture Model ---

mono_p = const_params(no_species = 1)

mono_config = ModelConfig(
    growth_fn = growth_rate,
    death_fn = death_rate,
    interaction = NoInteraction(),
    interaction_fn = growth_interaction,
    dilution = NoPooling(),
    mutation_fn = mutate_newborns!,
    crowding_fn = crowd_growth,
    drug_schedule = logistic_schedule,
    record_fn = record!,
    metrics = [Every(LineageSurvivalMetric(mono_p), 1.0), ], 
    params = mono_p
)

# --- Interaction Model --- 

int_co_p = const_params(no_species = 2)

int_co_config = ModelConfig(
    growth_fn = growth_rate,
    death_fn = death_rate,
    interaction = Interaction(),
    interaction_fn = growth_interaction,
    dilution = NoPooling(),
    mutation_fn = mutate_newborns!,
    crowding_fn = crowd_growth,
    drug_schedule = logistic_schedule,
    record_fn = record!,
    metrics = [Every(LineageSurvivalMetric(int_co_p), 1.0), ], 
    params = int_co_p
)

int_vol_p = const_params(no_species = 2, 
K = 96 * 1e8)

int_vol_config = ModelConfig(
    growth_fn = growth_rate,
    death_fn = death_rate,
    interaction = Interaction(),
    interaction_fn = growth_interaction,
    dilution = NoPooling(),
    mutation_fn = mutate_newborns!,
    crowding_fn = crowd_growth,
    drug_schedule = logistic_schedule,
    record_fn = record!,
    metrics = [Every(LineageSurvivalMetric(int_vol_p), 1.0), ], 
    params = int_vol_p
)

int_pool_p = const_params(no_species = 2, 
no_wells = 96)

int_pool_config = ModelConfig(
    growth_fn = growth_rate,
    death_fn = death_rate,
    interaction = Interaction(),
    interaction_fn = growth_interaction,
    dilution = Pooling(),
    mutation_fn = mutate_newborns!,
    crowding_fn = crowd_growth,
    drug_schedule = logistic_schedule,
    record_fn = record!,
    metrics = [Every(LineageSurvivalMetric(int_pool_p), 1.0), ], 
    params = int_pool_p
)

# --- Population Size Model ---

pop_co_p = const_params(no_species = 2, 
K = 1e7)

pop_co_config = ModelConfig(
    growth_fn = growth_rate,
    death_fn = death_rate,
    interaction = NoInteraction(),
    interaction_fn = growth_interaction,
    dilution = NoPooling(),
    mutation_fn = mutate_newborns!,
    crowding_fn = crowd_growth,
    drug_schedule = logistic_schedule,
    record_fn = record!,
    metrics = [Every(LineageSurvivalMetric(pop_co_p), 1.0), ], 
    params = pop_co_p
)

pop_vol_p = const_params(no_species = 2,
K = 96 * 1e7)

pop_vol_config = ModelConfig(
    growth_fn = growth_rate,
    death_fn = death_rate,
    interaction = NoInteraction(),
    interaction_fn = growth_interaction,
    dilution = NoPooling(),
    mutation_fn = mutate_newborns!,
    crowding_fn = crowd_growth,
    drug_schedule = logistic_schedule,
    record_fn = record!,
    metrics = [Every(LineageSurvivalMetric(pop_vol_p), 1.0), ], 
    params = pop_vol_p
)

pop_pool_p = const_params(no_species = 2, 
no_wells = 96,
K = 1e7)

pop_pool_config = ModelConfig(
    growth_fn = growth_rate,
    death_fn = death_rate,
    interaction = NoInteraction(),
    interaction_fn = growth_interaction,
    dilution = Pooling(),
    mutation_fn = mutate_newborns!,
    crowding_fn = crowd_growth,
    drug_schedule = logistic_schedule,
    record_fn = record!,
    metrics = [Every(LineageSurvivalMetric(pop_pool_p), 1.0), ], 
    params = pop_pool_p
)

models = Dict(
    :mono => (counts0 = mono_counts0, config = mono_config, repls = 384),
    :int_co => (counts0 = co_counts0, config = int_co_config, repls = 384), 
    :int_vol => (counts0 = co_counts0, config = int_vol_config, repls = 384),
    :int_pool => (counts0 = pool_counts0, config = int_pool_config, repls = 100),
    :pop_co => (counts0 = co_counts0, config = pop_co_config, repls = 384),
    :pop_vol => (counts0 = co_counts0, config = pop_vol_config, repls = 384),
    :pop_pool => (counts0 = pool_counts0, config = pop_pool_config, repls = 100)
)