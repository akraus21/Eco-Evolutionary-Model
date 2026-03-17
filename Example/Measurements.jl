include("config.jl")

using DataFrames, ProgressMeter

rows = NamedTuple[]

@showprogress for (name, model) in models
    nreps = model.repls
    Threads.@threads for r in 1:nreps

        results = run_simulation(model.counts0, model.config)
        row = (Model = name, Extinction_Time = results.extinction_time, P_est = p_est(results.LineageSurvivalMetric))

        push!(rows, merge(row, (; Replicate = r)))
    end
end

df = DataFrame(rows)

using JLD2
save_object(raw"C:\Users\SurfacePro8\Documents\Studium\Master Thesis\Theoretical Model\src\Model 2.0\Establishment Probability\Measurements.jld2", df)