# Random-search optimiser for cluster cohesion vs. rightward penetration.
#
# Samples random parameter sets within the ranges in
# data/sims/optimisation_parameters.toml, runs the simulation for each,
# scores them with the metrics in src/cluster_metrics.jl, and reports the
# Pareto front of non-dominated (tight cluster, deep penetration) sets.
#
# Run:  julia --project=. scripts/optimise_parameters.jl
#
# Outputs (in data/sims/):
#   optimisation_results.csv  — every candidate, with its parameters + metrics
#   optimisation_pareto.csv   — only the Pareto-optimal candidates

using DrWatson
@quickactivate "Cell_SE_Model"
using TOML
using DataFrames
using CSV
using Random
using Statistics

include(srcdir("simulation_logic.jl"))
include(srcdir("cell_model.jl"))
include(srcdir("cluster_metrics.jl"))

# ----------------------------------------------------------------- config
config = TOML.parsefile(datadir("sims", "optimisation_parameters.toml"))

base_file       = String(config["base_parameters_file"])
num_samples     = Int(config["num_samples"])
repetitions     = Int(config["repetitions"])
random_seed     = Int(config["random_seed"])
window_fraction = Float64(config["metric_window_fraction"])
ranges          = config["ranges"]

# Sort the keys so sampling order — and therefore the seeded RNG — is
# deterministic (Dict key iteration order is not guaranteed).
range_keys = sort(collect(keys(ranges)))

base_dict = TOML.parsefile(datadir("sims", base_file))

# `compute_state_changes!` reads a module-global `domain`. The domain is
# fixed for every candidate, so bind it once here from the base config.
domain = read_parameters(datadir("sims", base_file))[2]

Random.seed!(random_seed)

tmp_file = datadir("sims", "_opt_tmp.toml")

"Reorder so that d[lo_key] <= d[hi_key]."
function order_pair!(d, lo_key, hi_key)
    if d[lo_key] > d[hi_key]
        d[lo_key], d[hi_key] = d[hi_key], d[lo_key]
    end
end

# ------------------------------------------------------------- run search
samples     = Int[]
cohesion    = Float64[]
penetration = Float64[]
front_x     = Float64[]
sampled     = Dict(k => Float64[] for k in range_keys)

println("Random search: $num_samples candidates x $repetitions repetition(s) on $base_file")
t0 = time()

for s in 1:num_samples
    d = deepcopy(base_dict)
    for k in range_keys
        lo, hi = Float64(ranges[k][1]), Float64(ranges[k][2])
        d[k] = lo + rand() * (hi - lo)
    end
    # Constraint fix-ups so the candidate is physically well-formed.
    order_pair!(d, "run_to_tumble_rate_min", "run_to_tumble_rate_max")
    order_pair!(d, "tumble_to_run_rate_min", "tumble_to_run_rate_max")
    d["run_speed_max"] = d["run_speed_min"]   # run_speed_distribution = "constant"

    open(tmp_file, "w") do io
        TOML.print(io, d)
    end

    coh_reps = Float64[]
    pen_reps = Float64[]
    fro_reps = Float64[]
    try
        p, dom, mtop, mbot = read_parameters(tmp_file)
        for _ in 1:repetitions
            solution, p_used, _ = simulation_loop(p, dom, mtop, mbot)
            push!(coh_reps, compute_cluster_cohesion(solution, p_used; window_fraction = window_fraction))
            push!(pen_reps, compute_rightward_penetration(solution, p_used; window_fraction = window_fraction))
            push!(fro_reps, compute_front_position(solution, p_used; window_fraction = window_fraction))
        end
    catch err
        @warn "Candidate $s failed; recording NaN." exception = err
        push!(coh_reps, NaN); push!(pen_reps, NaN); push!(fro_reps, NaN)
    end

    push!(samples, s)
    push!(cohesion, mean(coh_reps))
    push!(penetration, mean(pen_reps))
    push!(front_x, mean(fro_reps))
    for k in range_keys
        push!(sampled[k], d[k])
    end

    println("  [$s/$num_samples]  cohesion(Rg) = $(round(cohesion[end], digits = 2))" *
            "   penetration = $(round(penetration[end], digits = 2))")
end

isfile(tmp_file) && rm(tmp_file)
println("Search finished in $(round(time() - t0, digits = 1)) s")

# --------------------------------------------------------------- results
pareto = pareto_mask(cohesion, penetration)

df = DataFrame()
df.sample = samples
for k in range_keys
    df[!, k] = sampled[k]
end
df.cohesion_rg = cohesion
df.penetration = penetration
df.front_x     = front_x
df.pareto      = pareto

results_path = datadir("sims", "optimisation_results.csv")
CSV.write(results_path, df)
println("\nAll $num_samples runs written to $results_path")

pareto_df = sort(df[df.pareto, :], :penetration, rev = true)
pareto_path = datadir("sims", "optimisation_pareto.csv")
CSV.write(pareto_path, pareto_df)
println("Pareto front ($(nrow(pareto_df)) sets) written to $pareto_path")

if nrow(pareto_df) == 0
    @warn "No Pareto-optimal candidates — every simulation produced a non-finite metric. " *
          "Try narrowing the ranges in optimisation_parameters.toml."
else
    println("\nPareto-optimal trade-offs (tighter cluster = lower Rg, deeper = higher penetration):")
    for row in eachrow(sort(pareto_df, :cohesion_rg))
        println("  sample $(Int(row.sample)):  Rg = $(round(row.cohesion_rg, digits = 2))" *
                "   penetration = $(round(row.penetration, digits = 2))" *
                "   front_x = $(round(row.front_x, digits = 2))")
    end
    println("\nInspect / animate any of these with scripts/inspect_optimised.jl")
end
