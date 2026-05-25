# Re-run and animate one Pareto-optimal parameter set found by
# scripts/optimise_parameters.jl.
#
# Run:  julia --project=. scripts/inspect_optimised.jl
#
# Choose which Pareto set to inspect by editing SELECTION below:
#   "knee"             -> best penetration / cohesion ratio (balanced trade-off)
#   "best_penetration" -> deepest rightward penetration
#   "best_cohesion"    -> tightest cluster (smallest radius of gyration)
#   <an Integer>       -> the `sample` id of a specific row in optimisation_pareto.csv
#
# Outputs (in data/sims/):
#   parameters_optimised.toml — full, reusable config for the chosen set
#   optimised_animation.mp4   — animation of the chosen set

using DrWatson
@quickactivate "Cell_SE_Model"
using TOML
using DataFrames
using CSV
using Statistics

include(srcdir("simulation_logic.jl"))
include(scriptsdir("animation_logic.jl"))
include(srcdir("cell_model.jl"))
include(srcdir("cluster_metrics.jl"))

const SELECTION = "knee"

# ----------------------------------------------------------------- config
config = TOML.parsefile(datadir("sims", "optimisation_parameters.toml"))
base_file       = String(config["base_parameters_file"])
window_fraction = Float64(config["metric_window_fraction"])
range_keys      = sort(collect(keys(config["ranges"])))

pareto_path = datadir("sims", "optimisation_pareto.csv")
isfile(pareto_path) || error("$pareto_path not found — run scripts/optimise_parameters.jl first.")
pareto_df = CSV.read(pareto_path, DataFrame)
nrow(pareto_df) == 0 && error("optimisation_pareto.csv is empty — the search found no valid candidates.")

# ------------------------------------------------------------- pick a row
function choose_row(df, selection)
    if selection == "knee"
        return argmax(df.penetration ./ df.cohesion_rg)
    elseif selection == "best_penetration"
        return argmax(df.penetration)
    elseif selection == "best_cohesion"
        return argmin(df.cohesion_rg)
    elseif selection isa Integer
        idx = findfirst(==(selection), df.sample)
        idx === nothing && error("sample $selection not found in optimisation_pareto.csv")
        return idx
    else
        error("Unknown SELECTION: $selection")
    end
end

row = pareto_df[choose_row(pareto_df, SELECTION), :]
println("Inspecting Pareto sample $(Int(row.sample))  (SELECTION = $SELECTION)")
println("  recorded:  cohesion Rg = $(round(row.cohesion_rg, digits = 2))" *
        "   penetration = $(round(row.penetration, digits = 2))")

# ------------------------------------------------- rebuild the parameter set
base_dict = TOML.parsefile(datadir("sims", base_file))
d = deepcopy(base_dict)
for k in range_keys
    d[k] = row[k]
end
d["run_speed_max"] = d["run_speed_min"]   # run_speed_distribution = "constant"

optimised_path = datadir("sims", "parameters_optimised.toml")
open(optimised_path, "w") do io
    TOML.print(io, d)
end
println("Full optimised config written to $optimised_path")

# --------------------------------------------------------- run + animate
# `read_parameters` rebuilds all structs; `domain` becomes the module-global
# read by `compute_state_changes!`.
p, domain, morse_top, morse_bottom = read_parameters(optimised_path)
solution, p, domain = simulation_loop(p, domain, morse_top, morse_bottom)

coh = compute_cluster_cohesion(solution, p; window_fraction = window_fraction)
pen = compute_rightward_penetration(solution, p; window_fraction = window_fraction)
fro = compute_front_position(solution, p; window_fraction = window_fraction)
println("  re-run:    cohesion Rg = $(round(coh, digits = 2))" *
        "   penetration = $(round(pen, digits = 2))   front_x = $(round(fro, digits = 2))")

animation_path = datadir("sims", "optimised_animation.mp4")
animate_solution(solution, p, domain, animation_path)
println("Animation written to $animation_path")
