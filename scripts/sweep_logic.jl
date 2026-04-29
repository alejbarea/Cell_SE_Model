using DrWatson
@quickactivate "Cell_SE_Model"
using TOML # Ensure TOML is loaded for this script
using DataFrames
using CSV

include(srcdir("simulation_logic.jl"))
include(srcdir("cell_model.jl"))
include(srcdir("run_analysis.jl"))

struct SweepParameters
    sweep_parameter_name::String
    sweep_init::Float64
    sweep_step::Float64
    sweep_stop::Float64
end

function read_sweep_parameters(filename::String)
    parameters = TOML.parsefile(filename)
    return SweepParameters(
        String(_toml_value(parameters, "sweep_parameter_name")),
        Float64(_toml_value(parameters, "sweep_init")),
        Float64(_toml_value(parameters, "sweep_step")),
        Float64(_toml_value(parameters, "sweep_stop"))
    )
end

# Optimization: Pass the dictionary directly so we only read the file once
function create_swept_parameters(base_dict::Dict, sweep_param::String, new_val)
    # Use deepcopy so we don't accidentally mutate the base template
    parameters_dict = deepcopy(base_dict)
    
    # Overwrite the specific parameter dynamically
    parameters_dict[sweep_param] = new_val
    
    # Rebuild ALL parameter structures, not just SimulationParameters
    simulation_parameters = SimulationParameters(
        Int(_toml_value(parameters_dict, "num_cells")),
        Float64(_toml_value(parameters_dict, "dt")),
        Float64(_toml_value(parameters_dict, "total_time")),
        Float64(_toml_value(parameters_dict, "D_angle")),
        Float64(_toml_value(parameters_dict, "D_position")),
        Float64(_toml_value(parameters_dict, "theta_alignment")),
        Float64(_toml_value(parameters_dict, "run_to_tumble_rate")),
        Float64(_toml_value(parameters_dict, "tumble_to_run_rate")),
        Float64(_toml_value(parameters_dict, "cell_hard_radius")),
        Float64(_toml_value(parameters_dict, "cell_soft_radius")),
        Float64(_toml_value(parameters_dict, "maximal_stretch")),
        Float64(_toml_value(parameters_dict, "run_speed")),
        Float64(_toml_value(parameters_dict, "top_internal_pull_strength")),
        Float64(_toml_value(parameters_dict, "soft_boundary_potential_strength")),
        Float64(_toml_value(parameters_dict, "bottom_internal_pull_strength")),
        Float64(_toml_value(parameters_dict, "cil_intensity"))
    )

    domain_specs = DomainSpecs(
        Float64(_toml_value(parameters_dict, "domain_width")),
        Float64(_toml_value(parameters_dict, "domain_height"))
    )

    morse_potential_top = MorsePotential(
        Float64(_toml_value(parameters_dict, "top_morse_u0")),
        Float64(_toml_value(parameters_dict, "top_morse_v0")),
        Float64(_toml_value(parameters_dict, "top_morse_xi1")),
        Float64(_toml_value(parameters_dict, "top_morse_xi2")),
        Float64(_toml_value(parameters_dict, "top_morse_cutoff"))
    )

    morse_potential_bottom = MorsePotential(
        Float64(_toml_value(parameters_dict, "bottom_morse_u0")),
        Float64(_toml_value(parameters_dict, "bottom_morse_v0")),
        Float64(_toml_value(parameters_dict, "bottom_morse_xi1")),
        Float64(_toml_value(parameters_dict, "bottom_morse_xi2")),
        Inf
    )

    return simulation_parameters, domain_specs, morse_potential_top, morse_potential_bottom
end

# 1. Parse the base TOML into a dictionary ONCE outside the loop (much faster)
base_toml_dict = TOML.parsefile(datadir("sims", "parameters.toml"))

# 2. Read your sweep logic
sweep_parameters = read_sweep_parameters(datadir("sims", "sweep_parameters.toml"))

# 3. Initialize a dictionary to store results
sweep_results = Dict{String, Vector{Float64}}()
sweep_results["p_value"] = Float64[]

# 4. Execute the sweep
for p_value in sweep_parameters.sweep_init:sweep_parameters.sweep_step:sweep_parameters.sweep_stop
    println("Running simulation with $(sweep_parameters.sweep_parameter_name) = $p_value")
    
    # Rebuild all 4 structs with the updated dictionary
    swept_p, swept_domain, swept_top, swept_bottom = create_swept_parameters(
        base_toml_dict, 
        sweep_parameters.sweep_parameter_name, 
        p_value
    )
    
    # Pass all 4 updated structs to the simulation loop
    solution, final_p, final_domain = simulation_loop(swept_p, swept_domain, swept_top, swept_bottom)
    
    # Run analysis
    analysis_results = run_analysis_from_solution(solution, final_p)
    
    # Store p_value
    push!(sweep_results["p_value"], p_value)
    
    # Compute mean for each field in AnalysisResults
    if !haskey(sweep_results, "mean_speeds_bottom")
        sweep_results["mean_speeds_bottom"] = Float64[]
    end
    if !haskey(sweep_results, "mean_speeds_top")
        sweep_results["mean_speeds_top"] = Float64[]
    end
    if !haskey(sweep_results, "group_speeds_bottom")
        sweep_results["group_speeds_bottom"] = Float64[]
    end
    if !haskey(sweep_results, "group_speeds_top")
        sweep_results["group_speeds_top"] = Float64[]
    end
    if !haskey(sweep_results, "mean_theta")
        sweep_results["mean_theta"] = Float64[]
    end
    
    push!(sweep_results["mean_speeds_bottom"], mean(analysis_results.mean_speeds_bottom))
    push!(sweep_results["mean_speeds_top"], mean(analysis_results.mean_speeds_top))
    push!(sweep_results["group_speeds_bottom"], mean(analysis_results.group_speeds_bottom))
    push!(sweep_results["group_speeds_top"], mean(analysis_results.group_speeds_top))
    push!(sweep_results["mean_theta"], mean(analysis_results.mean_theta))
end

# 5. Convert results to DataFrame and save to CSV
results_df = DataFrame(sweep_results)
filename = "analysis_$(sweep_parameters.sweep_parameter_name).csv"
filepath = datadir("sims", filename)
CSV.write(filepath, results_df)
println("Results saved to $filepath")