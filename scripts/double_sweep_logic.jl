using DrWatson
@quickactivate "Cell_SE_Model"
using TOML # Ensure TOML is loaded for this script
using DataFrames
using CSV

include(srcdir("simulation_logic.jl"))
include(srcdir("cell_model.jl"))
include(srcdir("run_analysis.jl"))

struct DoubleSweepParameters
    sweep_parameter_name_one::String
    sweep_init_one::Float64
    sweep_step_one::Float64
    sweep_stop_one::Float64
    sweep_parameter_name_two::String
    sweep_init_two::Float64
    sweep_step_two::Float64
    sweep_stop_two::Float64
    sweep_repetitions::Int
end

function read_sweep_parameters(filename::String)
    parameters = TOML.parsefile(filename)
    return DoubleSweepParameters(
        String(_toml_value(parameters, "sweep_parameter_name_one")),
        Float64(_toml_value(parameters, "sweep_init_one")),
        Float64(_toml_value(parameters, "sweep_step_one")),
        Float64(_toml_value(parameters, "sweep_stop_one")),
        String(_toml_value(parameters, "sweep_parameter_name_two")),
        Float64(_toml_value(parameters, "sweep_init_two")),
        Float64(_toml_value(parameters, "sweep_step_two")),
        Float64(_toml_value(parameters, "sweep_stop_two")),
        Int(_toml_value(parameters, "sweep_repetitions"))
    )
end

# Optimization: Pass the dictionary directly so we only read the file once
function create_swept_parameters(base_dict::Dict, sweep_param_one::String, new_val_one, sweep_param_two::String, new_val_two)
    # Use deepcopy so we don't accidentally mutate the base template
    parameters_dict = deepcopy(base_dict)
    
    # Overwrite the specific parameter dynamically
    parameters_dict[sweep_param_one] = new_val_one
    parameters_dict[sweep_param_two] = new_val_two
    
    # Rebuild ALL parameter structures, not just SimulationParameters

    simulation_parameters = SimulationParameters(
        Int(_toml_value(parameters_dict, "num_cells")),
        Float64(_toml_value(parameters_dict, "dt")),
        Float64(_toml_value(parameters_dict, "total_time")),
        Float64(_toml_value(parameters_dict, "D_angle")),
        Float64(_toml_value(parameters_dict, "D_position")),
        Float64(_toml_value(parameters_dict, "theta_alignment")),
        Float64(_toml_value(parameters_dict, "run_to_tumble_rate_min")),
        Float64(_toml_value(parameters_dict, "run_to_tumble_rate_max")),
        String(_toml_value(parameters_dict, "run_to_tumble_rate_function")),
        Float64(_toml_value(parameters_dict, "run_to_tumble_rate_allure")),
        Float64(_toml_value(parameters_dict, "tumble_to_run_rate_min")),
        Float64(_toml_value(parameters_dict, "tumble_to_run_rate_max")),
        String(_toml_value(parameters_dict, "tumble_to_run_rate_function")),
        Float64(_toml_value(parameters_dict, "tumble_to_run_rate_allure")),        
        Float64(_toml_value(parameters_dict, "cell_hard_radius")),
        Float64(_toml_value(parameters_dict, "cell_soft_radius")),
        Float64(_toml_value(parameters_dict, "maximal_stretch")),
        Float64(_toml_value(parameters_dict, "run_speed_min")),
        Float64(_toml_value(parameters_dict, "run_speed_max")),
        String(_toml_value(parameters_dict, "run_speed_distribution")),
        Float64(_toml_value(parameters_dict, "top_internal_pull_strength")),
        Float64(_toml_value(parameters_dict, "soft_boundary_potential_strength")),
        Float64(_toml_value(parameters_dict, "bottom_internal_pull_strength")),
        Float64(_toml_value(parameters_dict, "cil_intensity"))
    )

    domain_specs = DomainSpecs(
        Float64(_toml_value(parameters_dict, "domain_width")),
        Float64(_toml_value(parameters_dict, "domain_height")),
    )

    morse_potential_top = MorsePotential(
        Float64(_toml_value(parameters_dict, "top_morse_u0")),
        Float64(_toml_value(parameters_dict, "top_morse_v0")),
        Float64(_toml_value(parameters_dict, "top_morse_xi1")),
        Float64(_toml_value(parameters_dict, "top_morse_xi2")),
        Float64(_toml_value(parameters_dict, "top_morse_cutoff")),
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
sweep_parameters = read_sweep_parameters(datadir("sims", "double_sweep_parameters.toml"))

# 3. Initialize a dictionary to store results
sweep_results = Dict{String, Vector{Float64}}()
sweep_results["p_value_one"] = Float64[]
sweep_results["p_value_two"] = Float64[]

# 4. Execute the sweep
for p_value_one in sweep_parameters.sweep_init_one:sweep_parameters.sweep_step_one:sweep_parameters.sweep_stop_one
    println("Running simulation with $(sweep_parameters.sweep_parameter_name_one) = $p_value_one")
    for p_value_two in sweep_parameters.sweep_init_two:sweep_parameters.sweep_step_two:sweep_parameters.sweep_stop_two
        println("Running simulation with $(sweep_parameters.sweep_parameter_name_two) = $p_value_two")
        analysis_results_list = AnalysisResults[]
        for rep in 1:sweep_parameters.sweep_repetitions
            println("Repetition $rep of $(sweep_parameters.sweep_repetitions)")
            
            # Rebuild all 4 structs with the updated dictionary
            swept_p, swept_domain, swept_top, swept_bottom = create_swept_parameters(
                base_toml_dict, 
                sweep_parameters.sweep_parameter_name_one, 
                p_value_one, sweep_parameters.sweep_parameter_name_two, 
                p_value_two
            )
            
            # Pass all 4 updated structs to the simulation loop
            solution, final_p, final_domain = simulation_loop(swept_p, swept_domain, swept_top, swept_bottom)
            
            # Run analysis
            analysis_results = run_analysis_from_solution(solution, final_p)
            
            # Store results locally for averaging later
            push!(analysis_results_list, analysis_results)
        end

    # Store p_value
    push!(sweep_results["p_value_one"], p_value_one)
    push!(sweep_results["p_value_two"], p_value_two)
    
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
    
    push!(sweep_results["mean_speeds_bottom"], mean(vcat([r.mean_speeds_bottom for r in analysis_results_list]...)))
    push!(sweep_results["mean_speeds_top"], mean(vcat([r.mean_speeds_top for r in analysis_results_list]...)))
    push!(sweep_results["group_speeds_bottom"], mean(vcat([r.group_speeds_bottom for r in analysis_results_list]...)))
    push!(sweep_results["group_speeds_top"], mean(vcat([r.group_speeds_top for r in analysis_results_list]...)))
    push!(sweep_results["mean_theta"], mean(vcat([r.mean_theta for r in analysis_results_list]...)))
end
end

# 5. Convert results to DataFrame and save to CSV
results_df = DataFrame(sweep_results)
filename = "analysis_$(sweep_parameters.sweep_parameter_name_one)_$(sweep_parameters.sweep_parameter_name_two).csv"
filepath = datadir("sims", filename)
CSV.write(filepath, results_df)
println("Results saved to $filepath")