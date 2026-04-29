using DrWatson
include(srcdir("simulation_logic.jl"))




struct AnalysisResults
    mean_speeds_bottom::Vector{Float64}
    mean_speeds_top::Vector{Float64}
    group_speeds_bottom::Vector{Float64}
    group_speeds_top::Vector{Float64}
    mean_theta::Vector{Float64}
end


function run_analysis_from_simulation()
    solution, p, domain = simulation_loop()
    mean_speeds_bottom, mean_speeds_top = compute_mean_speed(solution, p)
    group_speeds_bottom, group_speeds_top = compute_group_speed(solution, p)
    mean_theta = compute_mean_theta(solution, p)

    return AnalysisResults(mean_speeds_bottom, mean_speeds_top, group_speeds_bottom, group_speeds_top, mean_theta)
end

function run_analysis_from_csv(solution_dir = datadir("sims"), solution_file = "solution.csv", parameters_file = "parameters.toml")
    solution = read_solution_from_csv(datadir(solution_dir, solution_file))
    p, dom = read_parameters(datadir(solution_dir, parameters_file))

    mean_speeds_bottom, mean_speeds_top = compute_mean_speed(solution, p)
    group_speeds_bottom, group_speeds_top = compute_group_speed(solution, p)
    mean_theta = compute_mean_theta(solution, p)

    return AnalysisResults(mean_speeds_bottom, mean_speeds_top, group_speeds_bottom, group_speeds_top, mean_theta)
end

function run_analysis_from_solution(solution::Solution, p::Parameters)
    mean_speeds_bottom, mean_speeds_top = compute_mean_speed(solution, p)
    group_speeds_bottom, group_speeds_top = compute_group_speed(solution, p)
    mean_theta = compute_mean_theta(solution, p)

    return AnalysisResults(mean_speeds_bottom, mean_speeds_top, group_speeds_bottom, group_speeds_top, mean_theta)
end



function compute_mean_speed(solution::Solution, p::Parameters)
    num_steps = length(solution.bottom)
    mean_speeds_bottom = zeros(num_steps)
    mean_speeds_top = zeros(num_steps)
    for step in 2:num_steps
        total_speed_bottom = 0.0
        total_speed_top = 0.0
        for i in 1:p.num_cells
            speed_bottom = norm(solution.bottom[step][i] - solution.bottom[step-1][i]) / p.dt
            speed_top = norm(solution.top[step][i] - solution.top[step-1][i]) / p.dt
            total_speed_bottom += speed_bottom
            total_speed_top += speed_top
        end
        mean_speeds_bottom[step] = total_speed_bottom / p.num_cells
        mean_speeds_top[step] = total_speed_top / p.num_cells
    end
    return mean_speeds_bottom, mean_speeds_top
end

function compute_group_speed(solution::Solution, p::Parameters)
    num_steps = length(solution.bottom)
    group_speeds_bottom = zeros(num_steps)
    group_speeds_top = zeros(num_steps)
    for step in 2:num_steps
        group_speed_bottom = norm(mean(solution.bottom[step]) - mean(solution.bottom[step-1])) / p.dt
        group_speed_top = norm(mean(solution.top[step]) - mean(solution.top[step-1])) / p.dt
        group_speeds_bottom[step] = group_speed_bottom
        group_speeds_top[step] = group_speed_top
    end
    return group_speeds_bottom, group_speeds_top
end

function compute_mean_theta(solution::Solution, p::Parameters)
    num_steps = length(solution.theta)
    mean_theta = zeros(num_steps)
    for step in 1:num_steps
        mean_theta[step] = mean(solution.theta[step])
    end
    return mean_theta
end