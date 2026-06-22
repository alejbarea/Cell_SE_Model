using DrWatson
include(srcdir("cell_model.jl"))

struct Solution
    bottom::Vector{Vector{SVector{2, Float64}}}
    top::Vector{Vector{SVector{2, Float64}}}
    state::Vector{Vector{Int}}
    theta::Vector{Vector{Float64}}
end

function simulation_init(p = Nothing, dom = Nothing, morse_potential_top = Nothing, morse_potential_bottom = Nothing)
    if p === Nothing || dom === Nothing || morse_potential_top === Nothing || morse_potential_bottom === Nothing
        p, dom, morse_potential_top, morse_potential_bottom = read_parameters(datadir("sims", "parameters.toml"))
    end
    collective = initialize_cell_collective(p, dom)
    return collective, p, dom, morse_potential_top, morse_potential_bottom
end

function simulation_loop(p = Nothing, dom = Nothing, morse_potential_top = Nothing, morse_potential_bottom = Nothing)
    collective, p, domain, morse_potential_top, morse_potential_bottom = simulation_init(p, dom, morse_potential_top, morse_potential_bottom)
    num_steps = Int(p.total_time / p.dt)
    solution = Solution(Vector{Vector{SVector{2, Float64}}}(undef, num_steps), Vector{Vector{SVector{2, Float64}}}(undef, num_steps), Vector{Vector{Int}}(undef, num_steps), Vector{Vector{Float64}}(undef, num_steps))
    for step in 1:num_steps
        solution.bottom[step] = Vector{SVector{2, Float64}}(undef, p.num_cells)
        solution.top[step] = Vector{SVector{2, Float64}}(undef, p.num_cells)
        solution.state[step] = Vector{Int}(undef, p.num_cells)
        solution.theta[step] = Vector{Float64}(undef, p.num_cells)
        compute_hard_interaction_forces!(collective, p)
        compute_soft_interaction_forces!(collective, p, morse_potential_top, morse_potential_bottom)
        compute_adhesion_forces!(collective, p)
        compute_stochastic_forces!(collective, p)
        compute_state_changes!(collective,p)
        update_cell_collective!(collective, p)
        compute_cell_cohesion_forces!(collective, p)
        compute_soft_wall_forces!(collective,domain,p)
        apply_hard_wall_boundary_conditions!(collective, domain, p)
        store_solution!(solution, collective, step,p.num_cells)
    end
    return solution, p, domain
end

function store_solution!(solution::Solution, collective::CellCollective, step::Int, num_cells::Int)
    for i in 1:num_cells
            solution.bottom[step][i] = collective.bottom[i]
            solution.top[step][i] = collective.top[i]
            solution.state[step][i] = collective.state[i]
            solution.theta[step][i] = collective.theta[i]
    end
end

