using DrWatson
include(srcdir("cell_model.jl"))

struct Solution
    bottom::Vector{Vector{SVector{2, Float64}}}
    top::Vector{Vector{SVector{2, Float64}}}
end

function simulation_init()
    p, dom = read_parameters(datadir("sims", "parameters.toml"))
    collective = initialize_cell_collective(p, dom)
    return collective, p, dom
end

function simulation_loop()
    collective, p, domain = simulation_init()
    num_steps = Int(p.total_time / p.dt)
    solution = Solution(Vector{Vector{SVector{2, Float64}}}(undef, num_steps), Vector{Vector{SVector{2, Float64}}}(undef, num_steps))
    for step in 1:num_steps
        compute_cell_cohesion_forces!(collective, p)
        compute_stochastic_forces!(collective, p)
        update_cell_collective!(collective, p)
        apply_hard_wall_boundary_conditions!(collective, domain, p)
        store_solution!(solution, collective, step)
    end
    return solution
end

function store_solution!(solution::Solution, collective::CellCollective, step::Int)
    solution.bottom[step] = collective.bottom
    solution.top[step] = collective.top
end

