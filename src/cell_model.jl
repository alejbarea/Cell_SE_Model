using StaticArrays, Distributions, TOML, LinearAlgebra

abstract type Parameters end

struct SimulationParameters <: Parameters
    num_cells::Int
    dt::Float64
    total_time::Float64
    D_angle::Float64
    D_position::Float64
    theta_alignment::Float64
    run_to_tumble_rate::Float64
    tumble_to_run_rate::Float64
    cell_hard_radius::Float64
    cell_soft_radius::Float64
    internal_cohesion_strength::Float64
    maximal_stretch::Float64
    run_speed::Float64
end

struct DomainSpecs <: Parameters
    domain_width::Float64
    domain_height::Float64
end


# Cell Collective Structure is a collection of all cell data for easy vectorization.
struct CellCollective
    bottom::Vector{SVector{2, Float64}} # Coordinates of cell centers at z = 0
    top::Vector{SVector{2, Float64}} # Coordinates of cell centers at z = z0
    theta::Vector{Float64} # Cell angles
    thetabar::Vector{Float64} # Cell desired polarity angles
    neighbours::Vector{Vector{Int}} # Neighbouring cell indices
    forcesbottom::Vector{SVector{2, Float64}} # Forces on each cell at z = 0
    forcestop::Vector{SVector{2, Float64}} # Forces on each cell at z = z0
    state::Vector{Int} # State of each cell (e.g., 0 for tumbling, 1 for running)
    state_timer::Vector{Float64} # Timer for state transitions
    run_speeds::Vector{Float64} # Speed of each cell when in the running state
end

function initialize_cell_collective(p::SimulationParameters, domain::DomainSpecs)
    bottom = Vector{SVector{2, Float64}}(undef, p.num_cells)
    top = Vector{SVector{2, Float64}}(undef, p.num_cells)
    theta = Vector{Float64}(undef, p.num_cells)
    thetabar = Vector{Float64}(undef, p.num_cells)
    neighbours = Vector{Vector{Int}}(undef, p.num_cells)
    forcesbottom = Vector{SVector{2, Float64}}(undef, p.num_cells)
    forcestop = Vector{SVector{2, Float64}}(undef, p.num_cells)
    state = Vector{Int}(undef, p.num_cells)
    state_timer = Vector{Float64}(undef, p.num_cells)
    run_speeds = Vector{Float64}(undef, p.num_cells)
    for i in 1:p.num_cells
        bottom[i] = SVector(rand() * domain.domain_width, rand() * domain.domain_height)
        top[i] = bottom[i] + SVector(0.0, p.cell_hard_radius) # Initial top position directly above bottom
        theta[i] = rand() * 2 * pi # Random initial angle
        thetabar[i] = theta[i] # Initial desired angle same as initial angle
        forcesbottom[i] = SVector(0.0, 0.0)
        forcestop[i] = SVector(0.0, 0.0)
        state[i] = rand(Bool) ? 1 : 0 # Random initial state
        state_timer[i] = rand(Exponential(1.0)) # Random initial timer for state transitions
        run_speeds[i] = p.run_speed # Random speed between 0 and p.run_speed for running state
    end
    return CellCollective(bottom, top, theta, thetabar, neighbours, forcesbottom, forcestop, state, state_timer, run_speeds)
end

function update_cell_collective!(collective::CellCollective, p::SimulationParameters)
    # Update bottom and top coordinates
    collective.bottom .+= p.dt .* collective.forcesbottom
    collective.top .+= p.dt .* collective.forcestop

    # Update theta only where state == 0
    mask = collective.state .== 0
    collective.theta[mask] .+= -p.theta_alignment .* sin.(collective.theta[mask] - collective.thetabar[mask]) .* p.dt + sqrt(2 * p.D_angle * p.dt) .* randn(sum(mask))
    collective.bottom[.!mask] .+= collective.run_speeds[.!mask] * p.dt .* SVector.(cos.(collective.theta[.!mask]), sin.(collective.theta[.!mask]))
end

function compute_stochastic_forces!(collective::CellCollective, p::SimulationParameters)
        for i in 1:p.num_cells
            collective.forcesbottom[i] += sqrt(2 * p.D_position * p.dt) .* randn(2)
            collective.forcestop[i] += sqrt(2 * p.D_position * p.dt) .* randn(2)
        end
        
end

function compute_state_changes!(collective::CellCollective, p::SimulationParameters)
    collective.state_timer .-= p.dt*((1 .- collective.state) .* p.tumble_to_run_rate .+ collective.state .* p.run_to_tumble_rate)
    # Check for state transitions
    mask = collective.state_timer .<= 0
    collective.state[mask] .= 1 .- collective.state[mask] # Toggle state
    # Reset timers for cells that changed state
    exponential_dist = Exponential(1.0)
    collective.state_timer[mask] .= exponential_dist.(sum(mask))
end

function apply_periodic_boundary_conditions!(collective::CellCollective, domain::DomainSpecs)
    collective.bottom .= mod.(collective.bottom, SVector(domain.domain_width, domain.domain_height))
    collective.top .= mod.(collective.top, SVector(domain.domain_width, domain.domain_height))
end

function apply_hard_wall_boundary_conditions!(collective::CellCollective, domain::DomainSpecs, p::SimulationParameters)
    # Reflect cells off the walls
    for i in eachindex(collective.bottom)
        if collective.bottom[i][1] < p.cell_hard_radius
            collective.bottom[i] = setindex(collective.bottom[i], p.cell_hard_radius, 1)
            collective.theta[i] = pi - collective.theta[i] # Reflect angle
        elseif collective.bottom[i][1] > domain.domain_width - p.cell_hard_radius
            collective.bottom[i] = setindex(collective.bottom[i], 2*(domain.domain_width - p.cell_hard_radius) - collective.bottom[i][1], 1)
            collective.theta[i] = pi - collective.theta[i] # Reflect angle
        end
        if collective.bottom[i][2] < p.cell_hard_radius
            collective.bottom[i] = setindex(collective.bottom[i], p.cell_hard_radius, 2)
            collective.theta[i] = -collective.theta[i] # Reflect angle
        elseif collective.bottom[i][2] > domain.domain_height - p.cell_hard_radius
            collective.bottom[i] = setindex(collective.bottom[i], 2*(domain.domain_height - p.cell_hard_radius) - collective.bottom[i][2], 2)
            collective.theta[i] = -collective.theta[i] # Reflect angle
        end
        if collective.top[i][1] < p.cell_hard_radius
            collective.top[i] = setindex(collective.top[i], p.cell_hard_radius, 1)
        elseif collective.top[i][1] > domain.domain_width - p.cell_hard_radius
            collective.top[i] = setindex(collective.top[i], 2*(domain.domain_width - p.cell_hard_radius) - collective.top[i][1], 1)
        end
        if collective.top[i][2] < p.cell_hard_radius
            collective.top[i] = setindex(collective.top[i], p.cell_hard_radius, 2)
        elseif collective.top[i][2] > domain.domain_height - p.cell_hard_radius
            collective.top[i] = setindex(collective.top[i], 2*(domain.domain_height - p.cell_hard_radius) - collective.top[i][2], 2)
        end
    end
end

function compute_cell_cohesion_forces!(collective::CellCollective, p::SimulationParameters)
    for i in 1:p.num_cells
        vector_to_top = collective.top[i] - collective.bottom[i]
        
        internal_distance = norm(vector_to_top)
        unit_vector_to_top = vector_to_top / internal_distance
        if internal_distance < 2 * p.cell_hard_radius
            collective.top[i] += (unit_vector_to_top * (2 * p.cell_hard_radius - internal_distance))/2
            collective.bottom[i] -= (unit_vector_to_top * (2 * p.cell_hard_radius - internal_distance))/2
        elseif internal_distance < p.maximal_stretch
            collective.forcesbottom[i] += p.internal_cohesion_strength * (internal_distance - 2 * p.cell_hard_radius) * unit_vector_to_top
            collective.forcestop[i] -= p.internal_cohesion_strength * (internal_distance - 2 * p.cell_hard_radius) * unit_vector_to_top
        else
            collective.top[i] -= (unit_vector_to_top * (internal_distance - p.maximal_stretch)) 
            #collective.bottom[i] += (unit_vector_to_top * (internal_distance - p.maximal_stretch)) / 2
        end
    end
end

function _toml_value(data, key::AbstractString)
    if data isa AbstractDict
        if haskey(data, key)
            return data[key]
        end
        for value in values(data)
            found_value = _toml_value(value, key)
            if found_value !== nothing
                return found_value
            end
        end
    end
    return nothing
end

function read_parameters(filename::String)
    parameters = TOML.parsefile(filename)

    simulation_parameters = SimulationParameters(
        Int(_toml_value(parameters, "num_cells")),
        Float64(_toml_value(parameters, "dt")),
        Float64(_toml_value(parameters, "total_time")),
        Float64(_toml_value(parameters, "D_angle")),
        Float64(_toml_value(parameters, "D_position")),
        Float64(_toml_value(parameters, "theta_alignment")),
        Float64(_toml_value(parameters, "run_to_tumble_rate")),
        Float64(_toml_value(parameters, "tumble_to_run_rate")),
        Float64(_toml_value(parameters, "cell_hard_radius")),
        Float64(_toml_value(parameters, "cell_soft_radius")),
        Float64(_toml_value(parameters, "internal_cohesion_strength")),
        Float64(_toml_value(parameters, "maximal_stretch")),
        Float64(_toml_value(parameters, "run_speed"))
    )

    domain_specs = DomainSpecs(
        Float64(_toml_value(parameters, "domain_width")),
        Float64(_toml_value(parameters, "domain_height")),
    )

    return simulation_parameters, domain_specs
end

