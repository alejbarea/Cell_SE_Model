using StaticArrays, Distributions, TOML, LinearAlgebra, CellListMap

abstract type Parameters end

struct SimulationParameters <: Parameters
    num_cells::Int
    dt::Float64
    total_time::Float64
    init_width::Float64
    init_height::Float64
    D_angle::Float64
    D_position::Float64
    theta_alignment::Float64
    run_to_tumble_rate_min::Float64
    run_to_tumble_rate_max::Float64
    run_to_tumble_rate_function::String
    run_to_tumble_rate_allure::Float64
    tumble_to_run_rate_min::Float64
    tumble_to_run_rate_max::Float64
    tumble_to_run_rate_function::String
    tumble_to_run_rate_allure::Float64
    cell_hard_radius::Float64
    cell_soft_radius::Float64
    maximal_stretch::Float64
    run_speed_min::Float64
    run_speed_max::Float64
    run_speed_dist::String
    top_internal_pull_strength::Float64
    soft_boundary_potential_strength::Float64
    hard_boundary_velocity_behaviour::String
    bottom_internal_pull_strength::Float64
    cil_intensity::Float64
    cil_wall_intensity::Float64
end

struct DomainSpecs <: Parameters
    domain_width::Float64
    domain_height::Float64
end

struct MorsePotential <: Parameters
    u0::Float64
    v0::Float64
    xi1::Float64
    xi2::Float64
    cutoff::Float64
end

# Cell Collective Structure is a collection of all cell data for easy vectorization.
struct CellCollective
    bottom::Vector{SVector{2, Float64}} # Coordinates of cell centers at z = 0
    top::Vector{SVector{2, Float64}} # Coordinates of cell centers at z = z0
    theta::Vector{Float64} # Cell angles
    thetabar::Vector{SVector{2, Float64}} # Cell desired polarity vectors
    neighboursbottom::Vector{Vector{Int}} # Neighbouring cell indices
    neighbourstop::Vector{Vector{Int}} # Neighbouring cell indices
    forcesbottom::Vector{SVector{2, Float64}} # Forces on each cell at z = 0
    forcestop::Vector{SVector{2, Float64}} # Forces on each cell at z = z0
    state::Vector{Int} # State of each cell (e.g., 0 for tumbling, 1 for running)
    state_timer::Vector{Float64} # Timer for state transitions
    run_speeds::Vector{Float64} # Speed of each cell when in the running state
    hitting_wall::Vector{Bool} # Whether each cell is currently hitting a wall
end

function sample_speed_distribution(p::SimulationParameters)
    if p.run_speed_dist == "constant"
        return fill(p.run_speed_min, p.num_cells)
    elseif p.run_speed_dist == "uniform"
        return rand(Uniform(p.run_speed_min, p.run_speed_max), p.num_cells)
    elseif p.run_speed_dist == "normal"
        return rand(Normal((p.run_speed_min + p.run_speed_max) / 2, (p.run_speed_max - p.run_speed_min) / 6), p.num_cells)
    else
        return fill(p.run_speed_min, p.num_cells)
    end
end

function initialize_cell_collective(p::SimulationParameters, domain::DomainSpecs)
    bottom = Vector{SVector{2, Float64}}(undef, p.num_cells)
    top = Vector{SVector{2, Float64}}(undef, p.num_cells)
    theta = Vector{Float64}(undef, p.num_cells)
    thetabar = Vector{SVector{2, Float64}}(undef, p.num_cells)
    neighboursbottom = Vector{Vector{Int}}(undef, p.num_cells)
    neighbourstop = Vector{Vector{Int}}(undef, p.num_cells)
    forcesbottom = Vector{SVector{2, Float64}}(undef, p.num_cells)
    forcestop = Vector{SVector{2, Float64}}(undef, p.num_cells)
    state = Vector{Int}(undef, p.num_cells)
    state_timer = Vector{Float64}(undef, p.num_cells)
    run_speeds = sample_speed_distribution(p)
    hitting_wall = Vector{Bool}(undef, p.num_cells)
    for i in 1:p.num_cells
        bottom[i] = SVector(rand() * p.init_width, rand() * p.init_height)
        theta[i] = rand() * 2 * pi # Random initial angle
        top[i] = bottom[i] - p.cell_hard_radius * SVector(cos(theta[i]), sin(theta[i])) # Initial top position directly above bottom
        thetabar[i] = SVector(cos(theta[i]), sin(theta[i])) # Initial desired angle same as initial angle
        forcesbottom[i] = SVector(0.0, 0.0)
        forcestop[i] = SVector(0.0, 0.0)
        rand_aux = rand()  
        time_quot = p.tumble_to_run_rate_min / (p.tumble_to_run_rate_min + p.run_to_tumble_rate_max)      
        state[i] = rand_aux < time_quot ? 1 : 0 # Random initial state
        state_timer[i] = rand(Exponential(1.0)) # Random initial timer for state transitions
        neighboursbottom[i] = []
        neighbourstop[i] = []
        hitting_wall[i] = false
    end
    return CellCollective(bottom, top, theta, thetabar, neighboursbottom, neighbourstop, forcesbottom, forcestop, state, state_timer, run_speeds, hitting_wall)
end



function update_cell_collective!(collective::CellCollective, p::SimulationParameters)
    # Update bottom and top coordinates
    collective.bottom .+= p.dt .* collective.forcesbottom
    collective.top .+= p.dt .* collective.forcestop
    neighbour_mask = (length.(collective.neighboursbottom) .== 0) .& (collective.hitting_wall .== false)
    collective.thetabar[neighbour_mask] .= SVector.(cos.(collective.theta[neighbour_mask]), sin.(collective.theta[neighbour_mask]))    # Update theta only where state == 0
    mask = collective.state .== 0
    thetabar_angle_conversion = atan.(getindex.(collective.thetabar, 2), getindex.(collective.thetabar, 1))
    collective.theta[mask] .+= -p.theta_alignment .* sin.(collective.theta[mask] - thetabar_angle_conversion[mask]) .* p.dt + sqrt(2 * p.D_angle * p.dt) .* randn(sum(mask))
    collective.bottom[.!mask] .+= collective.run_speeds[.!mask] * p.dt .* SVector.(cos.(collective.theta[.!mask]), sin.(collective.theta[.!mask])) 
    fill!(collective.forcesbottom, SVector(0.0, 0.0))
    fill!(collective.forcestop,    SVector(0.0, 0.0))
    empty!.(collective.neighboursbottom)
    empty!.(collective.neighbourstop)
    fill!(collective.thetabar, SVector(0.0, 0.0))
    fill!(collective.hitting_wall, false)
end


function morse_interaction_forces(potential::MorsePotential, unit_vector::SVector{2, Float64}, distance::Float64)
    if distance > eps()
        force_magnitude = (potential.u0 * exp(-distance / potential.xi1) - potential.v0 * exp(-distance / potential.xi2))
        return force_magnitude * unit_vector
    else
        return SVector(0.0, 0.0)
    end
end

function compute_interaction_forces!(collective::CellCollective, p::SimulationParameters, morse_potential_top::MorsePotential, morse_potential_bottom::MorsePotential)
    for i in 1:p.num_cells
        for j in i+1:p.num_cells
            distance_vector_bottom = collective.bottom[i] - collective.bottom[j]
            distance_vector_top = collective.top[i] - collective.top[j]
            distance_bottom = norm(distance_vector_bottom)
            distance_top = norm(distance_vector_top)
            unit_vector_bottom = distance_vector_bottom / distance_bottom
            unit_vector_top = distance_vector_top / distance_top
            if distance_bottom < 2 * p.cell_hard_radius
                push!(collective.neighboursbottom[i], j)
                push!(collective.neighboursbottom[j], i)
                force_magnitude = (2 * p.cell_hard_radius - distance_bottom)
                force_direction = distance_vector_bottom / distance_bottom
                collective.bottom[i] += force_magnitude * force_direction / 2
                collective.bottom[j] -= force_magnitude * force_direction / 2
            elseif distance_bottom < 2 * p.cell_soft_radius
                push!(collective.neighboursbottom[i], j)
                push!(collective.neighboursbottom[j], i)
            end
            if distance_top < 2 * p.cell_hard_radius
                push!(collective.neighbourstop[i], j)
                push!(collective.neighbourstop[j], i)
                force_magnitude = (2 * p.cell_hard_radius - distance_top)
                force_direction = distance_vector_top / distance_top
                collective.top[i] += force_magnitude * force_direction / 2
                collective.top[j] -= force_magnitude * force_direction / 2
                collective.thetabar[i] += unit_vector_top
                collective.thetabar[j] -= unit_vector_top
            elseif distance_top < 2 * p.cell_soft_radius
                push!(collective.neighbourstop[i], j)
                push!(collective.neighbourstop[j], i)
                collective.thetabar[i] += unit_vector_top
                collective.thetabar[j] -= unit_vector_top
            end
            if distance_top < morse_potential_top.cutoff && distance_top > 2 * p.cell_hard_radius
                force_top = morse_interaction_forces(morse_potential_top, unit_vector_top, distance_top - 2*p.cell_hard_radius)
                collective.forcestop[i] += force_top / 2
                collective.forcestop[j] -= force_top / 2
            end
            force_bottom = morse_interaction_forces(morse_potential_bottom, unit_vector_bottom, distance_bottom - 2*p.cell_soft_radius)
            collective.forcesbottom[i] += force_bottom / 2
            collective.forcesbottom[j] -= force_bottom / 2
        end
    end
end


function compute_stochastic_forces!(collective::CellCollective, p::SimulationParameters)
        for i in 1:p.num_cells
            xi = sqrt(2 * p.D_position * p.dt) .* randn(2)
            collective.forcesbottom[i] += xi
            collective.forcestop[i] += xi
        end
        
end

function compute_tumble_to_run_rate(x::Float64, p::SimulationParameters)
    rate_fct = p.tumble_to_run_rate_function
    allure = p.tumble_to_run_rate_allure
    if rate_fct == "constant"
        return p.tumble_to_run_rate_min
    elseif rate_fct == "sigmoid"
        return p.tumble_to_run_rate_min + (p.tumble_to_run_rate_max - p.tumble_to_run_rate_min) / (1 + exp(-allure * x))
    else
        return p.tumble_to_run_rate_min
    end
end

function compute_run_to_tumble_rate(x::Float64, p::SimulationParameters)
    rate_fct = p.run_to_tumble_rate_function
    allure = p.run_to_tumble_rate_allure
    if rate_fct == "constant"
        return p.run_to_tumble_rate_min
    elseif rate_fct == "sigmoid"
        return p.run_to_tumble_rate_max + (p.run_to_tumble_rate_min - p.run_to_tumble_rate_max) / (1 + exp(-allure * x))
    else
        return p.run_to_tumble_rate_min
    end
end

function compute_state_changes!(collective::CellCollective, p::SimulationParameters)
    

    #collective.state_timer .-= p.dt*((1 .- collective.state) .* tumble_to_run_rate .+ collective.state .* run_to_tumble_rate)
    for i in 1:p.num_cells
        run_to_tumble_rate = compute_run_to_tumble_rate(collective.bottom[i][1], p)
        tumble_to_run_rate = compute_tumble_to_run_rate(collective.bottom[i][1], p)
        if collective.state[i] == 1
            for j in collective.neighboursbottom[i]
                distance_vector_bottom = collective.bottom[j] - collective.bottom[i]
                distance_bottom = norm(distance_vector_bottom)
                unit_vector_bottom = distance_vector_bottom / distance_bottom
                cos35 = cos(35 * pi / 180)
                speed_direction = SVector(cos(collective.theta[i]), sin(collective.theta[i]))
                dot_product = dot(unit_vector_bottom, speed_direction)
                run_to_tumble_rate += p.cil_intensity * (atan((dot_product - cos35) / 0.00001) + pi/2)
            end

            # Check distance to walls
            dist_left = collective.bottom[i][1]
            dist_right = domain.domain_width - collective.bottom[i][1]
            dist_bottom = collective.bottom[i][2]
            dist_top = domain.domain_height - collective.bottom[i][2]
            cos85_wall = cos(85 * pi / 180)
            
            if dist_left < p.cell_soft_radius
                dot_product_wall = -cos(collective.theta[i]) # Dot product with left wall normal (-1,0)
                run_to_tumble_rate += p.cil_wall_intensity * (atan((dot_product_wall - cos85_wall) / 0.00001) + pi/2)
            end
            if dist_right < p.cell_soft_radius
                dot_product_wall = cos(collective.theta[i]) # Dot product with right wall normal (1,0)
                run_to_tumble_rate += p.cil_wall_intensity * (atan((dot_product_wall - cos85_wall) / 0.00001) + pi/2)
            end
            if dist_bottom < p.cell_soft_radius
                dot_product_wall = -sin(collective.theta[i]) # Dot product with bottom wall normal (0,-1)
                run_to_tumble_rate += p.cil_wall_intensity * (atan((dot_product_wall - cos85_wall) / 0.00001) + pi/2)
            end
            if dist_top < p.cell_soft_radius
                dot_product_wall = sin(collective.theta[i]) # Dot product with top wall normal (0,1)
                run_to_tumble_rate += p.cil_wall_intensity * (atan((dot_product_wall - cos85_wall) / 0.00001) + pi/2)
            end

            collective.state_timer[i] -= p.dt * run_to_tumble_rate
                
        else
            collective.state_timer[i] -= p.dt * tumble_to_run_rate
        end
    end
    
    # Check for state transitions
    mask = collective.state_timer .<= 0
    collective.state[mask] .= 1 .- collective.state[mask] # Toggle state
    # Reset timers for cells that changed state
    exponential_dist = Exponential(1.0)
    collective.state_timer[mask] .= rand(exponential_dist,sum(mask))
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
            if p.hard_boundary_velocity_behaviour == "bounce"
                collective.theta[i] = pi - collective.theta[i] # Reflect angle
            end
            collective.hitting_wall[i] = true
            collective.thetabar[i] += SVector(1,0) # Align to the right
        elseif collective.bottom[i][1] > domain.domain_width - p.cell_hard_radius
            collective.bottom[i] = setindex(collective.bottom[i], 2*(domain.domain_width - p.cell_hard_radius) - collective.bottom[i][1], 1)
            if p.hard_boundary_velocity_behaviour == "bounce"
                collective.theta[i] = pi - collective.theta[i] # Reflect angle
            end
            collective.hitting_wall[i] = true
            collective.thetabar[i] += SVector(-1,0) # Align to the left
        end
        if collective.bottom[i][2] < p.cell_hard_radius
            collective.bottom[i] = setindex(collective.bottom[i], p.cell_hard_radius, 2)
            if p.hard_boundary_velocity_behaviour == "bounce"
                collective.theta[i] = -collective.theta[i] # Reflect angle
            end
            collective.hitting_wall[i] = true
            collective.thetabar[i] += SVector(0,1) # Align to the top
        elseif collective.bottom[i][2] > domain.domain_height - p.cell_hard_radius
            collective.bottom[i] = setindex(collective.bottom[i], 2*(domain.domain_height - p.cell_hard_radius) - collective.bottom[i][2], 2)
            if p.hard_boundary_velocity_behaviour == "bounce"
                collective.theta[i] = -collective.theta[i] # Reflect angle
            end
            collective.hitting_wall[i] = true
            collective.thetabar[i] += SVector(0,-1) # Align to the bottom
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

function compute_soft_wall_forces!(collective::CellCollective, domain::DomainSpecs, p::SimulationParameters)
    for i in eachindex(collective.bottom)
        if collective.bottom[i][1] < p.cell_soft_radius
            force_magnitude = p.soft_boundary_potential_strength*(p.cell_soft_radius - collective.bottom[i][1])
            collective.forcesbottom[i] += SVector(force_magnitude, 0.0)
            collective.hitting_wall[i] = true
            collective.thetabar[i] += SVector(1,0) # Align to the right
            
        elseif collective.bottom[i][1] > domain.domain_width - p.cell_soft_radius
            force_magnitude = p.soft_boundary_potential_strength*(collective.bottom[i][1] - (domain.domain_width - p.cell_soft_radius))
            collective.forcesbottom[i] -= SVector(force_magnitude, 0.0)
            collective.hitting_wall[i] = true
            collective.thetabar[i] += SVector(-1,0)
        end
        if collective.bottom[i][2] < p.cell_soft_radius
            force_magnitude = p.soft_boundary_potential_strength*(p.cell_soft_radius - collective.bottom[i][2])
            collective.forcesbottom[i] += SVector(0.0, force_magnitude)
            collective.hitting_wall[i] = true
            collective.thetabar[i] += SVector(0,1)
        elseif collective.bottom[i][2] > domain.domain_height - p.cell_soft_radius
            force_magnitude = p.soft_boundary_potential_strength*(collective.bottom[i][2] - (domain.domain_height - p.cell_soft_radius))
            collective.forcesbottom[i] -= SVector(0.0, force_magnitude)
            collective.hitting_wall[i] = true
            collective.thetabar[i] += SVector(0,-1)
        end
        if collective.top[i][1] < p.cell_soft_radius
            force_magnitude = p.soft_boundary_potential_strength*(p.cell_soft_radius - collective.top[i][1])
            collective.forcestop[i] += SVector(force_magnitude, 0.0)
        elseif collective.top[i][1] > domain.domain_width - p.cell_soft_radius
            force_magnitude = p.soft_boundary_potential_strength*(collective.top[i][1] - (domain.domain_width - p.cell_soft_radius))
            collective.forcestop[i] -= SVector(force_magnitude, 0.0)
        end
        if collective.top[i][2] < p.cell_soft_radius
            force_magnitude = p.soft_boundary_potential_strength*(p.cell_soft_radius - collective.top[i][2])
            collective.forcestop[i] += SVector(0.0, force_magnitude)
        elseif collective.top[i][2] > domain.domain_height - p.cell_soft_radius
            force_magnitude = p.soft_boundary_potential_strength*(collective.top[i][2] - (domain.domain_height - p.cell_soft_radius))
            collective.forcestop[i] -= SVector(0.0, force_magnitude)
        end
    end
end

function compute_cell_cohesion_forces!(collective::CellCollective, p::SimulationParameters)
    for i in 1:p.num_cells
        vector_to_top = collective.top[i] - collective.bottom[i]
        
        internal_distance = norm(vector_to_top)
        if internal_distance > eps()
            unit_vector_to_top = vector_to_top / internal_distance
            if internal_distance < p.maximal_stretch
                #internal_cohesion_force = (p.top_internal_pull_strength/p.bottom_internal_pull_strength * exp(-(internal_distance-2*p.cell_hard_radius)/p.bottom_internal_pull_strength) - p.soft_boundary_potential_strength/p.morse_potential_xi2 * exp(-(internal_distance-p.cell_hard_radius)/p.morse_potential_xi2)) * unit_vector_to_top
                collective.forcestop[i] -= p.top_internal_pull_strength*internal_distance*unit_vector_to_top
                collective.forcesbottom[i] += p.bottom_internal_pull_strength*internal_distance*unit_vector_to_top
            else
                collective.top[i] -= (unit_vector_to_top * (internal_distance - p.maximal_stretch)) / 2
                collective.bottom[i] += (unit_vector_to_top * (internal_distance - p.maximal_stretch)) / 2
            end
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
        Float64(_toml_value(parameters, "init_width")),
        Float64(_toml_value(parameters, "init_height")),
        Float64(_toml_value(parameters, "D_angle")),
        Float64(_toml_value(parameters, "D_position")),
        Float64(_toml_value(parameters, "theta_alignment")),
        Float64(_toml_value(parameters, "run_to_tumble_rate_min")),
        Float64(_toml_value(parameters, "run_to_tumble_rate_max")),
        String(_toml_value(parameters, "run_to_tumble_rate_function")),
        Float64(_toml_value(parameters, "run_to_tumble_rate_allure")),
        Float64(_toml_value(parameters, "tumble_to_run_rate_min")),
        Float64(_toml_value(parameters, "tumble_to_run_rate_max")),
        String(_toml_value(parameters, "tumble_to_run_rate_function")),
        Float64(_toml_value(parameters, "tumble_to_run_rate_allure")),        
        Float64(_toml_value(parameters, "cell_hard_radius")),
        Float64(_toml_value(parameters, "cell_soft_radius")),
        Float64(_toml_value(parameters, "maximal_stretch")),
        Float64(_toml_value(parameters, "run_speed_min")),
        Float64(_toml_value(parameters, "run_speed_max")),
        String(_toml_value(parameters, "run_speed_distribution")),
        Float64(_toml_value(parameters, "top_internal_pull_strength")),
        Float64(_toml_value(parameters, "soft_boundary_potential_strength")),
        String(_toml_value(parameters, "hard_boundary_velocity_behaviour")),
        Float64(_toml_value(parameters, "bottom_internal_pull_strength")),
        Float64(_toml_value(parameters, "cil_intensity")),
        Float64(_toml_value(parameters, "cil_wall_intensity"))
    )

    domain_specs = DomainSpecs(
        Float64(_toml_value(parameters, "domain_width")),
        Float64(_toml_value(parameters, "domain_height")),
    )

    morse_potential_top = MorsePotential(
        Float64(_toml_value(parameters, "top_morse_u0")),
        Float64(_toml_value(parameters, "top_morse_v0")),
        Float64(_toml_value(parameters, "top_morse_xi1")),
        Float64(_toml_value(parameters, "top_morse_xi2")),
        Float64(_toml_value(parameters, "top_morse_cutoff")),
    )

    morse_potential_bottom = MorsePotential(
        Float64(_toml_value(parameters, "bottom_morse_u0")),
        Float64(_toml_value(parameters, "bottom_morse_v0")),
        Float64(_toml_value(parameters, "bottom_morse_xi1")),
        Float64(_toml_value(parameters, "bottom_morse_xi2")),
        Inf
    )

    return simulation_parameters, domain_specs, morse_potential_top, morse_potential_bottom
end

