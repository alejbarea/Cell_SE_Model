using GLMakie, GeometryBasics, Colors

function build_segments(bottom, theta, soft_radius)
    n = min(length(bottom), length(theta))
    seg = Vector{GeometryBasics.Point2f}(undef, 2 * n)
    for i in 1:n
        base = Point2f(bottom[i]...)
        dx = soft_radius * cos(theta[i])
        dy = soft_radius * sin(theta[i])
        seg[2i - 1] = base
        seg[2i] = base + Point2f(dx, dy)
    end
    return seg
end

function build_soft_colors(state)
    return [s == 1 ? RGBAf(0, 1, 0, 0.3) : RGBAf(1, 0, 0, 0.3) for s in state]
end

function animate_solution(solution, p, domain, output_path)
    fig = Figure(size = (1200,900))

    ax = Axis(fig[1, :],
        limits = (0, domain.domain_width, 0, domain.domain_height), aspect = DataAspect()
    )

    #=
    ax2 = Axis(fig[2, 1],
        limits = (0, domain.domain_width, 0, domain.domain_height)
    )
    ax3 = Axis(fig[2, 2],
        limits = (0, domain.domain_width, 0, domain.domain_height)
    )
        =#

    bottom_obs = Observable(solution.bottom[1])
    top_obs    = Observable(solution.top[1])
    state_obs  = Observable(solution.state[1])
    theta_obs  = Observable(solution.theta[1])

    soft_colors = build_soft_colors(state_obs[])
    soft_colors_obs = Observable(soft_colors)
    scatter!(ax, bottom_obs, color = soft_colors_obs, markersize = 2*p.cell_soft_radius, marker = Makie.Circle,markerspace=:data)
    scatter!(ax, top_obs, color = soft_colors_obs, markersize = 2*p.cell_soft_radius, marker = Makie.Circle,markerspace=:data)

    scatter!(ax, bottom_obs, color = :blue, markersize = 2*p.cell_hard_radius, marker = Makie.Circle,markerspace=:data)
    scatter!(ax, top_obs, color = :red, markersize = 2*p.cell_hard_radius, marker = Makie.Circle,markerspace=:data)
    
    # Add arrows showing direction (theta angle) at each cell's bottom position
    seg0 = build_segments(bottom_obs[], theta_obs[], p.cell_soft_radius)
    segs = Observable(seg0)
    linesegments!(ax, segs, color = :black, linewidth = 2)
    #arrows2d!(ax, bottom_obs, arrow_directions, shaftwidth = 2, shaftlength = 2, color = :black)
    #=
    scatter!(ax2, bottom_obs, color = soft_colors_obs, markersize = 2*p.cell_soft_radius, marker = Makie.Circle,markerspace=:data)
    scatter!(ax3, top_obs, color = soft_colors_obs, markersize = 2*p.cell_soft_radius, marker = Makie.Circle,markerspace=:data)

    scatter!(ax2, bottom_obs, color = :blue, markersize = 2*p.cell_hard_radius, marker = Makie.Circle,markerspace=:data)
    scatter!(ax3, top_obs, color = :red, markersize = 2*p.cell_hard_radius, marker = Makie.Circle,markerspace=:data)
    =#




    record(fig, output_path; framerate = 30) do io
        for step in eachindex(solution.bottom)
            bottom = solution.bottom[step]
            top = solution.top[step]
            state = solution.state[step]
            theta = solution.theta[step]

            # update data
            bottom_obs[] = bottom
            top_obs[] = top
            state_obs[] = state
            theta_obs[] = theta

            segs[] = build_segments(bottom, theta, p.cell_soft_radius)
            soft_colors_obs[] = build_soft_colors(state)
            recordframe!(io)
        end
    end

    return output_path
end