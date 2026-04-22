using GLMakie

function animate_solution(solution, p, domain, output_path)
    fig = Figure(size = (800,800))

    ax = Axis(fig[1, :],
        limits = (0, domain.domain_width, 0, domain.domain_height)
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

    color_mask = state_obs[] .== 1
    soft_colors = RGBAf.(0, 1, 0, 0.3) .* color_mask + RGBAf.(1, 0, 0, 0.3) .* .!color_mask
    soft_colors_obs = Observable(soft_colors)
    scatter!(ax, bottom_obs, color = soft_colors_obs, markersize = 2*p.cell_soft_radius, marker = Makie.Circle,markerspace=:data)
    scatter!(ax, top_obs, color = soft_colors_obs, markersize = 2*p.cell_soft_radius, marker = Makie.Circle,markerspace=:data)

    scatter!(ax, bottom_obs, color = :blue, markersize = 2*p.cell_hard_radius, marker = Makie.Circle,markerspace=:data)
    scatter!(ax, top_obs, color = :red, markersize = 2*p.cell_hard_radius, marker = Makie.Circle,markerspace=:data)
    #=
    scatter!(ax2, bottom_obs, color = soft_colors_obs, markersize = 2*p.cell_soft_radius, marker = Makie.Circle,markerspace=:data)
    scatter!(ax3, top_obs, color = soft_colors_obs, markersize = 2*p.cell_soft_radius, marker = Makie.Circle,markerspace=:data)

    scatter!(ax2, bottom_obs, color = :blue, markersize = 2*p.cell_hard_radius, marker = Makie.Circle,markerspace=:data)
    scatter!(ax3, top_obs, color = :red, markersize = 2*p.cell_hard_radius, marker = Makie.Circle,markerspace=:data)
    =#




    record(fig, output_path; framerate = 30) do io
        for step in eachindex(solution.bottom)

            # update data
            bottom_obs[] = solution.bottom[step]
            top_obs[]    = solution.top[step]
            state_obs[]  = solution.state[step]

             # update colors
            color_mask = state_obs[] .== 1
            soft_colors = RGBAf.(0, 1, 0, 0.3) .* color_mask + RGBAf.(1, 0, 0, 0.3) .* .!color_mask
            soft_colors_obs[] = soft_colors
            recordframe!(io)
        end
    end

    return output_path
end