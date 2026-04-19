using GLMakie

function animate_solution(solution, domain, output_path)
    fig = Figure(size = (800, 800))

    ax = Axis(fig[1, 1],
        limits = (0, domain.domain_width, 0, domain.domain_height)
    )

    bottom_obs = Observable(solution.bottom[1])
    top_obs    = Observable(solution.top[1])

    scatter!(ax, bottom_obs, color = :blue)
    scatter!(ax, top_obs, color = :red)

    record(fig, output_path; framerate = 30) do io
        for step in eachindex(solution.bottom)

            # update data
            bottom_obs[] = solution.bottom[step]
            top_obs[]    = solution.top[step]

            recordframe!(io)
        end
    end

    return output_path
end