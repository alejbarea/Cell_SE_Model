using GLMakie

include(srcdir("simulation_logic.jl"))

function animate_solution(solution::Solution, domain::DomainSpecs, output_path::AbstractString)
    mkpath(dirname(output_path))
    
    fig = Figure(size=(800, 800))
    ax = Axis(fig[1, 1], title="Cell Simulation", xlabel="X", ylabel="Y", limits=(0, domain.domain_width, 0, domain.domain_height))
    scatter_bottom = scatter!(ax, solution.bottom[1], color=:blue, label="Bottom")
    scatter_top = scatter!(ax, solution.top[1], color=:red, label="Top")
    axislegend(ax)

    record(fig, output_path; framerate=30) do io
        for step in eachindex(solution.bottom)
            scatter_bottom[1][] = solution.bottom[step]
            scatter_top[1][] = solution.top[step]
            sleep(0.01) # Adjust sleep for smoother animation
        end
    end

    return output_path
end
