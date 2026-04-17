include(srcdir("simulation_logic.jl"))

solution = simulation_loop()

function save_solution_csv(solution, output_path::AbstractString)
	mkpath(dirname(output_path))

	open(output_path, "w") do io
		println(io, "time,cell_id,bottom_x,bottom_y,top_x,top_y")
		for step in eachindex(solution.bottom)
			for cell_id in eachindex(solution.bottom[step])
				bottom = solution.bottom[step][cell_id]
				top = solution.top[step][cell_id]
				println(io, join((step, cell_id, bottom[1], bottom[2], top[1], top[2]), ","))
			end
		end
	end

	return output_path
end

save_solution_csv(solution, datadir("sims", "solution.csv"))

    