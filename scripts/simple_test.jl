using DrWatson
@quickactivate "Cell_SE_Model"

include(srcdir("simulation_logic.jl"))
include(scriptsdir("animation_logic.jl"))
include(srcdir("cell_model.jl"))

p, domain, morse_potential_top, morse_potential_bottom = read_parameters(datadir("sims", "parameters.toml"))



solution, p, domain = simulation_loop(p, domain, morse_potential_top, morse_potential_bottom)

animation_path = datadir("sims", "animation.mp4")
animate_solution(solution, p, domain, animation_path)