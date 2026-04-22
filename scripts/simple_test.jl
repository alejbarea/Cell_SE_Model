using DrWatson
@quickactivate "Cell_SE_Model"

include(srcdir("simulation_logic.jl"))
include(scriptsdir("animation_logic.jl"))

solution, p, domain = simulation_loop()

animation_path = datadir("sims", "animation.mp4")
animate_solution(solution, p, domain, animation_path)