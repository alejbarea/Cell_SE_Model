using DrWatson
@quickactivate "Cell_SE_Model"

include(srcdir("simulation_logic.jl"))
include(scriptsdir("animation_logic.jl"))

solution = simulation_loop()

animation_path = datadir("sims", "animation.mp4")
animate_solution(solution, DomainSpecs(1000.0, 1000.0), animation_path)