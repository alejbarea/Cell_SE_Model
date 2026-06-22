# Render the three crowding regimes (small / medium / big) from ONE physics
# block (data/sims/parameters.toml), overriding only count/domain/init/total_time.
# Outputs data/sims/animation_{small,medium,big}.mp4.
#
#   julia --project=. scripts/render_scenarios.jl [base_toml] [which]
#
# `which` (optional) = small | medium | big to render just one.
#
# NOTE: compute_state_changes! reads a module-global `domain`; we set it per run.

using DrWatson
@quickactivate "Cell_SE_Model"
include(srcdir("simulation_logic.jl"))
include(scriptsdir("animation_logic.jl"))

const BASE_TOML = length(ARGS) >= 1 ? ARGS[1] : datadir("sims", "parameters.toml")
const WHICH     = length(ARGS) >= 2 ? ARGS[2] : "all"

function render_scenario(name, base::AbstractDict, overrides)
    d = deepcopy(base)
    for (k, v) in overrides
        d[k] = v
    end
    tmp = datadir("sims", "_tmp_render.toml")
    open(tmp, "w") do io; TOML.print(io, d); end
    p, dom, mtop, mbot = read_parameters(tmp)
    global domain = dom
    sol, p, dom = simulation_loop(p, dom, mtop, mbot)
    out = datadir("sims", "animation_$(name).mp4")
    animate_solution(sol, p, dom, out)
    println("rendered $out  (N=$(p.num_cells))")
end

base = TOML.parsefile(BASE_TOML)

# total_time generous so cluster cohesion + travel are visible over the run.
scenarios = Dict(
    "small" => Dict("num_cells" => 10, "domain_width" => 800, "domain_height" => 100,
        "init_width" => 40, "init_height" => 60, "margin_width" => 5, "margin_height" => 20,
        "total_time" => 300),
    "medium" => Dict("num_cells" => 50, "domain_width" => 1000, "domain_height" => 120,
        "init_width" => 90, "init_height" => 100, "margin_width" => 5, "margin_height" => 10,
        "total_time" => 300),
    "big" => Dict("num_cells" => 300, "domain_width" => 1600, "domain_height" => 50,
        "init_width" => 600, "init_height" => 50, "margin_width" => 0, "margin_height" => 0,
        "total_time" => 500),
)

order = ["small", "medium", "big"]
for name in order
    (WHICH == "all" || WHICH == name) && render_scenario(name, base, scenarios[name])
end
