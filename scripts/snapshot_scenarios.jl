# Save final-frame PNG snapshots of the small/medium/big scenarios so packing /
# cohesion can be inspected without watching the full MP4.
#   julia --project=. scripts/snapshot_scenarios.jl [base_toml]

using DrWatson
@quickactivate "Cell_SE_Model"
include(srcdir("simulation_logic.jl"))
using GLMakie, Colors

const BASE_TOML = length(ARGS) >= 1 ? ARGS[1] : datadir("sims", "parameters_adhesion.toml")

function snapshot(name, base, overrides)
    d = deepcopy(base)
    for (k, v) in overrides; d[k] = v; end
    tmp = datadir("sims", "_tmp_snap.toml")
    open(tmp, "w") do io; TOML.print(io, d); end
    p, dom, mtop, mbot = read_parameters(tmp)
    global domain = dom
    sol, p, dom = simulation_loop(p, dom, mtop, mbot)

    bot = sol.bottom[end]; top = sol.top[end]; st = sol.state[end]
    fig = Figure(size = (900, 300))
    ax = Axis(fig[1, 1], title = "$name  N=$(p.num_cells)  (final frame)",
              limits = (0, dom.domain_width, 0, dom.domain_height), aspect = DataAspect())
    # hard-radius discs for bottom (leading) and top (trailing)
    scatter!(ax, [Point2f(b...) for b in top],  markersize = 2*p.cell_hard_radius,
             markerspace = :data, color = (:gray, 0.5))
    scatter!(ax, [Point2f(b...) for b in bot],  markersize = 2*p.cell_hard_radius,
             markerspace = :data, color = [s == 1 ? RGBAf(0,0.6,0,0.7) : RGBAf(0.8,0,0,0.7) for s in st])
    out = plotsdir("snap_$(name).png"); mkpath(dirname(out))
    save(out, fig); println("saved $out")
end

base = TOML.parsefile(BASE_TOML)
snapshot("small", base, Dict("num_cells"=>10,"domain_width"=>800,"domain_height"=>100,
    "init_width"=>40,"init_height"=>60,"margin_width"=>5,"margin_height"=>20,"total_time"=>300))
snapshot("medium", base, Dict("num_cells"=>50,"domain_width"=>1000,"domain_height"=>120,
    "init_width"=>90,"init_height"=>100,"margin_width"=>5,"margin_height"=>10,"total_time"=>300))
snapshot("big", base, Dict("num_cells"=>300,"domain_width"=>1600,"domain_height"=>50,
    "init_width"=>600,"init_height"=>50,"margin_width"=>0,"margin_height"=>0,"total_time"=>500))
