# Headless multi-scenario metrics runner (no rendering).
#
# Runs ONE physics block (a base TOML) against three crowding regimes — few(3),
# small(30), tissue(400) — by overriding only count/domain/init/margins/total_time,
# then prints diagnostics used to judge adhesion / collapse / vortex / drift.
#
#   julia --project=. scripts/scan_metrics.jl [base_toml]
#
# Default base_toml = data/sims/parameters_adhesion.toml
#
# NOTE: compute_state_changes! reads a module-global `domain`, so we set
# `global domain` before every simulation_loop call.

using DrWatson, Printf
@quickactivate "Cell_SE_Model"
include(srcdir("simulation_logic.jl"))   # pulls in cell_model.jl (TOML, StaticArrays, ...)

const BASE_TOML = length(ARGS) >= 1 ? ARGS[1] : datadir("sims", "parameters_adhesion.toml")

# ----------------------------------------------------------------------------- helpers
finite_positions(sol) = all(p -> isfinite(p[1]) && isfinite(p[2]),
                            vcat(sol.bottom[end], sol.top[end]))

com(points) = SVector(mean(getindex.(points, 1)), mean(getindex.(points, 2)))

function radius_of_gyration(points)
    c = com(points)
    sqrt(mean(sum.(abs2, (points .- Ref(c)))))
end

# mean nearest-neighbour distance over a set of 2D points (O(n^2))
function mean_nn(points)
    n = length(points)
    n < 2 && return NaN
    s = 0.0
    for i in 1:n
        best = Inf
        for j in 1:n
            i == j && continue
            d = norm(points[i] - points[j])
            d < best && (best = d)
        end
        s += best
    end
    s / n
end

# connected components with an edge whenever bottom distance < radius (union-find)
function n_components(points, radius)
    n = length(points)
    n == 0 && return 0
    parent = collect(1:n)
    find(x) = (while parent[x] != x; parent[x] = parent[parent[x]]; x = parent[x]; end; x)
    for i in 1:n, j in i+1:n
        if norm(points[i] - points[j]) < radius
            ri, rj = find(i), find(j)
            ri != rj && (parent[ri] = rj)
        end
    end
    length(unique(find(i) for i in 1:n))
end

# polar order |<(cosθ, sinθ)>|
function polar_order(thetas)
    isempty(thetas) && return NaN
    norm(SVector(mean(cos.(thetas)), mean(sin.(thetas))))
end

# signed rotational order about COM: <(r x v) / (|r||v|)>; |.|~1 => coherent vortex
function rotational_order(prev, curr)
    n = length(curr)
    n < 2 && return NaN
    c = com(curr)
    s = 0.0; cnt = 0
    for i in 1:n
        r = curr[i] - c
        v = curr[i] - prev[i]
        nr = norm(r); nv = norm(v)
        (nr < 1e-9 || nv < 1e-9) && continue
        s += (r[1]*v[2] - r[2]*v[1]) / (nr * nv)
        cnt += 1
    end
    cnt == 0 ? NaN : s / cnt
end

# fraction near top/bottom walls vs central band (|y-H/2| < 0.15H)
function wall_centre_frac(points, H)
    n = length(points)
    n == 0 && return (NaN, NaN)
    band = 0.15 * H
    wall = count(p -> min(p[2], H - p[2]) < band, points) / n
    centre = count(p -> abs(p[2] - H/2) < band, points) / n
    (wall, centre)
end

# preferred Morse spacing for reference (force-zero), shifted back to true distance
function preferred_spacing(m::MorsePotential, r_well)
    (m.xi1 == m.xi2 || m.v0 <= 0) && return NaN
    d_eff = log(m.u0 / m.v0) * m.xi1 * m.xi2 / (m.xi2 - m.xi1)
    r_well + d_eff
end

# ----------------------------------------------------------------------------- runner
function run_scenario(name, base::AbstractDict, overrides)
    d = deepcopy(base)
    for (k, v) in overrides
        d[k] = v
    end
    tmp = datadir("sims", "_tmp_scan.toml")
    open(tmp, "w") do io; TOML.print(io, d); end
    p, dom, mtop, mbot = read_parameters(tmp)
    global domain = dom                      # compute_state_changes! reads this global
    sol, p, dom = simulation_loop(p, dom, mtop, mbot)

    exploded = !finite_positions(sol)
    botf, topf = sol.bottom[end], sol.top[end]
    bot0 = sol.bottom[1]
    adh_range = 2 * p.adhesion_radius          # adhesion spring acts below this top distance
    rg0, rg = radius_of_gyration(bot0), radius_of_gyration(botf)
    drift = com(botf)[1] - com(bot0)[1]
    rot = rotational_order(sol.bottom[end-1], botf)
    wallf, centref = wall_centre_frac(botf, dom.domain_height)
    ncomp = n_components(botf, 2 * p.cell_soft_radius)

    @printf("%-7s N=%-4d | explode=%-5s top_nn=%6.2f bot_nn=%6.2f (adh<%4.1f) | Rg %6.2f->%6.2f | drift_x=%7.2f | polar=%.2f rot=%+.2f | wall=%.2f centre=%.2f | comps=%d\n",
            name, p.num_cells, string(exploded), mean_nn(topf), mean_nn(botf),
            adh_range, rg0, rg, drift, polar_order(sol.theta[end]), rot, wallf, centref, ncomp)
    return (; name, exploded, top_nn = mean_nn(topf), bot_nn = mean_nn(botf),
            adh_range, rg0, rg, drift, rot, wallf, centref, ncomp)
end

base = TOML.parsefile(BASE_TOML)

println("base = ", BASE_TOML)
println("="^140)

# All three are left-seeded channels so we can read cohesion AND cluster travel.
run_scenario("small", base, Dict(
    "num_cells" => 10, "domain_width" => 800, "domain_height" => 100,
    "init_width" => 40, "init_height" => 60, "margin_width" => 5, "margin_height" => 20,
    "total_time" => 150))

run_scenario("medium", base, Dict(
    "num_cells" => 50, "domain_width" => 1000, "domain_height" => 120,
    "init_width" => 90, "init_height" => 100, "margin_width" => 5, "margin_height" => 10,
    "total_time" => 150))

run_scenario("big", base, Dict(
    "num_cells" => 300, "domain_width" => 1600, "domain_height" => 50,
    "init_width" => 600, "init_height" => 50, "margin_width" => 0, "margin_height" => 0,
    "total_time" => 200))

println("="^140)
