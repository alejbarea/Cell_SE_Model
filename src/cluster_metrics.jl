# Cluster-cohesion and rightward-penetration metrics.
#
# These operate on the `Solution` object returned by `simulation_loop`
# (see src/simulation_logic.jl). The bottom surface is the authoritative
# cell position, so every metric is computed from `solution.bottom`.
#
# To damp run-to-run stochastic noise, each metric is averaged over the
# final `window_fraction` of timesteps.

using Statistics

"""
    window_range(num_steps, window_fraction)

Step indices forming the trailing `window_fraction` of the trajectory
(always at least the final step, never earlier than step 2).
"""
function window_range(num_steps::Int, window_fraction::Real)
    w = floor(Int, window_fraction * num_steps)
    return max(2, num_steps - w + 1):num_steps
end

"""
    compute_cluster_cohesion(solution, p; window_fraction = 0.1)

Radius of gyration of the bottom positions, averaged over the trailing
window. LOWER means a tighter, more cohesive cluster.
"""
function compute_cluster_cohesion(solution, p; window_fraction = 0.1)
    num_steps = length(solution.bottom)
    rgs = Float64[]
    for step in window_range(num_steps, window_fraction)
        pts = solution.bottom[step]
        centroid = mean(pts)
        push!(rgs, sqrt(mean(sum(abs2, pt - centroid) for pt in pts)))
    end
    return mean(rgs)
end

"""
    compute_rightward_penetration(solution, p; window_fraction = 0.1)

Advance of the cluster centroid in +x, from step 1 to the trailing
window. HIGHER means the cells penetrated further right in the domain.
"""
function compute_rightward_penetration(solution, p; window_fraction = 0.1)
    num_steps = length(solution.bottom)
    x_start = mean(pt[1] for pt in solution.bottom[1])
    x_end = mean(mean(pt[1] for pt in solution.bottom[step])
                 for step in window_range(num_steps, window_fraction))
    return x_end - x_start
end

"""
    compute_front_position(solution, p; window_fraction = 0.1, front_fraction = 0.1)

Mean x of the rightmost `front_fraction` of cells (the leading edge of
the group), averaged over the trailing window. Reported as an extra
diagnostic alongside the two optimisation objectives.
"""
function compute_front_position(solution, p; window_fraction = 0.1, front_fraction = 0.1)
    num_steps = length(solution.bottom)
    fronts = Float64[]
    for step in window_range(num_steps, window_fraction)
        xs = sort([pt[1] for pt in solution.bottom[step]])
        k = max(1, ceil(Int, front_fraction * length(xs)))
        push!(fronts, mean(@view xs[end-k+1:end]))
    end
    return mean(fronts)
end

"""
    pareto_mask(cohesion, penetration)

Boolean mask of the Pareto-optimal entries when MINIMISING `cohesion`
(radius of gyration) and MAXIMISING `penetration`. Non-finite entries
(e.g. blown-up simulations) are never marked optimal.
"""
function pareto_mask(cohesion::AbstractVector, penetration::AbstractVector)
    n = length(cohesion)
    mask = falses(n)
    for i in 1:n
        (isfinite(cohesion[i]) && isfinite(penetration[i])) || continue
        dominated = false
        for j in 1:n
            i == j && continue
            (isfinite(cohesion[j]) && isfinite(penetration[j])) || continue
            if cohesion[j] <= cohesion[i] && penetration[j] >= penetration[i] &&
               (cohesion[j] < cohesion[i] || penetration[j] > penetration[i])
                dominated = true
                break
            end
        end
        mask[i] = !dominated
    end
    return mask
end
