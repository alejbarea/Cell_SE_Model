# Cell_SE_Model — Codebase Reference

## Project Summary

Agent-based simulation of collective epithelial cell migration, written in Julia using the [DrWatson](https://juliadynamics.github.io/DrWatson.jl/stable/) reproducible-science framework. Each cell is modelled as a pseudo-3D object with two coupled 2D endpoints (bottom and top surface), run-tumble locomotion inherited from bacterial chemotaxis, Morse-potential cell-cell interactions, and contact inhibition of locomotion (CIL). The main output is time-series trajectories that can be analysed for collective speed, mean orientation, and group motion; animations are produced via GLMakie.

---

## Repository Structure

```
Cell_SE_Model/
├── src/
│   ├── cell_model.jl          # All structs, physics functions, parameter I/O  (481 lines)
│   ├── simulation_logic.jl    # Simulation loop, Solution struct               (48 lines)
│   └── run_analysis.jl        # Post-hoc analysis functions                   (84 lines)
├── scripts/
│   ├── run_sim_csv_save.jl    # Run one simulation and save trajectories as CSV
│   ├── simple_test.jl         # Run one simulation and render an MP4 animation
│   ├── sweep_logic.jl         # Single-parameter sweep
│   ├── double_sweep_logic.jl  # Two-parameter grid sweep
│   ├── animation_logic.jl     # GLMakie animation renderer
│   └── plotting_logic.jl      # Plotting helpers
├── data/sims/
│   ├── parameters.toml                        # Base simulation config
│   ├── sweep_parameters.toml                  # Single sweep config
│   ├── double_sweep_parameters.toml           # Double sweep config
│   ├── solution.csv                           # Output trajectories
│   ├── analysis_run_speed.csv                 # Sweep analysis output
│   ├── analysis_num_cells.csv
│   ├── analysis_num_cells_cell_hard_radius.csv
│   └── animation.mp4
├── test/runtests.jl           # Parameter-reading tests
├── Project.toml               # Julia package manifest
└── Manifest.toml              # Locked dependency versions
```

DrWatson path helpers used throughout: `srcdir(...)` → `src/`, `datadir(...)` → `data/`.

---

## Data Structures

All defined in [src/cell_model.jl](src/cell_model.jl).

### `SimulationParameters <: Parameters`
Immutable struct with 28 fields read from TOML:

| Field | Type | Role |
|---|---|---|
| `num_cells` | Int | Number of agents |
| `dt` | Float64 | Euler timestep |
| `total_time` | Float64 | Simulation duration (steps = total_time/dt) |
| `init_width` | Float64 | Width of initial cell placement region |
| `init_height` | Float64 | Height of initial cell placement region |
| `D_angle` | Float64 | Angular diffusion coefficient |
| `D_position` | Float64 | Positional diffusion coefficient |
| `theta_alignment` | Float64 | Restoring torque toward thetabar (tumbling only) |
| `run_to_tumble_rate_min` | Float64 | Min rate r→t (used as constant when function="constant") |
| `run_to_tumble_rate_max` | Float64 | Max rate r→t (sigmoid ceiling) |
| `run_to_tumble_rate_function` | String | `"constant"` or `"sigmoid"` |
| `run_to_tumble_rate_allure` | Float64 | Sigmoid steepness for r→t |
| `tumble_to_run_rate_min` | Float64 | Min rate t→r |
| `tumble_to_run_rate_max` | Float64 | Max rate t→r |
| `tumble_to_run_rate_function` | String | `"constant"` or `"sigmoid"` |
| `tumble_to_run_rate_allure` | Float64 | Sigmoid steepness for t→r |
| `cell_hard_radius` | Float64 | Excluded-volume radius; hard contact at 2×this |
| `cell_soft_radius` | Float64 | Soft contact / neighbour detection radius |
| `maximal_stretch` | Float64 | Max allowed bottom–top distance before hard constraint |
| `run_speed_min` | Float64 | Min run speed (constant dist uses this) |
| `run_speed_max` | Float64 | Max run speed |
| `run_speed_dist` | String | `"constant"`, `"uniform"`, or `"normal"` |
| `top_internal_pull_strength` | Float64 | Spring force pulling top toward bottom |
| `bottom_internal_pull_strength` | Float64 | Spring force pulling bottom toward top |
| `soft_boundary_potential_strength` | Float64 | Wall repulsion magnitude |
| `hard_boundary_velocity_behaviour` | String | `"bounce"`, `"tangential"`, or `"slide"` |
| `cil_intensity` | Float64 | CIL strength from cell neighbours |
| `cil_wall_intensity` | Float64 | CIL strength from walls |

### `DomainSpecs <: Parameters`
`domain_width`, `domain_height` — rectangular 2D domain boundaries.

### `MorsePotential <: Parameters`
`u0`, `v0`, `xi1`, `xi2`, `cutoff` — parameters for force law `F = (u0·exp(-d/ξ1) - v0·exp(-d/ξ2)) × unit_vec`.
Two instances created per simulation:
- **top** — applied at top surface between 2×hard_radius and `top_morse_cutoff`
- **bottom** — applied at bottom surface (input distance is shifted by 2×soft_radius), gated by `bottom_morse_cutoff`

### `CellCollective`
Twelve parallel arrays, one entry per cell:

| Field | Type | Contents |
|---|---|---|
| `bottom` | `Vector{SVector{2,Float64}}` | 2D position at z=0 (main surface) |
| `top` | `Vector{SVector{2,Float64}}` | 2D position at z=z0 (upper surface) |
| `theta` | `Vector{Float64}` | Orientation angle (radians) |
| `thetabar` | `Vector{SVector{2,Float64}}` | Accumulated desired-polarity vector |
| `neighboursbottom` | `Vector{Vector{Int}}` | Indices of bottom-surface neighbours |
| `neighbourstop` | `Vector{Vector{Int}}` | Indices of top-surface neighbours |
| `forcesbottom` | `Vector{SVector{2,Float64}}` | Accumulated forces at bottom |
| `forcestop` | `Vector{SVector{2,Float64}}` | Accumulated forces at top |
| `state` | `Vector{Int}` | 0 = tumbling, 1 = running |
| `state_timer` | `Vector{Float64}` | Countdown to next state switch (Exp(1) distributed) |
| `run_speeds` | `Vector{Float64}` | Per-cell speed when running |
| `hitting_wall` | `Vector{Bool}` | True when cell touches a wall |

### `Solution`
Four time-series, each `Vector` of length `num_steps`:
`bottom`, `top` (positions), `state`, `theta` — each element is a per-cell snapshot.

### `AnalysisResults`
Five `Vector{Float64}` of length `num_steps`:
`mean_speeds_bottom`, `mean_speeds_top`, `group_speeds_bottom`, `group_speeds_top`, `mean_theta`.

---

## Simulation Loop

Defined in [src/simulation_logic.jl](src/simulation_logic.jl). Per timestep, in this exact order:

```
1. compute_hard_interaction_forces!(collective, p)                              # overlap repulsion + neighbour lists + thetabar
2. compute_soft_interaction_forces!(collective, p, morse_top, morse_bottom)     # Morse forces
3. compute_stochastic_forces!(collective, p)
4. compute_state_changes!(collective, p)            # uses global `domain`
5. update_cell_collective!(collective, p)
6. compute_cell_cohesion_forces!(collective, p)
7. compute_soft_wall_forces!(collective, domain, p)
8. apply_hard_wall_boundary_conditions!(collective, domain, p)
9. store_solution!(solution, collective, step, num_cells)
```

`update_cell_collective!` resets forces, neighbours, thetabar, and hitting_wall at the end of each step — these arrays are scratch space rebuilt every timestep.

---

## Physics Details

### Cell Geometry
- `bottom` is the authoritative position. At init, `top = bottom − cell_hard_radius·(cos θ, sin θ)` — so `top` sits *behind* the heading by one hard radius (the cell is drawn as a leading bottom disc trailed by a top disc).
- The spring between them (`compute_cell_cohesion_forces!`) enforces internal integrity:
  - If `|top − bottom| < maximal_stretch`: linear pull `F = pull_strength · |top − bottom| · ûbot→top`, attractive on bottom and equal/opposite on top.
  - If `≥ maximal_stretch`: positions are snapped along their connecting axis so the separation is exactly `maximal_stretch` (not collapsed to a point).

### Hard Interaction Forces (`compute_hard_interaction_forces!`)
O(n²) double loop over all pairs. Position-level corrections + neighbour bookkeeping. For each pair (i, j):

**Bottom surface:**
- `d < 2·hard_radius` → snap both bottoms apart along the connecting axis by `(2·hard_radius − d)/2` each; add to `neighboursbottom`
- `d < 2·soft_radius` (and ≥ 2·hard_radius) → add to `neighboursbottom` only

**Top surface:**
- `d < 2·hard_radius` → snap apart + accumulate `thetabar` from repulsion direction
- `d < 2·soft_radius` → accumulate `thetabar` only (no position change)

### Soft Interaction Forces (`compute_soft_interaction_forces!`)
Runs after hard corrections, so it sees the snapped positions.

**Bottom surface:**
- Morse force into `forcesbottom` if `d < bottom_morse_cutoff`, with `d_eff = d − 2·soft_radius` (can be negative inside the soft zone — see the gotcha on the Morse curve).

**Top surface:**
- Morse force into `forcestop` if `2·hard_radius < d < top_morse_cutoff`, with `d_eff = d − 2·hard_radius`.

### Run-Tumble Dynamics (`compute_state_changes!`)
Rates are x-position-dependent (sigmoid) or constant:
- `compute_run_to_tumble_rate(x, p)` — sigmoid increases toward x=0 (left wall)
- `compute_tumble_to_run_rate(x, p)` — sigmoid also increases toward x=0

State timers count down at the current rate each dt. When timer ≤ 0: toggle state, reset timer to `rand(Exponential(1.0))`.

### CIL (Contact Inhibition of Locomotion)
Applied only to running cells (state=1):

**Cell neighbours:** For each neighbour j in `neighboursbottom[i]`, compute dot product of motion direction with unit vector toward j. If neighbour is in the forward 35° cone:
```
run_to_tumble_rate += cil_intensity * (atan((dot - cos(35°)) / 1e-5) + π/2)
```
This is a near-step-function approximation — the atan argument becomes ±π/2 very sharply at the threshold.

**Walls:** Same logic applied to each wall within `cell_soft_radius`, using 85° threshold and wall outward normals.

### Orientation Update (Tumbling cells, `update_cell_collective!`)
```
dθ = -theta_alignment · sin(θ - atan2(thetabar)) · dt + √(2·D_angle·dt) · N(0,1)
```
Running cells do not rotate via this equation — their θ is fixed during running and only indirectly updated by `thetabar` accumulation.

### Boundary Conditions
- **Soft walls** (`compute_soft_wall_forces!`): linear repulsion `F = strength · (soft_radius - dist)` once bottom or top enters soft zone. Sets `hitting_wall=true`.
- **Hard walls** (`apply_hard_wall_boundary_conditions!`): clamp position to `[hard_radius, domain_size - hard_radius]`. Velocity handled by `hard_boundary_velocity_behaviour`:
  - `"bounce"`: θ ← π - θ (lateral) or θ ← -θ (horizontal)
  - `"tangential"`: remove velocity component normal to wall
  - `"slide"`: no angle change (currently same as default)
- **Periodic** (`apply_periodic_boundary_conditions!`): defined but **not used** in the main loop.

---

## Configuration

### `data/sims/parameters.toml` — Current Values

```toml
# Domain
domain_width = 800
domain_height = 150

# Simulation
num_cells = 100
dt = 0.05
total_time = 400              # → 8000 timesteps

# Initialization region (cells drop uniformly into [0, init_width] × [0, init_height])
init_width = 200
init_height = 150

# Diffusion
D_angle = 0.05
D_position = 0.01

# Orientation
theta_alignment = 3

# State switching (x-dependent sigmoid rates)
run_to_tumble_rate_min = 0.03
run_to_tumble_rate_max = 0.4
run_to_tumble_rate_function = "sigmoid"
run_to_tumble_rate_allure = 0.01
tumble_to_run_rate_min = 0.5
tumble_to_run_rate_max = 2.0
tumble_to_run_rate_function = "sigmoid"
tumble_to_run_rate_allure = 0.01

# Cell geometry
cell_hard_radius = 4
cell_soft_radius = 4
maximal_stretch = 2000

# Cell speed
run_speed_min = 40
run_speed_max = 40
run_speed_distribution = "constant"

# Internal cohesion
top_internal_pull_strength = 1
bottom_internal_pull_strength = 1
internal_cohesion_strength = 1   # present in TOML but not read into struct

# Boundaries
soft_boundary_potential_strength = 2
hard_boundary_velocity_behaviour = "tangential"  # "bounce" | "tangential" | "slide"

# Morse potentials
top_morse_u0 = 100
top_morse_v0 = 60
top_morse_xi1 = 10
top_morse_xi2 = 40
top_morse_cutoff = 50
bottom_morse_u0 = 10
bottom_morse_v0 = 6.1
bottom_morse_xi1 = 10
bottom_morse_xi2 = 60
bottom_morse_cutoff = 60

# CIL
cil_intensity = 0.01
cil_wall_intensity = 0.02
```

Parameter reading is handled by `read_parameters(filename)` in [src/cell_model.jl](src/cell_model.jl), which uses the recursive helper `_toml_value` to locate keys anywhere in the TOML hierarchy. The parameter file does not need to be flat — TOML tables are supported.

---

## Parameter Reference & Impact

This section documents every parameter read from `parameters.toml`, what it controls in the equations of motion, and how the simulation responds to changes. Read alongside [src/cell_model.jl](src/cell_model.jl).

### Domain & Initialisation

| Parameter | Role | Impact of increasing |
|---|---|---|
| `domain_width` | Width of the rectangular domain. Used by every wall check (soft potential, hard clamp, CIL wall cone, sigmoid rate's x-axis). | More room to spread; sigmoid rate's "left vs right" gradient (rates depend on `x`) is stretched horizontally, so the run/tumble bias becomes shallower across the domain. |
| `domain_height` | Height of the rectangular domain. | More vertical room. With CIL walls disabled by geometry, you get less wall-induced tumbling. |
| `init_width`, `init_height` | Initial seeding rectangle. Cells are placed uniformly in `[0, init_width] × [0, init_height]`. | Larger → sparser start, less initial overlap (so hard-overlap snap events on step 1 are rarer); smaller → denser cluster, immediate jamming and a strong outward Morse "explosion". |

**Tuning note:** if `num_cells · (2·cell_hard_radius)² > init_width · init_height` the seeding region is over-saturated and cells start with many overlaps; the first few steps are dominated by hard-snap corrections.

### Simulation Control

| Parameter | Role | Impact |
|---|---|---|
| `num_cells` | Number of agents. | Quadratic cost in pair loops (`compute_hard_interaction_forces!`, `compute_soft_interaction_forces!`). Doubling cells ≈ 4× wall-clock per step. Above ~500 the O(n²) loops dominate. |
| `dt` | Forward-Euler timestep. Multiplies *every* force-to-position step (`bottom .+= dt·forces`), noise (`√(2·D·dt)`), state-timer decrement, run-speed displacement, and the orientation update. | Larger → faster but unstable: stiff Morse forces (`u0=100`) easily overshoot, cells tunnel through neighbours, the integrator blows up. Smaller → smoother but linearly slower. Empirical sweet spot: `dt ≤ xi1/(u0·dt) · safety` for the dominant Morse; current `dt=0.05` is conservative for `top_morse_u0=100`. |
| `total_time` | Simulated duration. `num_steps = total_time / dt`. | Linear cost in wall time and memory: `Solution` allocates `num_steps × num_cells × (2 positions + state + theta)`. |

### Stochastic Diffusion

| Parameter | Role | Impact |
|---|---|---|
| `D_position` | Translational noise. Each step adds `xi = √(2·D_position·dt)·𝒩(0, I₂)` to **both** bottom and top forces (same draw applied to both surfaces — this is intentional in [cell_model.jl:208-210](src/cell_model.jl#L208-L210)). | Higher → more jitter, diffusive component dominates ballistic runs; the cluster melts faster. Set to ~0 for clean deterministic CIL studies. |
| `D_angle` | Rotational noise on **tumbling cells only**. Enters `update_cell_collective!` as `√(2·D_angle·dt)·𝒩(0,1)` added to θ. | Higher → tumbling cells reorient more chaotically, less alignment with the polarity vector `thetabar`; the population becomes less polar. |

### Orientation Alignment

| Parameter | Role | Impact |
|---|---|---|
| `theta_alignment` | Restoring torque toward `atan2(thetabar)`. Applied to tumbling cells only: `dθ = −theta_alignment · sin(θ − θ̄) · dt`. | Higher → tumbling cells snap toward the locally accumulated polarity (contact direction, wall normals) very quickly. Together with `D_angle` it controls the alignment-vs-noise balance. Effective relaxation time ≈ `1 / theta_alignment`. |

### Run–Tumble State Machine

Each cell carries a `state ∈ {0=tumble, 1=run}` and an exponential timer that decrements at the current rate:
- If running: `timer −= dt · run_to_tumble_rate(x)` (boosted by CIL — see below)
- If tumbling: `timer −= dt · tumble_to_run_rate(x)`

Rate functions (over the cell's x-position):
- `"constant"` → returns `*_min` (regardless of `_max`)
- `"sigmoid"` → for `tumble_to_run`: `min + (max − min) / (1 + exp(−allure·x))` — *increases* with x
- `"sigmoid"` → for `run_to_tumble`: `max + (min − max) / (1 + exp(−allure·x))` — *decreases* with x

| Parameter | Role | Impact |
|---|---|---|
| `run_to_tumble_rate_min` / `_max` | Bounds on r→t rate. Constant mode uses `_min`; sigmoid sweeps between them. | Higher rate → shorter mean run duration (`E[run] = 1/r→t`). Set `max` ≪ `min` to make cells reach the right wall and tumble there; set `min` < `max` to make them tumble at the left wall. |
| `run_to_tumble_rate_function` | `"constant"` or `"sigmoid"`. | `"sigmoid"` makes runs x-dependent (gradient along the long axis); `"constant"` removes any spatial bias. |
| `run_to_tumble_rate_allure` | Sigmoid steepness. `0.01` gives a wide ramp across an 800-wide domain; `1.0` would saturate within a few units of x=0. | Higher → sharper switch at x=0; effectively makes a step function. |
| `tumble_to_run_rate_min` / `_max` | Bounds on t→r rate. | Higher → shorter tumbles (`E[tumble] = 1/t→r`). The **fraction running** in steady state ≈ `t→r / (t→r + r→t)`. With current defaults (`r→t∈[0.03,0.4]`, `t→r∈[0.5,2.0]`) most cells run most of the time. |
| `tumble_to_run_rate_function`, `_allure` | Same shape as their run-to-tumble counterparts. | Combine with `run_to_tumble_*` to engineer spatial preferences (e.g. "cells tend to run rightward and tumble at the right wall"). |

**Initial state distribution** uses a different ratio: `P(start = run) = tumble_to_run_rate_min / (tumble_to_run_rate_min + run_to_tumble_rate_max)`. The asymmetry (`_max` for r→t vs `_min` for t→r) is hard-coded in [cell_model.jl:98-99](src/cell_model.jl#L98-L99) — be aware when interpreting transients.

### Cell Geometry & Cohesion

| Parameter | Role | Impact |
|---|---|---|
| `cell_hard_radius` | Excluded-volume radius. Pair overlap test `d < 2·hard_radius` triggers position snap. Also used as Morse `d_eff` shift on the top surface, in the soft/hard wall clamps, and as the initial top-vs-bottom offset. | Bigger cells → more overlap snaps per step (jamming), higher effective Morse minimum on the top surface (since `d_eff = d − 2·hard_radius`), and the hard-wall clamp band gets thicker. |
| `cell_soft_radius` | Soft-contact / neighbour-detection radius. Used by the CIL wall trigger, the soft-wall potential range, and the bottom-surface Morse `d_eff` shift (`d_eff = d − 2·soft_radius`). | Wider neighbour halo → more CIL pairs counted per running cell → more tumbling. Note: if `soft_radius = hard_radius` (current default!) the "soft contact" band has zero width — only direct overlaps register as neighbours, and the bottom-surface Morse uses the same shift as the top would. |
| `maximal_stretch` | Cap on bottom–top separation. Below this, a linear pull-spring acts; at/above, both endpoints are snapped along their connecting axis until separation = `maximal_stretch`. | With the current default of `2000` (≫ any plausible separation), the cap is never hit and only the linear spring matters. Lower it to enforce a "cell length" constraint. |
| `top_internal_pull_strength` | Magnitude of the linear spring pulling top toward bottom. Force on top is `−strength · |Δ| · û_b→t`. Note this is *linear in distance*, not Hookean around an equilibrium — there's no rest length, so the spring always pulls together. | Higher → top snaps onto bottom faster; the pseudo-3D "extension" of the cell shrinks. Combined with run speed (which only displaces `bottom`), this controls how far top "drags" behind a moving cell. |
| `bottom_internal_pull_strength` | Same spring on the bottom side. | Higher → bottom is pulled toward top each step; if asymmetric with `top_internal_pull_strength`, the pair drifts as a whole. |

### Run Speed Distribution

| Parameter | Role | Impact |
|---|---|---|
| `run_speed_min`, `run_speed_max` | Bounds on per-cell run speed. Run displacement per step = `speed · dt · (cos θ, sin θ)` applied to `bottom` only (top follows via the cohesion spring). | Doubling speed roughly doubles the effective Péclet number (ballistic-vs-diffusive). With `dt=0.05`, a speed of 40 moves a cell 2 units/step — well below `cell_hard_radius=4`, so no tunnelling. |
| `run_speed_distribution` | `"constant"` (uses `min`), `"uniform"` (in `[min, max]`), or `"normal"` (μ=mean(min,max), σ=(max−min)/6). | Heterogeneous populations sort by speed: fast cells lead and slow cells trail in any directional bias. |

### Cell–Cell Morse Potentials

Force law: `F = (u0·exp(−d_eff/ξ1) − v0·exp(−d_eff/ξ2)) · unit_vec`. With `ξ1 < ξ2`, the first term dominates at short range (repulsion) and the second at long range (attraction) — the canonical Morse shape with a single minimum.

For each surface (top, bottom) the same five parameters control the curve, but with different gating:
- **Top:** active for `2·hard_radius < d < top_morse_cutoff`, with `d_eff = d − 2·hard_radius`.
- **Bottom:** active for `d < bottom_morse_cutoff`, with `d_eff = d − 2·soft_radius`. `d_eff` can be **negative** inside the soft band; with `ξ1=10` and `u0=10` the repulsion blows up but stays finite (clipped by the `eps()` guard at exactly zero distance).

| Parameter | Role | Impact |
|---|---|---|
| `top_morse_u0`, `bottom_morse_u0` | Short-range repulsion amplitude. | Higher → stiffer Morse "wall". Pair this with smaller `dt` or the integrator overshoots. |
| `top_morse_v0`, `bottom_morse_v0` | Long-range attraction amplitude. | Higher → stronger cohesion / clustering. The minimum well depth ≈ `v0·(ξ2/ξ1 − 1)` for typical parameters. |
| `top_morse_xi1`, `bottom_morse_xi1` | Repulsion length scale. | Larger → softer repulsion, lower barrier to overlap; minimum shifts outward. |
| `top_morse_xi2`, `bottom_morse_xi2` | Attraction length scale. | Larger → longer-range attraction, denser aggregates. |
| `top_morse_cutoff`, `bottom_morse_cutoff` | Hard cutoff distance. | Below the cutoff: Morse active. Above: zero. Smaller cutoffs save compute (already-cheap) and prevent the long tail of attraction from over-bunching distant cells. |

**Locating the well minimum** (helpful for tuning): solve `dF/d_eff = 0` → `d_eff* = (ξ1·ξ2)/(ξ2 − ξ1) · ln((u0·ξ2)/(v0·ξ1))`. For `top` defaults `(u0=100, v0=60, ξ1=10, ξ2=40)`: `d_eff* ≈ (400/30)·ln(100·40/(60·10)) ≈ 13.3·1.89 ≈ 25` ⇒ preferred *top* separation ≈ `25 + 2·hard_radius = 33`.

### Boundary Potential

| Parameter | Role | Impact |
|---|---|---|
| `soft_boundary_potential_strength` | Linear-spring wall potential. When a bottom or top endpoint enters within `cell_soft_radius` of a wall, force = `strength · (cell_soft_radius − distance_to_wall)` along the inward normal. Only bottom-surface wall contacts set `hitting_wall = true` and contribute to `thetabar`. | Higher → cells bounce off walls before reaching them. Too low (<<1) and the soft layer is ineffective; hard clamp does most of the work. Too high and you get an over-stiff repulsion that needs a smaller `dt`. |
| `hard_boundary_velocity_behaviour` | What happens to θ when a cell is clamped to the wall: `"bounce"` reflects θ (π − θ on lateral walls, −θ on horizontal), `"tangential"` removes the normal velocity component (cell slides along the wall), `"slide"` (or any unrecognised value) keeps θ unchanged — position is still clamped. | `"bounce"` produces ballistic-billiard reflections; `"tangential"` lets cells flow along walls (good for migration corridors); `"slide"` plus high `theta_alignment` lets the polarity machinery handle reorientation alone. |

### Contact Inhibition of Locomotion (CIL)

Applied only to *running* cells, inside `compute_state_changes!`. For each neighbour `j ∈ neighboursbottom[i]` (cell list built by `compute_hard_interaction_forces!`), if `j` lies in a forward cone:

```
run_to_tumble_rate += cil_intensity · (atan((dot(speed_dir, û_i→j) − cos(35°)) / 1e-5) + π/2)
```

The `atan(·/1e-5)` is a near-step approximation: contribution ≈ 0 outside the cone, ≈ `π · cil_intensity` inside. So each forward-cone neighbour adds roughly `π · cil_intensity` to the r→t rate.

| Parameter | Role | Impact |
|---|---|---|
| `cil_intensity` | Per-neighbour boost to r→t when that neighbour is in the forward 35° cone. | Higher → running cells tumble almost immediately on contact. Set to 0 to disable CIL between cells. With several neighbours, the effective r→t can dwarf the baseline (`run_to_tumble_rate_min/max`). |
| `cil_wall_intensity` | Same mechanism vs walls, using an 85° cone and the outward wall normals. Triggered only when the bottom endpoint is within `cell_soft_radius` of a wall. | Higher → cells tumble on wall contact (combined with `hard_boundary_velocity_behaviour`, controls how cells "escape" walls). |

### Cross-couplings worth remembering

- **`dt` × Morse `u0`** — stiff Morse + large `dt` = overshoot and explosion.
- **`cell_soft_radius` × CIL** — wider soft halo → more neighbours per cell → larger CIL boost.
- **`theta_alignment` × `D_angle`** — controls polar order during tumbling.
- **Sigmoid rates × `domain_width`** — `allure` should scale roughly as `1/domain_width` to keep the gradient shape comparable across domain sizes.
- **`run_speed × dt × cell_hard_radius`** — for `speed·dt < hard_radius` you're safe from tunnelling; otherwise pairs can pass through each other between snap steps.

---

## Scripts

| Script | Purpose |
|---|---|
| [scripts/run_sim_csv_save.jl](scripts/run_sim_csv_save.jl) | Run single simulation, save `solution.csv` (columns: time, cell_id, bottom_x, bottom_y, top_x, top_y) |
| [scripts/simple_test.jl](scripts/simple_test.jl) | Run single simulation, render `animation.mp4` via GLMakie |
| [scripts/sweep_logic.jl](scripts/sweep_logic.jl) | Single-parameter sweep; reads `sweep_parameters.toml`; saves `analysis_<param>.csv` |
| [scripts/double_sweep_logic.jl](scripts/double_sweep_logic.jl) | Two-parameter grid sweep; reads `double_sweep_parameters.toml` |
| [scripts/animation_logic.jl](scripts/animation_logic.jl) | `animate_solution(solution, p, domain)` — draws bottom/top circles (coloured by state), hard-radius circles, θ arrows; records MP4 at 30 fps |
| [scripts/plotting_logic.jl](scripts/plotting_logic.jl) | Plotting utilities |

All scripts start with `using DrWatson; @quickactivate "Cell_SE_Model"` which activates the project and enables path helpers.

---

## Analysis Functions

Defined in [src/run_analysis.jl](src/run_analysis.jl):

- **`compute_mean_speed`**: per step, average `‖pos(t) - pos(t-1)‖ / dt` over all cells (both surfaces)
- **`compute_group_speed`**: per step, `‖mean(pos(t)) - mean(pos(t-1))‖ / dt` — centre-of-mass speed
- **`compute_mean_theta`**: arithmetic mean of θ per step (not circular mean — note the naming)

Three entry points: `run_analysis_from_simulation()`, `run_analysis_from_csv()`, `run_analysis_from_solution(solution, p)`.

---

## Parameter Sweeps

Single sweep (`sweep_parameters.toml`): specify one parameter name, start, stop, step, and number of repetitions. The script deep-copies the base parameter dict, overrides the target field, runs the simulation, computes analysis, and aggregates statistics.

Double sweep (`double_sweep_parameters.toml`): same pattern but with two parameters forming a grid.

---

## Dependencies

| Package | Version | Role |
|---|---|---|
| DrWatson | 2.19.1 | `datadir`, `srcdir`, `@quickactivate`, reproducibility |
| StaticArrays | — | `SVector{2,Float64}` for all 2D vectors (performance) |
| Distributions | — | `Exponential`, `Uniform`, `Normal` |
| TOML | stdlib | Parameter file parsing |
| LinearAlgebra | stdlib | `norm`, `dot` |
| CellListMap | — | Imported in cell_model.jl but **not yet used** — interaction loop is O(n²) |
| Colors | — | Cell colouring in animations |
| FFMPEG | 0.4.5 | MP4 encoding |
| GLMakie | — | Interactive 2D animation |

---

## Known Gotchas

- **`compute_state_changes!` uses a global `domain`** — the variable is not passed as an argument but accessed from the enclosing scope in `simulation_loop`. This works but is fragile; be careful when calling this function in isolation.
- **`apply_periodic_boundary_conditions!` is dead code** — defined but never called; hard walls are always active.
- **`mean_theta` is not a circular mean** — arithmetic average of angles wraps incorrectly near ±π. Use `atan(mean(sin(θ)), mean(cos(θ)))` for a proper circular mean.
- **CellListMap imported but unused** — the interaction loops are plain O(n²) double `for`s (both hard and soft passes). For `num_cells` much larger than ~500, this becomes the bottleneck.
- **`internal_cohesion_strength` in TOML is unused** — present in `parameters.toml` but not read into `SimulationParameters`.
- **`run_speed_distribution` vs `run_speed_dist`** — TOML key is `run_speed_distribution` but the struct field is `run_speed_dist`. The `_toml_value` recursive lookup handles this transparently.
- **Stochastic noise is shared between bottom and top** — `compute_stochastic_forces!` draws one `xi` per cell and adds the *same* vector to both `forcesbottom[i]` and `forcestop[i]`. The two surfaces of a single cell are perfectly correlated in their Brownian kick; they're not independent.
- **Hard and soft interaction passes recompute distances independently** — after the split, `compute_soft_interaction_forces!` sees positions already corrected by `compute_hard_interaction_forces!`. Pairs that overlapped this step thus get a Morse force evaluated at the post-snap distance (≈ `2·hard_radius` for the top surface), not the pre-snap distance the original combined loop used. For non-stiff Morse parameters this is a minor difference.
- **`top` is initialised behind the heading** — `top[i] = bottom[i] − cell_hard_radius · (cos θ, sin θ)`. The top disc trails the bottom along θ; the cohesion spring is what pulls them together as the cell moves.
- **`maximal_stretch` does not snap to zero** — when exceeded, the two endpoints are moved along their connecting axis until separation = `maximal_stretch`, not until they coincide.
- **Initial state ratio is asymmetric** — `P(start = running) = tumble_to_run_rate_min / (tumble_to_run_rate_min + run_to_tumble_rate_max)` (note: `_min` for one side, `_max` for the other). With the current defaults this gives `0.5 / (0.5 + 0.4) ≈ 0.56`.
- **`compute_run_to_tumble_rate` sigmoid is decreasing in x** — opposite shape to `compute_tumble_to_run_rate`. If you set both to `"sigmoid"`, runs become long on the right side of the domain and short on the left; double-check direction before interpreting sweeps.

---

## How to Run

```julia
# Activate the project first (in Julia REPL)
using Pkg; Pkg.activate(".")

# Single simulation → CSV
include("scripts/run_sim_csv_save.jl")

# Single simulation → animation
include("scripts/simple_test.jl")

# Parameter sweep
include("scripts/sweep_logic.jl")
```

Or use `julia --project=. scripts/run_sim_csv_save.jl` from the shell.
