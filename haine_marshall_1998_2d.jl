using Oceananigans
using Oceananigans.TurbulenceClosures

#=
Haine & Marshall (1998) 2D / "2.5D" mixed-layer setup in Oceananigans.

This script follows the paper's experiment 2 as a practical Oceananigans
starting point:

  - nonhydrostatic Boussinesq dynamics on an f-plane
  - no x-dependence (implemented with Nx = 1 and periodic x)
  - along-front velocity u is retained, allowing thermal-wind shear
  - free-slip boundaries at the top, bottom, and side walls
  - initially motionless, uniformly stratified fluid
  - surface buoyancy loss varying across-channel as

      B(y) = B_half * (tanh(-2 * (y - Ly / 2) / Lf) + 1)

Paper details used here:
  - f   = 1.0e-4 s^-1
  - N   = 8.37e-4 s^-1
  - Ly  = 30 km
  - nominal along-front length = 50 km
  - forcing width Lf = 10 km
  - B_half = 1.96e-7 m^2 s^-3

The original paper used 250 m horizontal spacing and vertical spacing that
varied from about 40 m near the surface to about 400 m near the bottom.
We approximate that here with a stretched z grid.

Install dependencies with:

    julia -e 'using Pkg; Pkg.add("Oceananigans")'

Run with:

    julia haine_marshall_1998_2d.jl
=#

const Nx = 1
const Ny = 120                     # 250 m across-channel resolution over 30 km
const Nz = 13                      # ~40 m to ~400 m stretched vertical spacing

const Lx = 50e3                    # m
const Ly = 30e3                    # m

const N = 8.37e-4                  # s^-1
const N² = N^2
const f = 1.0e-4                   # s^-1

const B_half = 1.96e-7             # m^2 s^-3
const Lf = 10e3                    # m

const νh = 5.0                     # m^2 s^-1
const νz = 0.02                    # m^2 s^-1
const κh = 5.0                     # m^2 s^-1
const κz = 0.02                    # m^2 s^-1

parse_env_float(name, default) = haskey(ENV, name) ? parse(Float64, ENV[name]) : default

const simulation_days = parse_env_float("HM98_SIMULATION_DAYS", 8.0)
const stop_time = simulation_days * 86400

# Set to `nothing` to keep cooling on for the whole run.
const cooling_stop_days = haskey(ENV, "HM98_COOLING_STOP_DAYS") ? parse(Float64, ENV["HM98_COOLING_STOP_DAYS"]) : nothing
const cooling_stop_time = isnothing(cooling_stop_days) ? nothing : cooling_stop_days * 86400

const initial_dt = parse_env_float("HM98_INITIAL_DT", 60.0)
const max_dt = parse_env_float("HM98_MAX_DT", 600.0)

"""
    stretched_z_faces(Nz; Δz_top=40, Δz_bottom=400)

Return `Nz + 1` cell faces with geometric stretching such that the top cell is
`Δz_top` thick and the bottom cell is `Δz_bottom` thick.
"""
function stretched_z_faces(Nz; Δz_top = 40.0, Δz_bottom = 400.0)
    r = (Δz_bottom / Δz_top)^(1 / (Nz - 1))
    Δz = Δz_top .* r .^ (0:Nz-1)
    z_faces = zeros(Nz + 1)

    for k in 1:Nz
        z_faces[k + 1] = z_faces[k] + Δz[Nz - k + 1]
    end

    return z_faces .- maximum(z_faces)
end

const z_faces = stretched_z_faces(Nz)
const Lz = -minimum(z_faces)

grid = RectilinearGrid(
    size = (Nx, Ny, Nz),
    x = (0, Lx),
    y = (0, Ly),
    z = z_faces,
    topology = (Periodic, Bounded, Bounded)
)

closure = (
    HorizontalScalarDiffusivity(ν = νh, κ = κh),
    VerticalScalarDiffusivity(ν = νz, κ = κz)
)
buoyancy = BuoyancyTracer()
coriolis = FPlane(f = f)

free_slip = FieldBoundaryConditions(
    bottom = GradientBoundaryCondition(0.0),
    top = GradientBoundaryCondition(0.0)
)

side_walls = FieldBoundaryConditions(
    south = GradientBoundaryCondition(0.0),
    north = GradientBoundaryCondition(0.0),
    bottom = GradientBoundaryCondition(0.0),
    top = GradientBoundaryCondition(0.0)
)

@inline function surface_buoyancy_flux(x, y, t, p)
    forcing = p.B_half * (tanh(-2 * (y - p.Ly / 2) / p.Lf) + 1)
    return isnothing(p.cooling_stop_time) || t <= p.cooling_stop_time ? forcing : 0.0
end

b_bcs = FieldBoundaryConditions(
    south = GradientBoundaryCondition(0.0),
    north = GradientBoundaryCondition(0.0),
    bottom = GradientBoundaryCondition(0.0),
    top = FluxBoundaryCondition(surface_buoyancy_flux,
                                parameters = (; B_half = B_half,
                                                Ly = Ly,
                                                Lf = Lf,
                                                cooling_stop_time = cooling_stop_time))
)

background_buoyancy(x, y, z, t) = N² * z

model = NonhydrostaticModel(grid;
    tracers = :b,
    buoyancy,
    coriolis,
    closure,
    boundary_conditions = (
        u = side_walls,
        v = free_slip,
        b = b_bcs
    ),
    background_fields = (b = BackgroundField(background_buoyancy),)
)

u_initial(x, y, z) = 0.0
v_initial(x, y, z) = 1e-4 * randn()
w_initial(x, y, z) = 1e-4 * randn()
b_initial(x, y, z) = 1e-6 * randn()

set!(model, u = u_initial, v = v_initial, w = w_initial, b = b_initial)

simulation = Simulation(model, Δt = initial_dt, stop_time = stop_time)
wizard = TimeStepWizard(cfl = 0.4, max_change = 1.1, max_Δt = max_dt)
simulation.callbacks[:wizard] = Callback(wizard, IterationInterval(10))

function progress(sim)
    model = sim.model
    u, v, w = model.velocities
    b = model.tracers.b

    umax = maximum(abs, interior(u))
    vmax = maximum(abs, interior(v))
    wmax = maximum(abs, interior(w))
    bmax = maximum(abs, interior(b))

    @info "iteration=$(iteration(sim)) time_days=$(round(time(sim) / 86400, digits=3)) Δt=$(round(sim.Δt, digits=2)) max|u|=$(umax) max|v|=$(vmax) max|w|=$(wmax) max|b'|=$(bmax)"
end

simulation.callbacks[:progress] = Callback(progress, IterationInterval(50))

@info "Starting Haine-Marshall 1998 2D setup"
@info "grid = $(Nx) x $(Ny) x $(Nz)"
@info "domain = Lx=$(Lx / 1e3) km, Ly=$(Ly / 1e3) km, Lz=$(round(Lz, digits=1)) m"
@info "N = $(N) s^-1, f = $(f) s^-1"
@info "B_half = $(B_half) m^2 s^-3, Lf = $(Lf / 1e3) km"
cooling_stop_label = isnothing(cooling_stop_time) ? "none" : string(cooling_stop_time / 86400)
@info "cooling_stop_time_days = $(cooling_stop_label)"

run!(simulation)

@info "Finished at day $(time(simulation) / 86400) after $(iteration(simulation)) iterations."
