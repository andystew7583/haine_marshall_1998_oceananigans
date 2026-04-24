module HM98Common

using Oceananigans
using Oceananigans.TurbulenceClosures
using Oceananigans.OutputWriters: write_output!
using CairoMakie
using Statistics

export HM98Experiment,
       HM98Runtime,
       experiment_spec,
       load_runtime,
       run_case,
       plot_case,
       default_output_path,
       default_movie_path,
       default_xy_movie_path

# The HM98 scripts all use Float64 on CPU for consistency with the earlier runs.
const FT = Float64
const seconds_per_day = FT(86400)

# Shared domain scales used in the HM98 experiments.
const Lx_default = FT(50e3)
const Ly_default = FT(30e3)
const total_depth_default = FT(2100.0)

# Reference physical parameters from Haine & Marshall (1998).
const N = FT(8.37e-4)
const N² = N^2
const f = FT(1.0e-4)
const Q₀ = f * N²
const Lf_default = FT(10e3)

# Diffusivities/viscosities currently used in the refined 2D runs.
const νh_default = FT(1.25)
const νz_default = FT(0.005)
const κh_default = FT(1.25)
const κz_default = FT(0.005)
const νh_hm98 = FT(5.0)
const νz_hm98 = FT(0.02)
const κh_hm98 = FT(5.0)
const κz_hm98 = FT(0.02)

"""
    HM98Experiment

Container for case-specific choices: forcing amplitude, dimensionality, and plotting style.
"""
Base.@kwdef struct HM98Experiment
    id::Int
    name::String
    env_prefix::String
    output_prefix::String
    dimensions::Int
    forcing::Symbol
    B_half::FT
    Lf::FT = Lf_default
    cooling_stop_days::Union{Nothing, FT} = nothing
    plot_kind::Symbol
    pv_limits::Union{Nothing, Tuple{FT, FT}} = nothing
    contours_on_all_panels::Bool = false
end

"""
    HM98Runtime

Runtime and grid controls. Most values can be overridden through `HM98_EXP*_...`
environment variables.
"""
Base.@kwdef struct HM98Runtime
    simulation_days::FT = FT(10.0)
    snapshot_interval_days::FT = FT(0.05)
    initial_dt::FT = FT(60.0)
    max_dt::FT = FT(600.0)
    plot_depth::FT = FT(1000.0)
    Δx_target::FT = FT(62.5)
    Δy_target::FT = FT(62.5)
    use_stretched_vertical_grid::Bool = false
    Δz_uniform::FT = FT(0.0)
    Δz_upper::FT = FT(12.5)
    Δz_lower::FT = FT(50.0)
    upper_layer_depth::FT = FT(1000.0)
    surface_Δz::FT = FT(10.0)
    mid_Δz::FT = FT(30.0)
    deep_Δz::FT = FT(100.0)
    deep_transition_depth::FT = FT(2000.0)
    total_depth::FT = total_depth_default
    Lx::FT = Lx_default
    Ly::FT = Ly_default
    νh::FT = νh_default
    νz::FT = νz_default
    κh::FT = κh_default
    κz::FT = κz_default
    pv_slice_depth::FT = FT(671.0)
    buoyancy_slice_depth::FT = FT(65.0)
end

parse_env_float(name, default) = haskey(ENV, name) ? parse(Float64, ENV[name]) : default

"""
    experiment_spec(id)

Return the metadata for experiments 1–4. Experiments 3 and 4 are fully 3D in HM98.
Experiment 4 follows Table 1 with doubled surface buoyancy loss relative to experiment 3.
"""
function experiment_spec(id::Int)
    if id == 1
        return HM98Experiment(id = 1,
                              name = "experiment 1",
                              env_prefix = "HM98_EXP1",
                              output_prefix = "hm98_experiment1",
                              dimensions = 2,
                              forcing = :uniform,
                              B_half = FT(1.96e-7),
                              plot_kind = :diagnostics,
                              pv_limits = (-1.0, 1.0),
                              contours_on_all_panels = true)
    elseif id == 2
        return HM98Experiment(id = 2,
                              name = "experiment 2",
                              env_prefix = "HM98_EXP2",
                              output_prefix = "hm98_experiment2",
                              dimensions = 2,
                              forcing = :patterned,
                              B_half = FT(1.96e-7),
                              plot_kind = :diagnostics,
                              pv_limits = (-1.0, 1.0),
                              contours_on_all_panels = true)
    elseif id == 3
        return HM98Experiment(id = 3,
                              name = "experiment 3",
                              env_prefix = "HM98_EXP3",
                              output_prefix = "hm98_experiment3",
                              dimensions = 3,
                              forcing = :patterned,
                              B_half = FT(1.96e-7),
                              cooling_stop_days = FT(5.0),
                              plot_kind = :diagnostics,
                              pv_limits = (-1.0, 1.0),
                              contours_on_all_panels = true)
    elseif id == 4
        return HM98Experiment(id = 4,
                              name = "experiment 4",
                              env_prefix = "HM98_EXP4",
                              output_prefix = "hm98_experiment4",
                              dimensions = 3,
                              forcing = :patterned,
                              B_half = FT(3.92e-7),
                              plot_kind = :diagnostics,
                              pv_limits = (-1.0, 1.0),
                              contours_on_all_panels = true)
    end

    error("Unsupported HM98 experiment id: $id")
end

default_snapshot_interval_days(experiment::HM98Experiment) = experiment.id == 1 ? 0.1 : 0.05
default_dx(experiment::HM98Experiment) = experiment.id == 4 ? 250.0 : 62.5
default_dy(experiment::HM98Experiment) = experiment.id == 4 ? 250.0 : 62.5
default_dz_uniform(experiment::HM98Experiment) = 0.0
default_use_stretched_vertical_grid(experiment::HM98Experiment) = experiment.id == 4
default_νh(experiment::HM98Experiment) = experiment.id == 4 ? νh_hm98 : νh_default
default_νz(experiment::HM98Experiment) = experiment.id == 4 ? νz_hm98 : νz_default
default_κh(experiment::HM98Experiment) = experiment.id == 4 ? κh_hm98 : κh_default
default_κz(experiment::HM98Experiment) = experiment.id == 4 ? κz_hm98 : κz_default

"""
    load_runtime(experiment)

Read environment overrides for the requested experiment. Experiment 4 defaults to the
coarser 3D mesh requested here: 250 m horizontal spacing with a stretched vertical grid,
and the original HM98 Laplacian diffusivities/viscosities at that spacing.
"""
function load_runtime(experiment::HM98Experiment)
    prefix = experiment.env_prefix

    return HM98Runtime(
        simulation_days = FT(parse_env_float("$(prefix)_SIMULATION_DAYS", 10.0)),
        snapshot_interval_days = FT(parse_env_float("$(prefix)_SNAPSHOT_INTERVAL_DAYS", default_snapshot_interval_days(experiment))),
        initial_dt = FT(parse_env_float("$(prefix)_INITIAL_DT", 60.0)),
        max_dt = FT(parse_env_float("$(prefix)_MAX_DT", 600.0)),
        plot_depth = FT(parse_env_float("$(prefix)_PLOT_DEPTH", 1000.0)),
        Δx_target = FT(parse_env_float("$(prefix)_DX", default_dx(experiment))),
        Δy_target = FT(parse_env_float("$(prefix)_DY", default_dy(experiment))),
        use_stretched_vertical_grid = get(ENV, "$(prefix)_STRETCHED_Z", default_use_stretched_vertical_grid(experiment) ? "true" : "false") == "true",
        Δz_uniform = FT(parse_env_float("$(prefix)_DZ", default_dz_uniform(experiment))),
        Δz_upper = FT(parse_env_float("$(prefix)_DZ_UPPER", 12.5)),
        Δz_lower = FT(parse_env_float("$(prefix)_DZ_LOWER", 50.0)),
        upper_layer_depth = FT(parse_env_float("$(prefix)_UPPER_LAYER_DEPTH", 1000.0)),
        surface_Δz = FT(parse_env_float("$(prefix)_SURFACE_DZ", 10.0)),
        mid_Δz = FT(parse_env_float("$(prefix)_MID_DZ", 30.0)),
        deep_Δz = FT(parse_env_float("$(prefix)_DEEP_DZ", 100.0)),
        deep_transition_depth = FT(parse_env_float("$(prefix)_DEEP_DEPTH", 2000.0)),
        total_depth = FT(parse_env_float("$(prefix)_TOTAL_DEPTH", total_depth_default)),
        Lx = FT(parse_env_float("$(prefix)_LX", Lx_default)),
        Ly = FT(parse_env_float("$(prefix)_LY", Ly_default)),
        νh = FT(parse_env_float("$(prefix)_NUH", default_νh(experiment))),
        νz = FT(parse_env_float("$(prefix)_NUZ", default_νz(experiment))),
        κh = FT(parse_env_float("$(prefix)_KAPPAH", default_κh(experiment))),
        κz = FT(parse_env_float("$(prefix)_KAPPAZ", default_κz(experiment))),
        pv_slice_depth = FT(parse_env_float("$(prefix)_PV_SLICE_DEPTH", 671.0)),
        buoyancy_slice_depth = FT(parse_env_float("$(prefix)_BUOYANCY_SLICE_DEPTH", 65.0)),
    )
end

default_output_path(project_dir, experiment::HM98Experiment) =
    joinpath(project_dir, "outputs", "$(experiment.output_prefix)_fields.jld2")

default_yz_output_path(project_dir, experiment::HM98Experiment) =
    joinpath(project_dir, "outputs", "$(experiment.output_prefix)_yz_slice.jld2")

default_xy_output_path(project_dir, experiment::HM98Experiment) =
    joinpath(project_dir, "outputs", "$(experiment.output_prefix)_xy_slice.jld2")

default_xy_buoyancy_output_path(project_dir, experiment::HM98Experiment) =
    joinpath(project_dir, "outputs", "$(experiment.output_prefix)_xy_buoyancy_slice.jld2")

function default_movie_path(project_dir, experiment::HM98Experiment)
    suffix = experiment.plot_kind == :buoyancy ? "buoyancy" : "diagnostics"
    return joinpath(project_dir, "outputs", "$(experiment.output_prefix)_$(suffix).mp4")
end

default_xy_movie_path(project_dir, experiment::HM98Experiment) =
    joinpath(project_dir, "outputs", "$(experiment.output_prefix)_pv_xy.mp4")

default_xy_buoyancy_movie_path(project_dir, experiment::HM98Experiment) =
    joinpath(project_dir, "outputs", "$(experiment.output_prefix)_buoyancy_xy.mp4")

function uniform_z_faces(runtime::HM98Runtime)
    Nz = Int(round(runtime.total_depth / runtime.Δz_uniform))
    return collect(range(-runtime.total_depth, FT(0); length = Nz + 1))
end

function piecewise_z_faces(runtime::HM98Runtime)
    n_upper = Int(round(runtime.upper_layer_depth / runtime.Δz_upper))
    n_lower = Int(round((runtime.total_depth - runtime.upper_layer_depth) / runtime.Δz_lower))

    upper_faces = collect(range(-runtime.upper_layer_depth, FT(0); length = n_upper + 1))
    lower_faces = collect(range(-runtime.total_depth, -runtime.upper_layer_depth; length = n_lower + 1))

    return vcat(lower_faces, upper_faces[2:end])
end

"""
Construct a smoothly stretched vertical mesh using target cell thicknesses at the
surface, 1000 m, and 2000 m depths. The spacing is linearly interpolated between
those depth landmarks and held fixed below the deepest control point.
"""
function stretched_z_faces(runtime::HM98Runtime)
    control_depths = FT[0.0, runtime.upper_layer_depth, runtime.deep_transition_depth, runtime.total_depth]
    control_spacings = FT[runtime.surface_Δz, runtime.mid_Δz, runtime.deep_Δz, runtime.deep_Δz]

    faces = FT[-runtime.total_depth]
    z = -runtime.total_depth

    local_spacing(depth) = begin
        if depth <= control_depths[2]
            return linear_interpolate(depth, control_depths[1], control_depths[2],
                                      control_spacings[1], control_spacings[2])
        elseif depth <= control_depths[3]
            return linear_interpolate(depth, control_depths[2], control_depths[3],
                                      control_spacings[2], control_spacings[3])
        else
            return control_spacings[4]
        end
    end

    while z < 0
        depth = -z
        Δz = local_spacing(depth)
        z_next = min(z + Δz, FT(0))
        push!(faces, z_next)
        z = z_next
    end

    return faces
end

linear_interpolate(x, x0, x1, y0, y1) = y0 + (y1 - y0) * ((x - x0) / (x1 - x0))

function build_z_faces(runtime::HM98Runtime)
    if runtime.use_stretched_vertical_grid
        return stretched_z_faces(runtime)
    elseif runtime.Δz_uniform > 0
        return uniform_z_faces(runtime)
    else
        return piecewise_z_faces(runtime)
    end
end

minimum_spacing(z_faces) = (minimum(diff(z_faces)), maximum(diff(z_faces)))
random_buoyancy_noise() = FT(1e-6) * FT(randn())

"""
Weak near-surface perturbation used to break symmetry and trigger convection.
"""
initial_buoyancy_perturbation(z) = random_buoyancy_noise() * exp(z / FT(150))

function maybe_reduce_noise(z, runtime::HM98Runtime)
    if z < -runtime.upper_layer_depth
        return FT(0.1) * initial_buoyancy_perturbation(z)
    end

    return initial_buoyancy_perturbation(z)
end

function patterned_surface_buoyancy_flux(x, y, t, p)
    Bflux = p.B_half * (tanh(-FT(2) * (y - p.Ly / FT(2)) / p.Lf) + FT(1))

    if !isnothing(p.cooling_stop_time) && t > p.cooling_stop_time
        return zero(Bflux)
    end

    return Bflux
end

function buoyancy_boundary_conditions(experiment::HM98Experiment, runtime::HM98Runtime)
    if experiment.forcing == :uniform
        return FieldBoundaryConditions(
            south = GradientBoundaryCondition(FT(0)),
            north = GradientBoundaryCondition(FT(0)),
            bottom = GradientBoundaryCondition(FT(0)),
            top = FluxBoundaryCondition(experiment.B_half)
        )
    elseif experiment.forcing == :patterned
        cooling_stop_time = isnothing(experiment.cooling_stop_days) ? nothing : experiment.cooling_stop_days * seconds_per_day
        return FieldBoundaryConditions(
            south = GradientBoundaryCondition(FT(0)),
            north = GradientBoundaryCondition(FT(0)),
            bottom = GradientBoundaryCondition(FT(0)),
            top = FluxBoundaryCondition(patterned_surface_buoyancy_flux,
                                        parameters = (; B_half = experiment.B_half,
                                                       Ly = runtime.Ly,
                                                       Lf = experiment.Lf,
                                                       cooling_stop_time))
        )
    end

    error("Unsupported forcing style $(experiment.forcing)")
end

background_buoyancy(x, y, z, t) = N² * z
u_initial(x, y, z) = FT(0)
v_initial(x, y, z) = FT(1e-4) * FT(randn()) * exp(z / FT(200))
w_initial(x, y, z) = FT(1e-4) * FT(randn()) * exp(z / FT(200))

"""
    build_grid(experiment, runtime)

Construct either a quasi-2D grid (`Nx = 1`) or a fully 3D grid.
"""
function build_grid(experiment::HM98Experiment, runtime::HM98Runtime)
    z_faces = build_z_faces(runtime)
    Nz = length(z_faces) - 1
    Ny = Int(round(runtime.Ly / runtime.Δy_target))

    if experiment.dimensions == 2
        grid = RectilinearGrid(
            CPU(),
            FT,
            size = (1, Ny, Nz),
            x = (0, runtime.Lx),
            y = (0, runtime.Ly),
            z = z_faces,
            topology = (Periodic, Bounded, Bounded)
        )
    else
        Nx = Int(round(runtime.Lx / runtime.Δx_target))
        grid = RectilinearGrid(
            CPU(),
            FT,
            size = (Nx, Ny, Nz),
            x = (0, runtime.Lx),
            y = (0, runtime.Ly),
            z = z_faces,
            topology = (Periodic, Bounded, Bounded)
        )
    end

    return grid, z_faces
end

"""
Build the Oceananigans model for the selected HM98 experiment.
"""
function build_model(experiment::HM98Experiment, runtime::HM98Runtime, grid)
    closure = (
        HorizontalScalarDiffusivity(FT; ν = runtime.νh, κ = runtime.κh),
        VerticalScalarDiffusivity(FT; ν = runtime.νz, κ = runtime.κz)
    )

    free_slip = FieldBoundaryConditions(
        bottom = GradientBoundaryCondition(FT(0)),
        top = GradientBoundaryCondition(FT(0))
    )

    side_walls = FieldBoundaryConditions(
        south = GradientBoundaryCondition(FT(0)),
        north = GradientBoundaryCondition(FT(0)),
        bottom = GradientBoundaryCondition(FT(0)),
        top = GradientBoundaryCondition(FT(0))
    )

    model = NonhydrostaticModel(grid;
        advection = UpwindBiased(FT; order = 3),
        tracers = :b,
        buoyancy = BuoyancyTracer(),
        coriolis = FPlane(FT; f = f),
        closure,
        boundary_conditions = (
            u = side_walls,
            v = free_slip,
            b = buoyancy_boundary_conditions(experiment, runtime)
        ),
        background_fields = (b = BackgroundField(background_buoyancy),)
    )

    b_initial(x, y, z) = maybe_reduce_noise(z, runtime)
    set!(model, u = u_initial, v = v_initial, w = w_initial, b = b_initial)

    return model
end

function progress(sim)
    u, v, w = sim.model.velocities
    b = sim.model.tracers.b
    @info "iteration=$(iteration(sim)) time_days=$(round(time(sim) / seconds_per_day, digits=2)) Δt=$(round(sim.Δt, digits=2)) max|u|=$(maximum(abs, interior(u))) max|v|=$(maximum(abs, interior(v))) max|w|=$(maximum(abs, interior(w))) max|b'|=$(maximum(abs, interior(b)))"
end

"""
    normalized_pv_field(model)

Construct the full Ertel PV diagnostic normalized by the initial-condition PV `f N²`.
The calculation is based on the total buoyancy field `b + N² z`.
"""
function normalized_pv_field(model)
    u, v, w = model.velocities
    b = model.tracers.b
    B = b + model.background_fields.tracers.b

    ωx = ∂y(w) - ∂z(v)
    ωy = ∂z(u) - ∂x(w)
    ωz = ∂x(v) - ∂y(u)

    q = @at (Center, Center, Center) ((ωx * ∂x(B)) + (ωy * ∂y(B)) + ((f + ωz) * ∂z(B))) / Q₀
    return q
end

"""
Write the 2D experiments as a single quasi-2D field file.
"""
function add_output_writer_2d!(simulation, project_dir, experiment::HM98Experiment, runtime::HM98Runtime)
    filepath = default_output_path(project_dir, experiment)
    writer = JLD2Writer(simulation.model,
                        (; u = simulation.model.velocities.u,
                           w = simulation.model.velocities.w,
                           b = simulation.model.tracers.b);
                        filename = filepath,
                        schedule = TimeInterval(runtime.snapshot_interval_days * seconds_per_day),
                        overwrite_existing = true)

    simulation.output_writers[:fields] = writer
    write_output!(writer, simulation)

    return (fields = filepath,)
end

"""
Write only the slices needed for experiment 4 movies to keep output sizes manageable.
"""
function add_output_writer_3d!(simulation, project_dir, experiment::HM98Experiment, runtime::HM98Runtime)
    model = simulation.model
    grid = model.grid
    x_mid = Int(cld(size(grid, 1), 2))
    z_nodes = vec(znodes(grid, Center()))
    pv_k = argmin(abs.(z_nodes .+ runtime.pv_slice_depth))
    buoyancy_k = argmin(abs.(z_nodes .+ runtime.buoyancy_slice_depth))
    q = normalized_pv_field(model)

    yz_path = default_yz_output_path(project_dir, experiment)
    xy_path = default_xy_output_path(project_dir, experiment)
    xy_buoyancy_path = default_xy_buoyancy_output_path(project_dir, experiment)

    yz_writer = JLD2Writer(model,
                           (; u = model.velocities.u,
                              w = model.velocities.w,
                              b = model.tracers.b,
                              q = q);
                           filename = yz_path,
                           indices = (x_mid, :, :),
                           schedule = TimeInterval(runtime.snapshot_interval_days * seconds_per_day),
                           overwrite_existing = true)

    xy_writer = JLD2Writer(model,
                           (; q = q);
                           filename = xy_path,
                           indices = (:, :, pv_k),
                           schedule = TimeInterval(runtime.snapshot_interval_days * seconds_per_day),
                           overwrite_existing = true)

    xy_buoyancy_writer = JLD2Writer(model,
                                    (; b = model.tracers.b);
                                    filename = xy_buoyancy_path,
                                    indices = (:, :, buoyancy_k),
                                    schedule = TimeInterval(runtime.snapshot_interval_days * seconds_per_day),
                                    overwrite_existing = true)

    simulation.output_writers[:yz_slice] = yz_writer
    simulation.output_writers[:xy_slice] = xy_writer
    simulation.output_writers[:xy_buoyancy_slice] = xy_buoyancy_writer
    write_output!(yz_writer, simulation)
    write_output!(xy_writer, simulation)
    write_output!(xy_buoyancy_writer, simulation)

    return (yz = yz_path, xy = xy_path, xy_buoyancy = xy_buoyancy_path,
            x_mid = x_mid, pv_k = pv_k, buoyancy_k = buoyancy_k)
end

function add_output_writers!(simulation, project_dir, experiment::HM98Experiment, runtime::HM98Runtime)
    if experiment.dimensions == 2
        return add_output_writer_2d!(simulation, project_dir, experiment, runtime)
    else
        return add_output_writer_3d!(simulation, project_dir, experiment, runtime)
    end
end

function run_case(project_dir, experiment::HM98Experiment, runtime::HM98Runtime)
    output_dir = joinpath(project_dir, "outputs")
    mkpath(output_dir)

    grid, z_faces = build_grid(experiment, runtime)
    model = build_model(experiment, runtime, grid)
    simulation = Simulation(model, Δt = runtime.initial_dt, stop_time = runtime.simulation_days * seconds_per_day)

    simulation.callbacks[:wizard] = Callback(TimeStepWizard(cfl = 0.4, max_change = 1.1, max_Δt = runtime.max_dt),
                                             IterationInterval(10))
    simulation.callbacks[:progress] = Callback(progress, IterationInterval(200))

    output_info = add_output_writers!(simulation, project_dir, experiment, runtime)
    min_Δz, max_Δz = minimum_spacing(z_faces)

    @info "Starting Haine-Marshall 1998 $(experiment.name)"
    @info "architecture = CPU()"
    @info "grid = $(size(grid, 1)) x $(size(grid, 2)) x $(size(grid, 3))"
    @info "domain = Lx=$(runtime.Lx / 1e3) km, Ly=$(runtime.Ly / 1e3) km, Lz=$(round(-minimum(z_faces), digits = 1)) m"
    @info "resolution = Δx≈$(round(runtime.Lx / size(grid, 1), digits = 2)) m, Δy≈$(round(runtime.Ly / size(grid, 2), digits = 2)) m, minΔz=$(round(min_Δz, digits = 2)) m, maxΔz=$(round(max_Δz, digits = 2)) m"
    @info "N = $(N) s^-1, f = $(f) s^-1, B_half = $(experiment.B_half) m^2 s^-3, Lf = $(experiment.Lf / 1e3) km"
    @info "νh = $(runtime.νh) m^2 s^-1, κh = $(runtime.κh) m^2 s^-1, νz = $(runtime.νz) m^2 s^-1, κz = $(runtime.κz) m^2 s^-1"
    @info "simulation_days = $(runtime.simulation_days), snapshot_interval_days = $(runtime.snapshot_interval_days)"

    wall_seconds = @elapsed run!(simulation)
    sim_days = time(simulation) / seconds_per_day
    projected_ten_day_hours = sim_days > 0 ? 10 * wall_seconds / sim_days / 3600 : Inf

    @info "Finished simulation at day $(sim_days) after $(iteration(simulation)) iterations."
    @info "Wall-clock runtime = $(round(wall_seconds, digits = 2)) s"
    @info "Projected 10-day wall-clock runtime (linear estimate) = $(round(projected_ten_day_hours, digits = 2)) h"

    if experiment.dimensions == 2
        @info "Wrote fields to $(output_info.fields)"
    else
        @info "Wrote yz slices to $(output_info.yz)"
        @info "Wrote xy PV slices to $(output_info.xy)"
        @info "Wrote xy buoyancy slices to $(output_info.xy_buoyancy)"
    end

    return (; wall_seconds, projected_ten_day_hours, output_info)
end

collapse_quasi_2d(data) = ndims(data) == 2 ? data :
                          ndims(data) == 3 && size(data, 1) == 1 ? dropdims(data; dims = 1) :
                          ndims(data) == 3 && size(data, 1) == 2 ? @views(0.5 .* (data[1, :, :] .+ data[2, :, :])) :
                          ndims(data) == 3 && size(data, 3) == 1 ? dropdims(data; dims = 3) :
                          error("Expected a quasi-2D field with singleton x-dimension, got size $(size(data)).")

centered_field_array(field) = collapse_quasi_2d(Array(interior(field)))

function centered_vertical_velocity_array(field)
    data = centered_field_array(field)
    return size(data, 2) > 1 ? @views(0.5 .* (data[:, 1:end-1] .+ data[:, 2:end])) : data
end

function yz_plot_geometry(timeseries, runtime::HM98Runtime)
    grid = timeseries.grid
    y_nodes = vec(ynodes(grid, Center()))
    z_nodes = vec(znodes(grid, Center()))
    z_plot_mask = z_nodes .>= -runtime.plot_depth
    z_plot_nodes = z_nodes[z_plot_mask]
    return y_nodes, z_nodes, z_plot_mask, z_plot_nodes
end

function xy_plot_geometry(timeseries)
    grid = timeseries.grid
    x_nodes = vec(xnodes(grid, Center()))
    y_nodes = vec(ynodes(grid, Center()))
    return x_nodes, y_nodes
end

function total_buoyancy_snapshot(b_field, z_nodes, z_plot_mask)
    b′ = centered_field_array(b_field)
    b_total = similar(b′)

    for k in axes(b_total, 2)
        b_total[:, k] .= b′[:, k] .+ N² * z_nodes[k]
    end

    return b_total[:, z_plot_mask]
end

function total_horizontal_buoyancy_snapshot(b_field, z_nodes)
    b′ = centered_field_array(b_field)
    z_center = ndims(b′) == 2 ? z_nodes[1] : error("Expected a 2D horizontal slice, got size $(size(b′)).")
    return b′ .+ N² * z_center
end

field_limits(snapshots; symmetric = false) = symmetric ?
    let values = reduce(vcat, vec.(snapshots)); bound = maximum(abs, values); (-bound, bound) end :
    extrema(reduce(vcat, vec.(snapshots)))

"""
Compute a centered finite-difference derivative along the first dimension of a
2D `y-z` array, falling back to one-sided differences at the boundaries.
"""
function first_derivative_y(field, y_nodes)
    derivative = similar(field)

    for j in axes(field, 1)
        if j == firstindex(y_nodes)
            Δy = y_nodes[j + 1] - y_nodes[j]
            @views derivative[j, :] .= (field[j + 1, :] .- field[j, :]) ./ Δy
        elseif j == lastindex(y_nodes)
            Δy = y_nodes[j] - y_nodes[j - 1]
            @views derivative[j, :] .= (field[j, :] .- field[j - 1, :]) ./ Δy
        else
            Δy = y_nodes[j + 1] - y_nodes[j - 1]
            @views derivative[j, :] .= (field[j + 1, :] .- field[j - 1, :]) ./ Δy
        end
    end

    return derivative
end

"""
Compute a centered finite-difference derivative along the second dimension of a
2D `y-z` array, falling back to one-sided differences at the boundaries.
"""
function first_derivative_z(field, z_nodes)
    derivative = similar(field)

    for k in axes(field, 2)
        if k == firstindex(z_nodes)
            Δz = z_nodes[k + 1] - z_nodes[k]
            @views derivative[:, k] .= (field[:, k + 1] .- field[:, k]) ./ Δz
        elseif k == lastindex(z_nodes)
            Δz = z_nodes[k] - z_nodes[k - 1]
            @views derivative[:, k] .= (field[:, k] .- field[:, k - 1]) ./ Δz
        else
            Δz = z_nodes[k + 1] - z_nodes[k - 1]
            @views derivative[:, k] .= (field[:, k + 1] .- field[:, k - 1]) ./ Δz
        end
    end

    return derivative
end

function plot_buoyancy_case(project_dir, experiment::HM98Experiment, runtime::HM98Runtime)
    output_path = default_output_path(project_dir, experiment)
    movie_path = default_movie_path(project_dir, experiment)

    b_timeseries = FieldTimeSeries(output_path, "b")
    y_nodes, z_nodes, z_plot_mask, z_plot_nodes = yz_plot_geometry(b_timeseries, runtime)
    snapshot_times_days = b_timeseries.times ./ seconds_per_day
    buoyancy_snapshots = [total_buoyancy_snapshot(b_timeseries[n], z_nodes, z_plot_mask) for n in eachindex(snapshot_times_days)]

    snapshot_obs = Observable(first(buoyancy_snapshots))
    title_obs = Observable("HM98 $(experiment.name) buoyancy, day $(round(first(snapshot_times_days), digits = 2))")

    fig = Figure(size = (1100, 700))
    ax = Axis(fig[1, 1], xlabel = "y (km)", ylabel = "z (m)", title = title_obs)
    hm = heatmap!(ax, y_nodes ./ 1e3, z_plot_nodes, snapshot_obs;
                  colorrange = field_limits(buoyancy_snapshots),
                  colormap = :thermal)
    Colorbar(fig[1, 2], hm, label = "total buoyancy (m s^-2)")

    record(fig, movie_path, eachindex(buoyancy_snapshots); framerate = 10) do n
        snapshot_obs[] = buoyancy_snapshots[n]
        title_obs[] = "HM98 $(experiment.name) buoyancy, day $(round(snapshot_times_days[n], digits = 2))"
    end

    @info "Wrote movie to $(movie_path)"
    return movie_path
end

function plot_diagnostics_case_2d(project_dir, experiment::HM98Experiment, runtime::HM98Runtime)
    output_path = default_output_path(project_dir, experiment)
    movie_path = default_movie_path(project_dir, experiment)

    u_timeseries = FieldTimeSeries(output_path, "u")
    w_timeseries = FieldTimeSeries(output_path, "w")
    b_timeseries = FieldTimeSeries(output_path, "b")

    y_nodes, z_nodes, z_plot_mask, z_plot_nodes = yz_plot_geometry(b_timeseries, runtime)
    snapshot_times_days = b_timeseries.times ./ seconds_per_day

    u_snapshots = [centered_field_array(u_timeseries[n])[:, z_plot_mask] for n in eachindex(snapshot_times_days)]
    w_snapshots = [centered_vertical_velocity_array(w_timeseries[n])[:, z_plot_mask] for n in eachindex(snapshot_times_days)]
    buoyancy_snapshots = [total_buoyancy_snapshot(b_timeseries[n], z_nodes, z_plot_mask) for n in eachindex(snapshot_times_days)]
    pv_snapshots = [begin
        u = u_snapshots[n]
        B = buoyancy_snapshots[n]
        ∂u_∂y = first_derivative_y(u, y_nodes)
        ∂u_∂z = first_derivative_z(u, z_plot_nodes)
        ∂b_∂y = first_derivative_y(B, y_nodes)
        ∂b_∂z = first_derivative_z(B, z_plot_nodes)
        @. ((f - ∂u_∂y) * ∂b_∂z + ∂u_∂z * ∂b_∂y) / Q₀
    end for n in eachindex(snapshot_times_days)]

    render_yz_diagnostics(movie_path, experiment, snapshot_times_days, y_nodes, z_plot_nodes,
                          u_snapshots, pv_snapshots, buoyancy_snapshots, w_snapshots;
                          pv_limits = isnothing(experiment.pv_limits) ? field_limits(pv_snapshots; symmetric = true) : experiment.pv_limits,
                          contours_on_all_panels = experiment.contours_on_all_panels)
end

function render_yz_diagnostics(movie_path, experiment, snapshot_times_days, y_nodes, z_plot_nodes,
                               u_snapshots, pv_snapshots, buoyancy_snapshots, w_snapshots;
                               pv_limits, contours_on_all_panels)
    u_limits = field_limits(u_snapshots; symmetric = true)
    w_limits = field_limits(w_snapshots; symmetric = true)
    buoyancy_limits = field_limits(buoyancy_snapshots)
    contour_levels = range(buoyancy_limits[1], buoyancy_limits[2]; length = 12)

    u_obs = Observable(first(u_snapshots))
    pv_obs = Observable(first(pv_snapshots))
    buoyancy_obs = Observable(first(buoyancy_snapshots))
    w_obs = Observable(first(w_snapshots))
    title_obs = Observable("HM98 $(experiment.name) diagnostics, day $(round(first(snapshot_times_days), digits = 2))")

    fig = Figure(size = (1700, 1050))
    ax_u = Axis(fig[1, 1], xlabel = "y (km)", ylabel = "z (m)", title = title_obs)
    ax_pv = Axis(fig[1, 3], xlabel = "y (km)", ylabel = "z (m)", title = "PV / (f N²)")
    ax_b = Axis(fig[2, 1], xlabel = "y (km)", ylabel = "z (m)", title = "Buoyancy")
    ax_w = Axis(fig[2, 3], xlabel = "y (km)", ylabel = "z (m)", title = "w")

    u_hm = heatmap!(ax_u, y_nodes ./ 1e3, z_plot_nodes, u_obs; colorrange = u_limits, colormap = :balance)
    pv_hm = heatmap!(ax_pv, y_nodes ./ 1e3, z_plot_nodes, pv_obs; colorrange = pv_limits, colormap = :balance)
    b_hm = heatmap!(ax_b, y_nodes ./ 1e3, z_plot_nodes, buoyancy_obs; colorrange = buoyancy_limits, colormap = :thermal)
    w_hm = heatmap!(ax_w, y_nodes ./ 1e3, z_plot_nodes, w_obs; colorrange = w_limits, colormap = :balance)

    contour!(ax_u, y_nodes ./ 1e3, z_plot_nodes, buoyancy_obs; levels = contour_levels, color = (:black, 0.7), linewidth = 1)

    if contours_on_all_panels
        contour!(ax_pv, y_nodes ./ 1e3, z_plot_nodes, buoyancy_obs; levels = contour_levels, color = (:black, 0.7), linewidth = 1)
        contour!(ax_b, y_nodes ./ 1e3, z_plot_nodes, buoyancy_obs; levels = contour_levels, color = (:black, 0.7), linewidth = 1)
        contour!(ax_w, y_nodes ./ 1e3, z_plot_nodes, buoyancy_obs; levels = contour_levels, color = (:black, 0.7), linewidth = 1)
    end

    Colorbar(fig[1, 2], u_hm, label = "u (m s^-1)")
    Colorbar(fig[1, 4], pv_hm, label = "PV / (f N²)")
    Colorbar(fig[2, 2], b_hm, label = "buoyancy (m s^-2)")
    Colorbar(fig[2, 4], w_hm, label = "w (m s^-1)")

    record(fig, movie_path, eachindex(snapshot_times_days); framerate = 10) do n
        u_obs[] = u_snapshots[n]
        pv_obs[] = pv_snapshots[n]
        buoyancy_obs[] = buoyancy_snapshots[n]
        w_obs[] = w_snapshots[n]
        title_obs[] = "HM98 $(experiment.name) diagnostics, day $(round(snapshot_times_days[n], digits = 2))"
    end

    @info "Wrote movie to $(movie_path)"
    return movie_path
end

function plot_diagnostics_case_3d(project_dir, experiment::HM98Experiment, runtime::HM98Runtime)
    yz_path = default_yz_output_path(project_dir, experiment)
    xy_path = default_xy_output_path(project_dir, experiment)
    xy_buoyancy_path = default_xy_buoyancy_output_path(project_dir, experiment)
    movie_path = default_movie_path(project_dir, experiment)
    xy_movie_path = default_xy_movie_path(project_dir, experiment)
    xy_buoyancy_movie_path = default_xy_buoyancy_movie_path(project_dir, experiment)

    u_timeseries = FieldTimeSeries(yz_path, "u")
    w_timeseries = FieldTimeSeries(yz_path, "w")
    b_timeseries = FieldTimeSeries(yz_path, "b")
    q_timeseries = FieldTimeSeries(yz_path, "q")
    qxy_timeseries = FieldTimeSeries(xy_path, "q")
    bxy_timeseries = FieldTimeSeries(xy_buoyancy_path, "b")

    y_nodes, z_nodes, z_plot_mask, z_plot_nodes = yz_plot_geometry(b_timeseries, runtime)
    x_nodes, y_xy_nodes = xy_plot_geometry(qxy_timeseries)
    snapshot_times_days = b_timeseries.times ./ seconds_per_day

    u_snapshots = [centered_field_array(u_timeseries[n])[:, z_plot_mask] for n in eachindex(snapshot_times_days)]
    w_snapshots = [centered_vertical_velocity_array(w_timeseries[n])[:, z_plot_mask] for n in eachindex(snapshot_times_days)]
    buoyancy_snapshots = [total_buoyancy_snapshot(b_timeseries[n], z_nodes, z_plot_mask) for n in eachindex(snapshot_times_days)]
    pv_snapshots = [centered_field_array(q_timeseries[n])[:, z_plot_mask] for n in eachindex(snapshot_times_days)]

    render_yz_diagnostics(movie_path, experiment, snapshot_times_days, y_nodes, z_plot_nodes,
                          u_snapshots, pv_snapshots, buoyancy_snapshots, w_snapshots;
                          pv_limits = isnothing(experiment.pv_limits) ? field_limits(pv_snapshots; symmetric = true) : experiment.pv_limits,
                          contours_on_all_panels = experiment.contours_on_all_panels)

    qxy_snapshots = [centered_field_array(qxy_timeseries[n]) for n in eachindex(snapshot_times_days)]
    qxy_obs = Observable(first(qxy_snapshots))
    title_obs = Observable("HM98 $(experiment.name) PV, z ≈ -$(round(runtime.pv_slice_depth, digits = 0)) m, day $(round(first(snapshot_times_days), digits = 2))")
    fig = Figure(size = (1400, 840))
    ax = Axis(fig[1, 1], xlabel = "x (km)", ylabel = "y (km)", title = title_obs, aspect = DataAspect())
    hm = heatmap!(ax, x_nodes ./ 1e3, y_xy_nodes ./ 1e3, qxy_obs;
                  colorrange = isnothing(experiment.pv_limits) ? field_limits(qxy_snapshots; symmetric = true) : experiment.pv_limits,
                  colormap = :balance)
    Colorbar(fig[1, 2], hm, label = "PV / (f N²)")

    record(fig, xy_movie_path, eachindex(snapshot_times_days); framerate = 10) do n
        qxy_obs[] = qxy_snapshots[n]
        title_obs[] = "HM98 $(experiment.name) PV, z ≈ -$(round(runtime.pv_slice_depth, digits = 0)) m, day $(round(snapshot_times_days[n], digits = 2))"
    end

    @info "Wrote movie to $(xy_movie_path)"

    bxy_z_nodes = vec(znodes(bxy_timeseries.grid, Center()))
    bxy_snapshots = [total_horizontal_buoyancy_snapshot(bxy_timeseries[n], bxy_z_nodes) for n in eachindex(snapshot_times_days)]
    bxy_limits = field_limits(bxy_snapshots)
    bxy_obs = Observable(first(bxy_snapshots))
    bxy_title_obs = Observable("HM98 $(experiment.name) buoyancy, z ≈ -$(round(runtime.buoyancy_slice_depth, digits = 0)) m, day $(round(first(snapshot_times_days), digits = 2))")
    bxy_fig = Figure(size = (1400, 840))
    bxy_ax = Axis(bxy_fig[1, 1], xlabel = "x (km)", ylabel = "y (km)", title = bxy_title_obs, aspect = DataAspect())
    bxy_hm = heatmap!(bxy_ax, x_nodes ./ 1e3, y_xy_nodes ./ 1e3, bxy_obs; colorrange = bxy_limits, colormap = :thermal)
    Colorbar(bxy_fig[1, 2], bxy_hm, label = "buoyancy (m s^-2)")

    record(bxy_fig, xy_buoyancy_movie_path, eachindex(snapshot_times_days); framerate = 10) do n
        bxy_obs[] = bxy_snapshots[n]
        bxy_title_obs[] = "HM98 $(experiment.name) buoyancy, z ≈ -$(round(runtime.buoyancy_slice_depth, digits = 0)) m, day $(round(snapshot_times_days[n], digits = 2))"
    end

    @info "Wrote movie to $(xy_buoyancy_movie_path)"
    return (yz_movie = movie_path, xy_movie = xy_movie_path, xy_buoyancy_movie = xy_buoyancy_movie_path)
end

function plot_case(project_dir, experiment::HM98Experiment, runtime::HM98Runtime)
    if experiment.plot_kind == :buoyancy
        return plot_buoyancy_case(project_dir, experiment, runtime)
    elseif experiment.dimensions == 2
        return plot_diagnostics_case_2d(project_dir, experiment, runtime)
    else
        return plot_diagnostics_case_3d(project_dir, experiment, runtime)
    end
end

end
