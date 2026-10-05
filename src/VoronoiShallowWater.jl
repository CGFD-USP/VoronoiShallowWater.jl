module VoronoiShallowWater

using TOML
using Zeros
using LinearAlgebra
using TensorsLite
using Polyester
using NCDatasets, ReadVTK
using VoronoiMeshes, VoronoiOperators

using VoronoiOperators: mytmap!

export SWModel, TriskSWModel, PeixotoSWModel
export create_initial_condition
export RK4
export SWSimulation
export run_simulation

const Rₑ = 6371220.0
const Ωₑ = 7.292e-5
const g = 9.80616

include("utils.jl")

include("rhs_computation.jl")

include("model_structs.jl")

include("initial_conditions.jl")

include("time_integration_methods.jl")

include("callbacks.jl")

struct SWSimulation{SW<:SWModel, TI, Tinit, TF, TFields, TCall<:Tuple}
    model::SW
    time_integration::TI
    initial_condition::Tinit
    output_callback::WriteOutputCallBack{TF, TFields}
    callbacks::TCall
    dt::Float64
    run_duration::Float64
    current_dt::Base.RefValue{Float64}
    current_time::Base.RefValue{Float64}
    current_iteration::Base.RefValue{Int}
end

SWSimulation(conf_file::AbstractString) = SWSimulation(TOML.parsefile(conf_file))

function SWSimulation(d::Dict)
    model = SWModel(d)
    initial_condition = create_initial_condition(d, model.mesh)
    time_integration = RK4(initial_condition)
    output_callback = WriteOutputCallBack(model.mesh, d["output"]::Dict{String, Any})
    callbacks = create_callbacks(model.mesh, d)
    dt = create_dt(model, d)
    current_dt = Ref(dt)
    run_duration = parse_seconds(d["duration"]::String)
    current_time = Ref(0.0)
    current_iteration = Ref(1)

    return SWSimulation(
        model,
        time_integration,
        initial_condition,
        output_callback,
        callbacks,
        dt,
        run_duration,
        current_dt,
        current_time,
        current_iteration
    )
end

function create_dt(model::SWModel, d::Dict)
    dt = d["time_step"]
    if dt isa Number
        return Float64(dt)::Float64
    elseif dt isa Dict
        cfl = Float64(dt["CFL"]::Number)::Float64
        wave_speed = sqrt(model.constants.g*model.constants.H)
        u_speed = Float64(dt["maximum_velocity"]::Number)::Float64
        Δx = minimum(model.mesh.edges.lengthDual)

        return cfl * (Δx / (wave_speed + u_speed))
    end
end

include("run_simulation.jl")

run_simulation(config_file::String) = run_simulation(SWSimulation(config_file))

end # module VoronoiShallowWater
