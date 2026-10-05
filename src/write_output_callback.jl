abstract type OutputField end

struct HOutput <: OutputField
end

HOutput(::AbstractVoronoiMesh) = HOutput()

struct UOutput <: OutputField
end

UOutput(::AbstractVoronoiMesh) = UOutput()

struct PVOutput <: OutputField
end

PVOutput(::AbstractVoronoiMesh) = PVOutput()

struct KineticEnergyCellOutput <: OutputField
end

KineticEnergyCellOutput(::AbstractVoronoiMesh) = KineticEnergyCellOutput()

struct VelocityCellOutput{TA, T} <: OutputField
    cell_velocity::TA
    cell_velocity_reconstuction::T
end

VelocityCellOutput(m::AbstractVoronoiMesh) =
    VelocityCellOutput(
        VecArray(x = zeros(m.cells.n), y = zeros(m.cells.n)),
        CellVelocityReconstructionLSq2(m)
    )

struct WriteOutputCallBack{T, TFields} <: AbstractCallBack
    file_name::String
    time::Float64 
    precision::Type{T}
    next_time::Base.RefValue{Float64}
    fields::TFields
end

function WriteOutputCallBack(m::AbstractVoronoiMesh, d::Dict)

    if !haskey(d, "file_name")
        @error "A `file_name` string in the `[output]` section must be specified in the configuration file. Aborting simulation" 
    end

    file_name = d["file_name"]::String

    if !haskey(d, "time")
        @error "A `time` value (string or number)  in the `[output]` section must be specified in the configuration file. Aborting simulation" 
    end

    d_time = d["time"]

    time = d_time isa Number ? Float64(d_time)::Float64 : parse_seconds(d_time)

    prec = if haskey(d, "precision")
        p_s = lowercase(d["precision"]::String)
        if (p_s == "double" || p_s == "f64" || p_s == "float64")
            Float64
        elseif (p_s == "single" || p_s == "f32" || p_s == "float32")
            Float32
        elseif (p_s == "half" || p_s == "f16" || p_s == "float16")
            Float16
        else
            @warn "Unreconizided `precision` option '$(d["precision"])' in the `[output]` section. Using `precision = Float64` instead"
            Float64
        end
    else
        Float64
    end

    next_time = Ref(time)

    variables = haskey(d, "variables") ? parse_output_variables(d["variables"]::Vector{String}, m) : (HOutput(m), UOutput(m))
    
    return WriteOutputCallBack(file_name, time, prec, next_time, variables)
end

function parse_output_variables(var::Vector{String}, m::AbstractVoronoiMesh)

    t = "h" in var ? (HOutput(m),) :  ()

    "u" in var && (t = (t..., UOutput(m)))

    "velocity" in var && (t = (t..., VelocityCellOutput(m)))

    "pv" in var && (t = (t..., PVOutput(m)))

    "ke" in var && (t = (t..., KineticEnergyCellOutput(m)))

    return t::Tuple
end

function (cb::WriteOutputCallBack{T})(yⁿ⁺¹, yⁿ, model, current_time, current_iteration) where {T}

    if current_iteration == 1

        create_output_file(cb, yⁿ⁺¹, yⁿ, model)

        cb.next_time[] = cb.time

    else

        if (current_time ≈ cb.next_time[])

            cb.next_time[] += cb.time

            NCDataset(cb.file_name, "a") do ds

                tdata = ds["time"]
                itime = length(tdata) + 1
                tdata[itime] =  current_time

                for var in cb.fields
                    write_variable_netcdf(ds, var, yⁿ⁺¹, model, itime)
                end

            end

        end
        
    end

    return nothing
end

function create_output_file(cb::WriteOutputCallBack{TF}, yⁿ⁺¹, yⁿ, model) where {TF}


    NCDataset(cb.file_name, "c"; format=:netcdf5_64bit_data) do ds

        mesh = model.mesh

        ds.dim["nCells"] = mesh.cells.n 
        ds.dim["nVertices"] = mesh.vertices.n 
        ds.dim["nEdges"] = mesh.edges.n 
        ds.dim["Time"] = Inf

        ds.attrib["configuration"] = toml_string(model.configuration)

        T = defVar(ds, "time", Float64, ("Time",); attrib = ["units" => "s", "long_name" => "Simulation time"])

        T[1] = 0.0

        H = defVar(ds, "H", Float64, (); attrib = ["units" => "m", "long_name" => "Mean fluid Height"]) 

        H[] = model.constants.H

        g = defVar(ds, "g", Float64, (); attrib = ["units" => "m*s⁻²", "long_name" => "Gravity acceleration"]) 

        g[] = model.constants.g

        if model.constants.f isa ConstArray
            f = defVar(ds, "f", Float64, (); attrib = ["units" => "rad*s⁻¹", "long_name" => "Coriolis parameter"]) 
            f[] = model.constants.f[1]
        end

        for var in cb.fields
            create_variable_netcdf(ds, var, cb)
            write_variable_netcdf(ds, var, yⁿ⁺¹, model, 1)
        end

    end
end

function create_variable_netcdf(ds::NCDataset, ::HOutput, ::WriteOutputCallBack{TF}) where {TF}
    defVar(ds, "h", TF, ("nCells", "Time"); attrib = ["units" => "m", "long_name" => "Fluid height deviation from mean height"])
end

function write_variable_netcdf(ds::NCDataset, ::HOutput, y, model, itime::Integer)
    data = ds["h"]
    data[:, itime] = y[1]
end

function create_variable_netcdf(ds::NCDataset, ::UOutput, ::WriteOutputCallBack{TF}) where {TF}
    defVar(ds, "u", TF, ("nEdges", "Time"); attrib = ["units" => "m*s⁻¹", "long_name" => "Normal velocity component at edge"])
end

function write_variable_netcdf(ds::NCDataset, ::UOutput, y, model, itime::Integer)
    data = ds["u"]
    data[:, itime] = y[2]
end

function create_variable_netcdf(ds::NCDataset, ::VelocityCellOutput, ::WriteOutputCallBack{TF}) where {TF}
    defVar(ds, "zonalVelocity", TF, ("nCells", "Time"); attrib = ["units" => "m*s⁻¹", "long_name" => "Reconstructed zonal velocity at cell center"])

    defVar(ds, "meridionalVelocity", TF, ("nCells", "Time"); attrib = ["units" => "m*s⁻¹", "long_name" => "Reconstructed meridional velocity at cell center"])
end

function write_variable_netcdf(ds::NCDataset, uR::VelocityCellOutput, y, model, itime::Integer)
    uR.cell_velocity_reconstuction(uR.cell_velocity, y[2])

    datax = ds["zonalVelocity"]
    datax[:, itime] = uR.cell_velocity.x

    datay = ds["meridionalVelocity"]
    datay[:, itime] = uR.cell_velocity.y
end

function create_variable_netcdf(ds::NCDataset, ::PVOutput, ::WriteOutputCallBack{TF}) where {TF}
        defVar(ds, "pv", TF, ("nEdges", "Time"); attrib = ["units" => "s⁻¹*m⁻¹", "long_name" => "Potential vorticity at edge"])
end

function write_variable_netcdf(ds::NCDataset, ::PVOutput, y, model, itime::Integer)
    data = ds["pv"]

    compute_potential_vorticity!(model.diagnostics, y[1], y[2],
                                 model.constants,
                                 model.operators)

    data[:, itime] = model.diagnostics.pv_e
end

function create_variable_netcdf(ds::NCDataset, ::KineticEnergyCellOutput, ::WriteOutputCallBack{TF}) where {TF}
        defVar(ds, "ke", TF, ("nCells", "Time"); attrib = ["units" => "m²*s⁻²", "long_name" => "Kinetic energy per unit mass at cell center"])
end

function write_variable_netcdf(ds::NCDataset, ::KineticEnergyCellOutput, y, model, itime::Integer)
    data = ds["ke"]

    k_c = model.diagnostics.k_c

    model.operators.k_reconstruction(k_c, y[2])

    data[:, itime] = k_c
end
