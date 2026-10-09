abstract type SWModel<:Function end

function SWModel(d::Dict)

    !haskey(d, "grid_file") && @error "No `grid_file` provided in the simulation configuration. Aborting simulation."

    m = VoronoiMesh(d["grid_file"]::String)

    if haskey(d, "grid_scale_factor")
        factor = Float64(d["grid_scale_factor"]::Number)::Float64
        scale!(m, factor)
    end

    model = if haskey(d, "model")
        lowercase(d["model"]::String)
    else
        @warn "No `model` specified in the simulation configuration. Using the \"TRiSK\" model."
        "trisk"
    end

    constants = get_constants(d, m) # constants given by the initial condition.

    if model == "trisk"
        return TriskSWModel(m, constants, d)
    elseif model == "peixoto"
        return PeixotoSWModel(m, constants, d)
    elseif model == "mpas"
        return MPASSWModel(m, constants, d)
    elseif model == "consistent"
        return ConsistentMPASSWModel(m, constants, d)
    elseif model == "wconsistent"
        return WeightedConsistentMPASSWModel(m, constants, d)
    elseif model == "lsq2"
        return LSq2MPASSWModel(m, constants, d)
    elseif model == "wlsq2"
        return WeightedLSq2MPASSWModel(m, constants, d)
    end
end

include("trisk_model_struct.jl")

include("peixoto_model_struct.jl")

include("mpas_model_struct.jl")

include("consistent_model_struct.jl")

include("wconsistent_model_struct.jl")

include("lsq2_model_struct.jl")

include("wlsq2_model_struct.jl")

