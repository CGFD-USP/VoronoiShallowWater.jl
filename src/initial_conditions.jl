function get_constants(d::Dict, m::AbstractVoronoiMesh)

    di = d["initial_condition"]

    init = if di isa String
        lowercase(di::String)
    else#if isa Dict
        lowercase(di["name"]::String)
    end

    return _get_constants(d, m, Val(Symbol(init)))
end

function create_initial_condition(d::Dict, m::AbstractVoronoiMesh)

    di = d["initial_condition"]

    init = if di isa String
        lowercase(di::String)
    else#if isa Dict
        lowercase(di["name"]::String)
    end

    return _create_initial_condition(d, m, Val(Symbol(init)))
end

include("zonal_geostrophic_balance.jl")

