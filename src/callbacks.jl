abstract type AbstractCallBack end

include("write_output_callback.jl")

function create_callbacks(m::AbstractVoronoiMesh, d::Dict)
    cbs = ()

    if haskey(d, "summary")

        ds = d["summary"]::Dict{String, Any}

        turned_on = if haskey(ds, "enabled")
            ds["enabled"]::Bool
        else
            true
        end

        if turned_on
            cbs = (cbs..., SummaryCallBack(m, ds))
        end
    end

    return cbs
end

include("summary_callback.jl")
