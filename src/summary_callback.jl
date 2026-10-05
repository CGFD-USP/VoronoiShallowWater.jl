struct SummaryCallBack <: AbstractCallBack
    n_steps::Int
    timing::Base.RefValue{Float64}
    min_max::Bool
    energy_terms::Bool
end

function SummaryCallBack(::AbstractVoronoiMesh, d::Dict)
    n_steps = if haskey(d, "steps")
        d["steps"]::Int
    else
        10
    end

    max_min = if haskey(d, "min_max")
        d["min_max"]::Bool
    else
        true
    end

    energy_terms = if haskey(d, "energy_terms")
        d["energy_terms"]::Bool
    else
        false
    end

    return SummaryCallBack(n_steps, Ref(0.0), max_min, energy_terms)
end


function (cb::SummaryCallBack)(yⁿ⁺¹, yⁿ, model, current_time::Float64, current_iteration::Int)
    if current_iteration == 1

        str = """

        Current iteration: $current_iteration.
        Current simulation time: $(format_seconds(current_time)).
        """ 
        if cb.min_max
            h_min_max = tminmax(yⁿ⁺¹[1])
            u_min_max = tminmax(yⁿ⁺¹[2])
            str *= """

            h minimun and maximum values: $h_min_max
            u minimum and maximum values: $u_min_max
            """
        end

        if cb.energy_terms
            At = surface_area(model.mesh)
            pe = model.constants.g*(surface_integral(model.mesh, yⁿ⁺¹[1]) / At)
            ke = surface_integral(model.mesh, square, yⁿ⁺¹[2]) / At
            str *= """

            Mean potential energy ((1/Aₜ)*∫(g*h)): $(pe)
            Mean kinetic energy ((1/Aₜ)*∫(uₑ²)): $(ke)
            """
        end

        @info "$str"

        cb.timing[] = time()
    else
        #Time to print summary
        if rem(current_iteration, cb.n_steps) == 0

            time_per_timestep = (time() - cb.timing[]) / cb.n_steps

            str = """

            Current iteration: $current_iteration.
            Current simulation time: $(format_seconds(current_time)).
            Time step timing: ~$(time_per_timestep) seconds per timestep.
            """ 
            if cb.min_max
                h_min_max = tminmax(yⁿ⁺¹[1])
                u_min_max = tminmax(yⁿ⁺¹[2])
                str *= """

                h minimun and maximum values: $h_min_max
                u minimum and maximum values: $u_min_max
                """
            end

            if cb.energy_terms
                At = surface_area(model.mesh)
                pe = model.constants.g*(surface_integral(model.mesh, yⁿ⁺¹[1]) / At)
                ke = surface_integral(model.mesh, square, yⁿ⁺¹[2]) / At
                str *= """

                Mean potential energy ((1/Aₜ)*∫(g*h)): $(pe)
                Mean kinetic energy ((1/Aₜ)*∫(uₑ²)): $(ke)
                """
            end

            @info "$str"

            cb.timing[] = time()
        end
    end
end
